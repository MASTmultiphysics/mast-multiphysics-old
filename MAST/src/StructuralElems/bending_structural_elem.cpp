//
//  bending_structural_elem_base.cpp
//  MAST
//
//  Created by Manav Bhatia on 11/27/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//


// MAST includes
#include "StructuralElems/bending_structural_elem.h"
#include "PropertyCards/element_property_card_2D.h"
#include "Numerics/fem_operator_matrix.h"
#include "ThermalElems/temperature_function.h"
#include "BoundaryConditions/boundary_condition.h"
#include "StructuralElems/structural_element_1D.h"
#include "StructuralElems/structural_element_2D.h"


MAST::BendingStructuralElem::BendingStructuralElem(System& sys,
                                                   const Elem& elem,
                                                   const MAST::ElementPropertyCardBase& p):
MAST::StructuralElementBase(sys, elem, p)
{
    _local_elem.reset(MAST::build_local_elem(*this).release());
    
    _init_fe_and_qrule(_local_elem->local_elem());

    // initialize the bending operator
    MAST::BendingOperatorType bending_model =
    _property.bending_model(elem, _fe->get_fe_type());

    _bending_operator.reset(MAST::build_bending_operator(bending_model, *this).release());
}



bool
MAST::BendingStructuralElem::internal_force (bool request_jacobian,
                                           DenseVector<Real>& f,
                                           DenseMatrix<Real>& jac)
{
    FEMOperatorMatrix Bmat_mem, Bmat_bend, Bmat_vk;
    
    const std::vector<Real>& JxW = _fe->get_JxW();
    const std::vector<Point>& xyz = _fe->get_xyz();
    const unsigned int n_phi = (unsigned int)_fe->get_phi().size();
    const unsigned int n1= this->n_direct_strain_components(), n2=6*n_phi,
    n3 = this->n_von_karman_strain_components();
    DenseMatrix<Real> material_A_mat, material_B_mat, material_D_mat,
    tmp_mat1_n1n2, tmp_mat2_n2n2, tmp_mat3,
    tmp_mat4_n3n2, vk_dwdxi_mat, stress, local_jac;
    DenseVector<Real>  tmp_vec1_n1, tmp_vec2_n1, tmp_vec3_n2,
    tmp_vec4_n3, tmp_vec5_n3, local_f;
    
    tmp_mat1_n1n2.resize(n1, n2); tmp_mat2_n2n2.resize(n2, n2);
    tmp_mat4_n3n2.resize(n3, n2); local_jac.resize(n2, n2);
    vk_dwdxi_mat.resize(n1,n3); stress.resize(2,2); local_f.resize(n2);
    tmp_vec1_n1.resize(n1); tmp_vec2_n1.resize(n1);
    tmp_vec3_n2.resize(n2); tmp_vec4_n3.resize(n3); tmp_vec5_n3.resize(n3);
    
    
    Bmat_mem.reinit(n1, _system.n_vars(), n_phi); // three stress-strain components
    Bmat_bend.reinit(n1, _system.n_vars(), n_phi);
    Bmat_vk.reinit(n3, _system.n_vars(), n_phi); // only dw/dx and dw/dy
    
    bool if_vk = (_property.strain_type() == MAST::VON_KARMAN_STRAIN),
    if_bending = (_property.bending_model(_elem, _fe->get_fe_type()) != MAST::NO_BENDING);
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        // if temperature is specified, the initialize it to the current location
        if (_temperature)
            _temperature->initialize(xyz[qp]);
        
        // get the material matrix
        _property.calculate_matrix(_elem,
                                   MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX,
                                   material_A_mat);
        if (if_bending) {
            _property.calculate_matrix(_elem,
                                       MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_B_MATRIX,
                                       material_B_mat);
            _property.calculate_matrix(_elem,
                                       MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_D_MATRIX,
                                       material_D_mat);
        }
        
        // now calculte the quantity for these matrices
        _internal_force_operation(if_bending, if_vk, n2, qp, JxW,
                                  request_jacobian, local_f, local_jac,
                                  Bmat_mem, Bmat_bend, Bmat_vk,
                                  stress, vk_dwdxi_mat, material_A_mat,
                                  material_B_mat, material_D_mat, tmp_vec1_n1,
                                  tmp_vec2_n1, tmp_vec3_n2, tmp_vec4_n3,
                                  tmp_vec5_n3, tmp_mat1_n1n2, tmp_mat2_n2n2,
                                  tmp_mat3, tmp_mat4_n3n2);
        
    }
    
    
    // now calculate the transverse shear contribution if appropriate for the
    // element
    if (if_bending && _bending_operator->include_transverse_shear_energy())
        _bending_operator->calculate_transverse_shear_force(request_jacobian,
                                                            local_f, local_jac,
                                                            NULL);
    
    
    // now transform to the global coorodinate system
    transform_to_global_system(local_f, tmp_vec3_n2);
    f.add(1., tmp_vec3_n2);
    if (request_jacobian) {
        // for 2D elements
        if (_elem.dim() == 2) {
            // add small values to the diagonal of the theta_z dofs
            for (unsigned int i=0; i<n_phi; i++)
                local_jac(5*n_phi+i, 5*n_phi+i) = 1.0e-6;
        }
        transform_to_global_system(local_jac, tmp_mat2_n2n2);
        jac.add(1., tmp_mat2_n2n2);
    }
    
    return request_jacobian;
}





bool
MAST::BendingStructuralElem::internal_force_sensitivity (bool request_jacobian,
                                                       DenseVector<Real>& f,
                                                       DenseMatrix<Real>& jac)
{
    // this should be true if the function is called
    libmesh_assert(this->sensitivity_params);
    libmesh_assert(!this->sensitivity_params->shape_sensitivity()); // this is not implemented for now
    libmesh_assert(this->sensitivity_params->total_order() == 1); // only first order sensitivity
    
    // check if the material property or the provided exterior
    // values, like temperature, are functions of the sensitivity parameter
    bool calculate = false;
    if (_temperature)
        calculate = calculate || _temperature->depends_on(*(this->sensitivity_params));
    calculate = calculate || _property.depends_on(*(this->sensitivity_params));
    
    // nothing to be calculated if the element does not depend on the
    // sensitivity parameter.
    if (!calculate)
        return false;
    
    FEMOperatorMatrix Bmat_mem, Bmat_bend, Bmat_vk;
    
    const std::vector<Real>& JxW = _fe->get_JxW();
    const std::vector<Point>& xyz = _fe->get_xyz();
    const unsigned int n_phi = (unsigned int)_fe->get_phi().size();
    const unsigned int n1= this->n_direct_strain_components(), n2=6*n_phi,
    n3 = this->n_von_karman_strain_components();
    DenseMatrix<Real> material_A_mat, material_B_mat, material_D_mat,
    material_trans_shear_mat, tmp_mat1_n1n2, tmp_mat2_n2n2, tmp_mat3,
    tmp_mat4_n3n2, vk_dwdxi_mat, stress, local_jac;
    DenseVector<Real>  tmp_vec1_n1, tmp_vec2_n1, tmp_vec3_n2,
    tmp_vec4_n3, tmp_vec5_n3, local_f;
    
    tmp_mat1_n1n2.resize(n1, n2); tmp_mat2_n2n2.resize(n2, n2);
    tmp_mat4_n3n2.resize(n3, n2); local_jac.resize(n2, n2);
    vk_dwdxi_mat.resize(n1,n3); stress.resize(2,2); local_f.resize(n2);
    tmp_vec1_n1.resize(n1); tmp_vec2_n1.resize(n1);
    tmp_vec3_n2.resize(n2); tmp_vec4_n3.resize(n3); tmp_vec5_n3.resize(n3);
    
    
    Bmat_mem.reinit(n1, _system.n_vars(), n_phi); // three stress-strain components
    Bmat_bend.reinit(n1, _system.n_vars(), n_phi);
    Bmat_vk.reinit(n3, _system.n_vars(), n_phi); // only dw/dx and dw/dy
    
    bool if_vk = (_property.strain_type() == MAST::VON_KARMAN_STRAIN),
    if_bending = (_property.bending_model(_elem, _fe->get_fe_type()) != MAST::NO_BENDING);
    
    
    // first calculate the sensitivity due to the parameter
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        // if temperature is specified, the initialize it to the current location
        if (_temperature)
            _temperature->initialize(xyz[qp]);
        
        // get the material matrix
        _property.calculate_matrix_sensitivity
        (_elem,
         MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX,
         material_A_mat, *(this->sensitivity_params));
        if (if_bending) {
            _property.calculate_matrix_sensitivity(_elem,
                                                   MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_B_MATRIX,
                                                   material_B_mat,
                                                   *(this->sensitivity_params));
            
            _property.calculate_matrix_sensitivity(_elem,
                                                   MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_D_MATRIX,
                                                   material_D_mat,
                                                   *(this->sensitivity_params));
        }
        
        // now calculte the quantity for these matrices
        _internal_force_operation(if_bending, if_vk, n2, qp, JxW,
                                  request_jacobian, local_f, local_jac,
                                  Bmat_mem, Bmat_bend, Bmat_vk,
                                  stress, vk_dwdxi_mat, material_A_mat,
                                  material_B_mat, material_D_mat, tmp_vec1_n1,
                                  tmp_vec2_n1, tmp_vec3_n2, tmp_vec4_n3,
                                  tmp_vec5_n3, tmp_mat1_n1n2, tmp_mat2_n2n2,
                                  tmp_mat3, tmp_mat4_n3n2);
    }
    
    // now transform to the global coorodinate system
    transform_to_global_system(local_f, tmp_vec3_n2);
    f.add(1., tmp_vec3_n2);
    if (request_jacobian) {
        transform_to_global_system(local_jac, tmp_mat2_n2n2);
        jac.add(1., tmp_mat2_n2n2);
    }
    
    // now calculate the transverse shear contribution if appropriate for the
    // element
    if (if_bending && _bending_operator->include_transverse_shear_energy())
        _bending_operator->calculate_transverse_shear_force(request_jacobian,
                                                            local_f, local_jac,
                                                            this->sensitivity_params);
    
    
    // next, calculate the sensitivity due to temperature, if it is provided
    if (!_temperature)
        return request_jacobian;
    
    local_f.zero();
    local_jac.zero();
    MAST::SensitivityParameters temp_param;
    temp_param.add_parameter(_temperature, 1);
    
    // first calculate the sensitivity due to the parameter
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        // if temperature is specified, the initialize it to the current location
        if (_temperature)
            _temperature->initialize(xyz[qp]);
        
        // get the material matrix
        _property.calculate_matrix_sensitivity
        (_elem,
         MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX,
         material_A_mat, temp_param);
        if (if_bending) {
            _property.calculate_matrix_sensitivity(_elem,
                                                   MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_B_MATRIX,
                                                   material_B_mat,
                                                   temp_param);
            
            _property.calculate_matrix_sensitivity(_elem,
                                                   MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_D_MATRIX,
                                                   material_D_mat,
                                                   temp_param);
        }
        
        // now calculte the quantity for these matrices
        _internal_force_operation(if_bending, if_vk, n2, qp, JxW,
                                  request_jacobian, local_f, local_jac,
                                  Bmat_mem, Bmat_bend, Bmat_vk,
                                  stress, vk_dwdxi_mat, material_A_mat,
                                  material_B_mat, material_D_mat, tmp_vec1_n1,
                                  tmp_vec2_n1, tmp_vec3_n2, tmp_vec4_n3,
                                  tmp_vec5_n3, tmp_mat1_n1n2, tmp_mat2_n2n2,
                                  tmp_mat3, tmp_mat4_n3n2);
    }
    
    // now transform to the global coorodinate system
    Real dTemp_dparam = _temperature->sensitivity(*(this->sensitivity_params));
    transform_to_global_system(local_f, tmp_vec3_n2);
    f.add(dTemp_dparam, tmp_vec3_n2);
    if (request_jacobian) {
        transform_to_global_system(local_jac, tmp_mat2_n2n2);
        jac.add(dTemp_dparam, tmp_mat2_n2n2);
    }
    
    return request_jacobian;
}




bool
MAST::BendingStructuralElem::prestress_force (bool request_jacobian,
                                            DenseVector<Real>& f,
                                            DenseMatrix<Real>& jac)
{
    if (!_property.if_prestressed())
        return false;
    
    FEMOperatorMatrix Bmat_mem, Bmat_bend, Bmat_vk;
    
    const std::vector<Real>& JxW = _fe->get_JxW();
    const unsigned int n_phi = (unsigned int)_fe->get_phi().size();
    const unsigned int n1= this->n_direct_strain_components(), n2=6*n_phi,
    n3 = this->n_von_karman_strain_components();
    DenseMatrix<Real> tmp_mat2_n2n2, tmp_mat3, vk_dwdxi_mat, local_jac,
    prestress_mat_A;
    DenseVector<Real> tmp_vec2_n1, tmp_vec3_n2, tmp_vec4_n3, tmp_vec5_n3,
    local_f, prestress_vec_A, prestress_vec_B;
    
    tmp_mat2_n2n2.resize(n2, n2); local_jac.resize(n2, n2);
    vk_dwdxi_mat.resize(n1,n3); local_f.resize(n2); tmp_vec2_n1.resize(n1);
    tmp_vec3_n2.resize(n2); tmp_vec4_n3.resize(n3); tmp_vec5_n3.resize(n3);
    
    
    Bmat_mem.reinit(n1, _system.n_vars(), n_phi); // three stress-strain components
    Bmat_bend.reinit(n1, _system.n_vars(), n_phi);
    Bmat_vk.reinit(n3, _system.n_vars(), n_phi); // only dw/dx and dw/dy
    
    bool if_vk = (_property.strain_type() == MAST::VON_KARMAN_STRAIN),
    if_bending = (_property.bending_model(_elem, _fe->get_fe_type()) != MAST::NO_BENDING);
    
    
    // get the element prestress
    _property.prestress_matrix(SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX,
                               prestress_mat_A);
    _property.prestress_vector(SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX,
                               prestress_vec_A);
    
    _property.prestress_vector(SECTION_INTEGRATED_MATERIAL_STIFFNESS_B_MATRIX,
                               prestress_vec_B);
    // transform to the local coordinate system
    
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        this->initialize_direct_strain_operator(qp, Bmat_mem);
        
        // get the bending strain operator if needed
        tmp_vec2_n1.zero(); // used to store vk strain, if applicable
        if (if_bending) {
            _bending_operator->initialize_bending_strain_operator(qp, Bmat_bend);
            
            if (if_vk)  // get the vonKarman strain operator if needed
                this->initialize_von_karman_strain_operator(qp,
                                                            tmp_vec2_n1,
                                                            vk_dwdxi_mat,
                                                            Bmat_vk);
        }
        
        // first handle constant throught the thickness stresses: membrane and vonKarman
        // multiply this with the constant through the thickness strain
        // membrane strain
        Bmat_mem.vector_mult_transpose(tmp_vec3_n2, prestress_vec_A);
        local_f.add(-JxW[qp], tmp_vec3_n2); // epsilon_mem * sigma_0
        
        if (if_bending) {
            if (if_vk) {
                // von Karman strain
                vk_dwdxi_mat.vector_mult_transpose(tmp_vec4_n3, prestress_vec_A);
                Bmat_vk.vector_mult_transpose(tmp_vec3_n2, tmp_vec4_n3);
                local_f.add(-JxW[qp], tmp_vec3_n2); // epsilon_vk * sigma_0
            }
            
            // now coupling with the bending strain
            Bmat_bend.vector_mult_transpose(tmp_vec3_n2, prestress_vec_B);
            local_f.add(-JxW[qp], tmp_vec3_n2); // epsilon_bend * sigma_0
        }
        
        if (request_jacobian) {
            if (if_bending && if_vk) {
                tmp_mat3.resize(2, n2);
                Bmat_vk.left_multiply(tmp_mat3, prestress_mat_A);
                Bmat_vk.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat3);
                local_jac.add(-JxW[qp], tmp_mat2_n2n2);
            }
        }
    }
    
    // now transform to the global coorodinate system
    transform_to_global_system(local_f, tmp_vec3_n2);
    f.add(1., tmp_vec3_n2);
    if (request_jacobian && if_vk) {
        transform_to_global_system(local_jac, tmp_mat2_n2n2);
        jac.add(1., tmp_mat2_n2n2);
    }
    
    // only the nonlinear strain returns a Jacobian for prestressing
    return (request_jacobian && if_vk);
}





bool
MAST::BendingStructuralElem::prestress_force_sensitivity (bool request_jacobian,
                                                        DenseVector<Real>& f,
                                                        DenseMatrix<Real>& jac)
{
    if (!_property.if_prestressed())
        return false;
    
    libmesh_error(); // revisit to include prestress matrix sensitivity
    
    FEMOperatorMatrix Bmat_mem, Bmat_bend, Bmat_vk;
    
    const std::vector<Real>& JxW = _fe->get_JxW();
    const unsigned int n_phi = (unsigned int)_fe->get_phi().size();
    const unsigned int n1= this->n_direct_strain_components(), n2=6*n_phi,
    n3 = this->n_von_karman_strain_components();
    DenseMatrix<Real> tmp_mat2_n2n2, tmp_mat3, vk_dwdxi_mat, local_jac,
    prestress_mat_A;
    DenseVector<Real> tmp_vec2_n1, tmp_vec3_n2, tmp_vec4_n3, tmp_vec5_n3,
    local_f, prestress_vec_A, prestress_vec_B;
    
    tmp_mat2_n2n2.resize(n2, n2); local_jac.resize(n2, n2);
    vk_dwdxi_mat.resize(n1,n3); local_f.resize(n2); tmp_vec2_n1.resize(n1);
    tmp_vec3_n2.resize(n2); tmp_vec4_n3.resize(n3); tmp_vec5_n3.resize(n3);
    
    
    Bmat_mem.reinit(n1, _system.n_vars(), n_phi); // three stress-strain components
    Bmat_bend.reinit(n1, _system.n_vars(), n_phi);
    Bmat_vk.reinit(n3, _system.n_vars(), n_phi); // only dw/dx and dw/dy
    
    bool if_vk = (_property.strain_type() == MAST::VON_KARMAN_STRAIN),
    if_bending = (_property.bending_model(_elem, _fe->get_fe_type()) != MAST::NO_BENDING);
    
    
    // get the element prestress
    _property.prestress_matrix(SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX,
                               prestress_mat_A);
    _property.prestress_vector(SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX,
                               prestress_vec_A);
    
    _property.prestress_vector(SECTION_INTEGRATED_MATERIAL_STIFFNESS_B_MATRIX,
                               prestress_vec_B);
    // transform to the local coordinate system
    
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        this->initialize_direct_strain_operator(qp, Bmat_mem);
        
        // get the bending strain operator if needed
        tmp_vec2_n1.zero(); // used to store vk strain, if applicable
        if (if_bending) {
            _bending_operator->initialize_bending_strain_operator(qp, Bmat_bend);
            
            if (if_vk)  // get the vonKarman strain operator if needed
                this->initialize_von_karman_strain_operator(qp,
                                                            tmp_vec2_n1,
                                                            vk_dwdxi_mat,
                                                            Bmat_vk);
        }
        
        // first handle constant throught the thickness stresses: membrane and vonKarman
        // multiply this with the constant through the thickness strain
        // membrane strain
        Bmat_mem.vector_mult_transpose(tmp_vec3_n2, prestress_vec_A);
        local_f.add(-JxW[qp], tmp_vec3_n2); // epsilon_mem * sigma_0
        
        if (if_bending) {
            if (if_vk) {
                // von Karman strain
                vk_dwdxi_mat.vector_mult_transpose(tmp_vec4_n3, prestress_vec_A);
                Bmat_vk.vector_mult_transpose(tmp_vec3_n2, tmp_vec4_n3);
                local_f.add(-JxW[qp], tmp_vec3_n2); // epsilon_vk * sigma_0
            }
            
            // now coupling with the bending strain
            Bmat_bend.vector_mult_transpose(tmp_vec3_n2, prestress_vec_B);
            local_f.add(-JxW[qp], tmp_vec3_n2); // epsilon_bend * sigma_0
        }
        
        if (request_jacobian) {
            if (if_bending && if_vk) {
                tmp_mat3.resize(2, n2);
                Bmat_vk.left_multiply(tmp_mat3, prestress_mat_A);
                Bmat_vk.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat3);
                local_jac.add(-JxW[qp], tmp_mat2_n2n2);
            }
        }
    }
    
    // now transform to the global coorodinate system
    transform_to_global_system(local_f, tmp_vec3_n2);
    f.add(1., tmp_vec3_n2);
    if (request_jacobian && if_vk) {
        transform_to_global_system(local_jac, tmp_mat2_n2n2);
        jac.add(1., tmp_mat2_n2n2);
    }
    
    // only the nonlinear strain returns a Jacobian for prestressing
    return (request_jacobian && if_vk);
}





bool
MAST::BendingStructuralElem::thermal_force (bool request_jacobian,
                                          DenseVector<Real>& f,
                                          DenseMatrix<Real>& jac)
{
    if (!_temperature) // only if a temperature load is specified
        return false;
    
    libmesh_error(); // to be implemented
    
    FEMOperatorMatrix Bmat;
    
    const std::vector<Real>& JxW = _fe->get_JxW();
    const std::vector<Point>& xyz = _fe->get_xyz();
    const unsigned int n_phi = (unsigned int)JxW.size();
    const unsigned int n1=6, n2=6*n_phi;
    DenseMatrix<Real> material_mat, expansion_mat;
    DenseVector<Real>  phi, temperature, tmp_vec1_n1, tmp_vec2_n1,
    tmp_vec3_n2, local_f;
    
    phi.resize(n_phi); tmp_vec1_n1.resize(n1); tmp_vec2_n1.resize(n1);
    tmp_vec3_n2.resize(n2); temperature.resize(6); local_f.resize(n2);
    
    Bmat.reinit(3, _system.n_vars(), n_phi); // three stress-strain components
    
    Real temperature_value, ref_temperature;
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        this->initialize_direct_strain_operator(qp, Bmat);
        
        // set the temperature vector to the value at this point
        _temperature->initialize(xyz[qp]);
        temperature_value = (*_temperature)();
        ref_temperature   = _temperature->reference();
        
        // this is moved inside the domain since
        _property.calculate_matrix(_elem,
                                   MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX,
                                   material_mat);
        _property.calculate_matrix(_elem,
                                   MAST::SECTION_INTEGRATED_MATERIAL_THERMAL_EXPANSION_MATRIX,
                                   expansion_mat);
        
        // calculate the strain
        expansion_mat.vector_mult(tmp_vec2_n1, tmp_vec1_n1);
        
        // calculate the stress
        material_mat.vector_mult(tmp_vec1_n1, tmp_vec2_n1);
        
        // now calculate the internal force vector
        Bmat.vector_mult_transpose(tmp_vec3_n2, tmp_vec1_n1);
        local_f.add(JxW[qp], tmp_vec3_n2);
    }
    
    // now transform to the global system and add
    transform_to_global_system(local_f, tmp_vec3_n2);
    f.add(1., tmp_vec3_n2);
    
    return false;
}




bool
MAST::BendingStructuralElem::thermal_force_sensitivity (bool request_jacobian,
                                                      DenseVector<Real>& f,
                                                      DenseMatrix<Real>& jac)
{
    if (!_temperature) // only if a temperature load is specified
        return false;
    
    libmesh_error(); // to be implemented
    
    
    return false;
}






void
MAST::BendingStructuralElem::_internal_force_operation
(bool if_bending,
 bool if_vk,
 const unsigned int n2,
 const unsigned int qp,
 const std::vector<Real>& JxW,
 bool request_jacobian,
 DenseVector<Real>& local_f,
 DenseMatrix<Real>& local_jac,
 FEMOperatorMatrix& Bmat_mem,
 FEMOperatorMatrix& Bmat_bend,
 FEMOperatorMatrix& Bmat_vk,
 DenseMatrix<Real>& stress,
 DenseMatrix<Real>& vk_dwdxi_mat,
 DenseMatrix<Real>& material_A_mat,
 DenseMatrix<Real>& material_B_mat,
 DenseMatrix<Real>& material_D_mat,
 DenseVector<Real>& tmp_vec1_n1,
 DenseVector<Real>& tmp_vec2_n1,
 DenseVector<Real>& tmp_vec3_n2,
 DenseVector<Real>& tmp_vec4_2,
 DenseVector<Real>& tmp_vec5_2,
 DenseMatrix<Real>& tmp_mat1_n1n2,
 DenseMatrix<Real>& tmp_mat2_n2n2,
 DenseMatrix<Real>& tmp_mat3,
 DenseMatrix<Real>& tmp_mat4_2n2)
{
    this->initialize_direct_strain_operator(qp, Bmat_mem);
    
    // get the bending strain operator if needed
    tmp_vec2_n1.zero(); // used to store vk strain, if applicable
    if (if_bending) {
        _bending_operator->initialize_bending_strain_operator(qp, Bmat_bend);

        if (if_vk)  // get the vonKarman strain operator if needed
            this->initialize_von_karman_strain_operator(qp,
                                                        tmp_vec2_n1, // epsilon_vk
                                                        vk_dwdxi_mat,
                                                        Bmat_vk);
    }
    
    
    // first handle constant throught the thickness stresses: membrane and vonKarman
    Bmat_mem.vector_mult(tmp_vec1_n1, local_solution);
    tmp_vec2_n1.add(1., tmp_vec1_n1);  // epsilon_mem + epsilon_vk
    
    
    // multiply this with the constant-through-the-thickness strain
    // membrane strain
    material_A_mat.vector_mult(tmp_vec1_n1, tmp_vec2_n1); // stress
    Bmat_mem.vector_mult_transpose(tmp_vec3_n2, tmp_vec1_n1);
    local_f.add(-JxW[qp], tmp_vec3_n2);
    // copy the stress values to a matrix if needed    
    stress(0,0) = tmp_vec1_n1(0); // sigma_xx
    if (_elem.dim() == 2) { // this is not needed for 1D element
        stress(0,1) = tmp_vec1_n1(2); // sigma_xy
        stress(1,0) = tmp_vec1_n1(2); // sigma_yx
        stress(1,1) = tmp_vec1_n1(1); // sigma_yy
    }
    
    if (if_bending) {
        if (if_vk) {
            // von Karman strain
            vk_dwdxi_mat.vector_mult_transpose(tmp_vec4_2, tmp_vec1_n1);
            Bmat_vk.vector_mult_transpose(tmp_vec3_n2, tmp_vec4_2);
            local_f.add(-JxW[qp], tmp_vec3_n2);
        }
        
        // now coupling with the bending strain
        material_B_mat.vector_mult(tmp_vec1_n1, tmp_vec2_n1);
        Bmat_bend.vector_mult_transpose(tmp_vec3_n2, tmp_vec1_n1);
        local_f.add(-JxW[qp], tmp_vec3_n2);
        
        // now bending stress and its coupling
        Bmat_bend.vector_mult(tmp_vec2_n1, local_solution);
        
        // now get its projection onto the constant through thickness
        // and membrane operators
        material_B_mat.vector_mult(tmp_vec1_n1, tmp_vec2_n1);
        Bmat_mem.vector_mult_transpose(tmp_vec3_n2, tmp_vec1_n1);
        local_f.add(-JxW[qp], tmp_vec3_n2);
        
        if (if_vk) {
            // von Karman strain
            vk_dwdxi_mat.vector_mult_transpose(tmp_vec4_2, tmp_vec1_n1);
            Bmat_vk.vector_mult_transpose(tmp_vec3_n2, tmp_vec4_2);
            local_f.add(-JxW[qp], tmp_vec3_n2);
        }
        
        // and the bending operator
        material_D_mat.vector_mult(tmp_vec1_n1, tmp_vec2_n1);
        Bmat_bend.vector_mult_transpose(tmp_vec3_n2, tmp_vec1_n1);
        local_f.add(-JxW[qp], tmp_vec3_n2);
    }
    
    if (request_jacobian) {
        // membrane - membrane
        Bmat_mem.left_multiply(tmp_mat1_n1n2, material_A_mat);
        Bmat_mem.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat1_n1n2);
        local_jac.add(-JxW[qp], tmp_mat2_n2n2);
        
        if (if_bending) {
            if (if_vk) {
                // membrane - vk
                tmp_mat3.resize(vk_dwdxi_mat.m(), n2);
                Bmat_vk.left_multiply(tmp_mat3, vk_dwdxi_mat);
                tmp_mat3.left_multiply(material_A_mat);
                Bmat_mem.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat3);
                local_jac.add(-JxW[qp], tmp_mat2_n2n2);
                
                // vk - membrane
                Bmat_mem.left_multiply(tmp_mat1_n1n2, material_A_mat);
                vk_dwdxi_mat.get_transpose(tmp_mat3);
                tmp_mat3.right_multiply(tmp_mat1_n1n2);
                Bmat_vk.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat3);
                local_jac.add(-JxW[qp], tmp_mat2_n2n2);
                
                // vk - vk
                tmp_mat3.resize(2, n2);
                Bmat_vk.left_multiply(tmp_mat3, stress);
                Bmat_vk.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat3);
                local_jac.add(-JxW[qp], tmp_mat2_n2n2);
                
                tmp_mat3.resize(vk_dwdxi_mat.m(), n2);
                Bmat_vk.left_multiply(tmp_mat3, vk_dwdxi_mat);
                tmp_mat3.left_multiply(material_A_mat);
                tmp_mat3.left_multiply_transpose(vk_dwdxi_mat);
                Bmat_vk.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat3);
                local_jac.add(-JxW[qp], tmp_mat2_n2n2);
                
                // bending - vk
                tmp_mat3.resize(vk_dwdxi_mat.m(), n2);
                Bmat_vk.left_multiply(tmp_mat3, vk_dwdxi_mat);
                tmp_mat3.left_multiply(material_B_mat);
                Bmat_bend.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat3);
                local_jac.add(-JxW[qp], tmp_mat2_n2n2);
                
                // vk - bending
                Bmat_bend.left_multiply(tmp_mat1_n1n2, material_B_mat);
                vk_dwdxi_mat.get_transpose(tmp_mat3);
                tmp_mat3.right_multiply(tmp_mat1_n1n2);
                Bmat_vk.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat3);
                local_jac.add(-JxW[qp], tmp_mat2_n2n2);
            }
            
            // bending - membrane
            Bmat_mem.left_multiply(tmp_mat1_n1n2, material_B_mat);
            Bmat_bend.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat1_n1n2);
            local_jac.add(-JxW[qp], tmp_mat2_n2n2);
            
            // membrane - bending
            Bmat_bend.left_multiply(tmp_mat1_n1n2, material_B_mat);
            Bmat_mem.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat1_n1n2);
            local_jac.add(-JxW[qp], tmp_mat2_n2n2);
            
            // bending - bending
            Bmat_bend.left_multiply(tmp_mat1_n1n2, material_D_mat);
            Bmat_bend.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat1_n1n2);
            local_jac.add(-JxW[qp], tmp_mat2_n2n2);
        }
    }
}




bool
MAST::BendingStructuralElem::surface_pressure_force(bool request_jacobian,
                                                  DenseVector<Real> &f,
                                                  DenseMatrix<Real> &jac,
                                                  const unsigned int side,
                                                  MAST::BoundaryCondition &p) {
#ifndef LIBMESH_USE_COMPLEX_NUMBERS
    libmesh_assert(!follower_forces); // not implemented yet for follower forces
    
    FEMOperatorMatrix Bmat;
    
    // get the function from this boundary condition
    libMesh::FunctionBase<Number>& func = p.function();
    std::auto_ptr<FEBase> fe;
    std::auto_ptr<QBase> qrule;
    _get_side_fe_and_qrule(_local_elem->local_elem(), side, fe, qrule);
    
    const std::vector<Real> &JxW = fe->get_JxW();
    
    // Physical location of the quadrature points
    const std::vector<Point>& qpoint = fe->get_xyz();
    const std::vector<std::vector<Real> >& phi = fe->get_phi();
    const unsigned int n_phi = (unsigned int)phi.size();
    const unsigned int n1=3, n2=6*n_phi;
    
    // boundary normals
    const std::vector<Point>& face_normals = fe->get_normals();
    Real press;
    
    DenseVector<Real> phi_vec, force, local_f, tmp_vec_n2;
    phi_vec.resize(n_phi); force.resize(2*n1); local_f.resize(n2);
    tmp_vec_n2.resize(n2);
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++)
    {
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        Bmat.reinit(2*n1, phi_vec);
        
        // get pressure value
        press = func(qpoint[qp], _system.time);
        
        // calculate force
        for (unsigned int i_dim=0; i_dim<n1; i_dim++)
            force(i_dim) = press * face_normals[qp](i_dim);
        
        Bmat.vector_mult_transpose(tmp_vec_n2, force);
        
        local_f.add(-JxW[qp], tmp_vec_n2);
    }
    
    // now transform to the global system and add
    transform_to_global_system(local_f, tmp_vec_n2);
    f.add(1., tmp_vec_n2);

#endif // LIBMESH_USE_COMPLEX_NUMBERS
    return (request_jacobian && follower_forces);
}



bool
MAST::BendingStructuralElem::surface_pressure_force(bool request_jacobian,
                                                    DenseVector<Real> &f,
                                                    DenseMatrix<Real> &jac,
                                                    MAST::BoundaryCondition &p) {
#ifndef LIBMESH_USE_COMPLEX_NUMBERS
    libmesh_assert(!follower_forces); // not implemented yet for follower forces
    
    FEMOperatorMatrix Bmat;
    
    // get the function from this boundary condition
    libMesh::FunctionBase<Number>& func = p.function();
    const std::vector<Real> &JxW = _fe->get_JxW();
    
    // Physical location of the quadrature points
    const std::vector<Point>& qpoint = _fe->get_xyz();
    const std::vector<std::vector<Real> >& phi = _fe->get_phi();
    const unsigned int n_phi = (unsigned int)phi.size();
    const unsigned int n1=3, n2=6*n_phi;
    
    // normal for face integration
    Point normal;
    // direction of pressure assumed to be normal (along local z-axis)
    // to the element face for 2D and along local y-axis for 1D element.
    normal(_elem.dim()) = 1.;
    
    Real press;
    
    DenseVector<Real> phi_vec, force, local_f, tmp_vec_n2;
    phi_vec.resize(n_phi); force.resize(2*n1); local_f.resize(n2);
    tmp_vec_n2.resize(n2);
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++)
    {
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        Bmat.reinit(2*n1, phi_vec);
        
        // get pressure value
        press = func(qpoint[qp], _system.time);
        
        // calculate force
        for (unsigned int i_dim=0; i_dim<n1; i_dim++)
            force(i_dim) = press * normal(i_dim);
        
        Bmat.vector_mult_transpose(tmp_vec_n2, force);
        
        local_f.add(-JxW[qp], tmp_vec_n2);
    }
    
    // now transform to the global system and add
    transform_to_global_system(local_f, tmp_vec_n2);
    f.add(1., tmp_vec_n2);
    
#endif // LIBMESH_USE_COMPLEX_NUMBERS
    return (request_jacobian && follower_forces);
}


bool
MAST::BendingStructuralElem::surface_pressure_force_sensitivity(bool request_jacobian,
                                                              DenseVector<Real> &f,
                                                              DenseMatrix<Real> &jac,
                                                              const unsigned int side,
                                                              MAST::BoundaryCondition &p) {
    libmesh_assert(!follower_forces); // not implemented yet for follower forces
    
    // currently not implemented.
    return false;
    
    // this should be true if the function is called
    libmesh_assert(this->sensitivity_params);
    libmesh_assert(!this->sensitivity_params->shape_sensitivity()); // this is not implemented for now
    libmesh_assert(this->sensitivity_params->total_order() == 1); // only first order sensitivity
    
    // check if the material property or the provided exterior
    // values, like temperature, are functions of the sensitivity parameter
    bool calculate = false;
    if (_temperature)
        calculate = calculate || _temperature->depends_on(*(this->sensitivity_params));
    calculate = calculate || _property.depends_on(*(this->sensitivity_params));
    
    // nothing to be calculated if the element does not depend on the
    // sensitivity parameter.
    if (!calculate)
        return false;
    
    
    return (request_jacobian && follower_forces);
}



bool
MAST::BendingStructuralElem::surface_pressure_force_sensitivity(bool request_jacobian,
                                                              DenseVector<Real> &f,
                                                              DenseMatrix<Real> &jac,
                                                              MAST::BoundaryCondition &p) {
    libmesh_assert(!follower_forces); // not implemented yet for follower forces
    
    // currently not implemented.
    return false;
    
    // this should be true if the function is called
    libmesh_assert(this->sensitivity_params);
    libmesh_assert(!this->sensitivity_params->shape_sensitivity()); // this is not implemented for now
    libmesh_assert(this->sensitivity_params->total_order() == 1); // only first order sensitivity
    
    // check if the material property or the provided exterior
    // values, like temperature, are functions of the sensitivity parameter
    bool calculate = false;
    if (_temperature)
        calculate = calculate || _temperature->depends_on(*(this->sensitivity_params));
    calculate = calculate || _property.depends_on(*(this->sensitivity_params));
    
    // nothing to be calculated if the element does not depend on the
    // sensitivity parameter.
    if (!calculate)
        return false;
    
    
    return (request_jacobian && follower_forces);
}



std::auto_ptr<MAST::LocalElemBase>
MAST::build_local_elem(MAST::StructuralElementBase &elem) {
    
    std::auto_ptr<MAST::LocalElemBase> rval;
    switch (elem.elem().dim()) {
        case 1:
            rval.reset(new MAST::Local1DElem(elem.elem()));
            break;
            
        case 2:
            rval.reset(new MAST::Local2DElem(elem.elem()));
            break;
            
        default:
            libmesh_error();
    }
    
    return rval;
}

