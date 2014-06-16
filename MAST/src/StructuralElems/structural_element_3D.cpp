//
//  structural_element_3D.cpp
//  MAST
//
//  Created by Manav Bhatia on 11/14/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

// MAST includes
#include "StructuralElems/structural_element_3D.h"
#include "PropertyCards/element_property_card_base.h"
#include "Numerics/fem_operator_matrix.h"
#include "BoundaryConditions/temperature.h"


bool
MAST::StructuralElement3D::internal_force (bool request_jacobian,
                                           DenseRealVector& f,
                                           DenseRealMatrix& jac,
                                           bool if_ignore_ho_jac)
{
    FEMOperatorMatrix Bmat;
    
    const std::vector<libMesh::Real>& JxW = _fe->get_JxW();
    const std::vector<Point>& xyz = _fe->get_xyz();
    const unsigned int n_phi = (unsigned int)JxW.size();
    const unsigned int n1=6, n2=6*n_phi;
    DenseRealMatrix material_mat, tmp_mat1_n1n2, tmp_mat2_n2n2;
    DenseRealVector tmp_vec1_n1, tmp_vec2_n2;
    
    tmp_mat1_n1n2.resize(n1, n2); tmp_mat2_n2n2.resize(n2, n2);
    tmp_vec1_n1.resize(n1); tmp_vec2_n2.resize(n2);
    
    
    Bmat.reinit(6, _system.n_vars(), _elem.n_nodes()); // six stress-strain components
    
    std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > > mat_stiff
    (_property.get_property(MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX,
                            *this).release());
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        this->initialize_strain_operator(qp, Bmat);
        
        // get the material matrix
        (*mat_stiff)(xyz[qp], _system.time, material_mat);
        
        // calculate the stress
        Bmat.left_multiply(tmp_mat1_n1n2, material_mat);
        tmp_mat1_n1n2.vector_mult(tmp_vec1_n1, local_solution); // this is stress
        
        // now calculate the internal force vector
        Bmat.vector_mult_transpose(tmp_vec2_n2, tmp_vec1_n1);
        f.add(-JxW[qp], tmp_vec2_n2);
        
        if (request_jacobian) {
            
            Bmat.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat1_n1n2);
            jac.add(-JxW[qp], tmp_mat2_n2n2);
        }
    }
    
    return request_jacobian;
}



bool
MAST::StructuralElement3D::internal_force_sensitivity (bool request_jacobian,
                                                       DenseRealVector& f,
                                                       DenseRealMatrix& jac,
                                                       bool if_ignore_ho_jac)
{
    // this should be true if the function is called
    libmesh_assert(this->sensitivity_param);
    libmesh_assert(!this->sensitivity_param->is_shape_parameter()); // this is not implemented for now
    
    
    // check if the material property or the provided exterior
    // values, like temperature, are functions of the sensitivity parameter
    bool calculate = false;
    calculate = calculate || _property.depends_on(*(this->sensitivity_param));
    
    // nothing to be calculated if the element does not depend on the
    // sensitivity parameter.
    if (!calculate)
        return false;
    
    FEMOperatorMatrix Bmat;
    
    const std::vector<libMesh::Real>& JxW = _fe->get_JxW();
    const std::vector<Point>& xyz = _fe->get_xyz();
    const unsigned int n_phi = (unsigned int)JxW.size();
    const unsigned int n1=6, n2=6*n_phi;
    DenseRealMatrix material_mat, tmp_mat1_n1n2, tmp_mat2_n2n2;
    DenseRealVector tmp_vec1_n1, tmp_vec2_n2;
    
    tmp_mat1_n1n2.resize(n1, n2); tmp_mat2_n2n2.resize(n2, n2);
    tmp_vec1_n1.resize(n1); tmp_vec2_n2.resize(n2);
    
    
    Bmat.reinit(6, _system.n_vars(), _elem.n_nodes()); // six stress-strain components

    std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > > mat_stiff
    (_property.get_property(MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX,
                            *this).release());

    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        this->initialize_strain_operator(qp, Bmat);
        
        // get the material matrix
        mat_stiff->total(*this->sensitivity_param,
                         xyz[qp], _system.time, material_mat);

        // calculate the stress
        Bmat.left_multiply(tmp_mat1_n1n2, material_mat);
        tmp_mat1_n1n2.vector_mult(tmp_vec1_n1, local_solution); // this is stress
        
        // now calculate the internal force vector
        Bmat.vector_mult_transpose(tmp_vec2_n2, tmp_vec1_n1);
        f.add(-JxW[qp], tmp_vec2_n2);
        
        if (request_jacobian) {
            
            Bmat.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat1_n1n2);
            jac.add(-JxW[qp], tmp_mat2_n2n2);
        }
    }
    
    return request_jacobian;
}




bool
MAST::StructuralElement3D::prestress_force (bool request_jacobian,
                                           DenseRealVector& f,
                                           DenseRealMatrix& jac)
{
    FEMOperatorMatrix Bmat;
    
    const std::vector<libMesh::Real>& JxW = _fe->get_JxW();
    const std::vector<Point>& xyz = _fe->get_xyz();
    const unsigned int n_phi = (unsigned int)JxW.size();
    const unsigned int n1=6, n2=6*n_phi;
    DenseRealMatrix prestress_mat_A, tmp_mat1_n1n2, tmp_mat2_n2n2;
    DenseRealVector  tmp_vec1_n1, tmp_vec2_n2, prestress_vec_A;
    
    tmp_mat1_n1n2.resize(n1, n2); tmp_mat2_n2n2.resize(n2, n2);
    tmp_vec1_n1.resize(n1); tmp_vec2_n2.resize(n2);
    
    Bmat.reinit(6, _system.n_vars(), _elem.n_nodes()); // six stress-strain components
    
    std::auto_ptr<MAST::SectionIntegratedPrestressMatrixBase>
    prestress_A
    (dynamic_cast<MAST::SectionIntegratedPrestressMatrixBase*>
     (_property.get_property(MAST::SECTION_INTEGRATED_PRESTRESS_A_MATRIX,
                             *this).release()));
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        this->initialize_strain_operator(qp, Bmat);
        
        // get the material matrix
        (*prestress_A)(xyz[qp], _system.time, prestress_mat_A);
        prestress_A->convert_to_vector(prestress_mat_A, prestress_vec_A);
        
        // now calculate the force vector
        Bmat.vector_mult_transpose(tmp_vec2_n2, prestress_vec_A);
        f.add(-JxW[qp], tmp_vec2_n2);
        
        if (request_jacobian) {
            // nothing to be done here unless a nonlinear strain is included
        }
    }
    
    return request_jacobian;
}



bool
MAST::StructuralElement3D::prestress_force_sensitivity (bool request_jacobian,
                                                        DenseRealVector& f,
                                                        DenseRealMatrix& jac)
{
    // this should be true if the function is called
    libmesh_assert(this->sensitivity_param);
    libmesh_assert(!this->sensitivity_param->is_shape_parameter()); // this is not implemented for now
    
    
    // check if the material property or the provided exterior
    // values, like temperature, are functions of the sensitivity parameter
    bool calculate = false;
    calculate = calculate || _property.depends_on(*(this->sensitivity_param));
    
    // nothing to be calculated if the element does not depend on the
    // sensitivity parameter.
    if (!calculate)
        return false;
    
    FEMOperatorMatrix Bmat;
    
    const std::vector<libMesh::Real>& JxW = _fe->get_JxW();
    const std::vector<Point>& xyz = _fe->get_xyz();
    const unsigned int n_phi = (unsigned int)JxW.size();
    const unsigned int n1=6, n2=6*n_phi;
    DenseRealMatrix prestress_mat_A, tmp_mat1_n1n2, tmp_mat2_n2n2;
    DenseRealVector prestress_vec_A, tmp_vec1_n1, tmp_vec2_n2;
    
    tmp_mat1_n1n2.resize(n1, n2); tmp_mat2_n2n2.resize(n2, n2);
    tmp_vec1_n1.resize(n1); tmp_vec2_n2.resize(n2);
    
    
    Bmat.reinit(6, _system.n_vars(), _elem.n_nodes()); // six stress-strain components
    
    std::auto_ptr<MAST::SectionIntegratedPrestressMatrixBase>
    prestress_A
    (dynamic_cast<MAST::SectionIntegratedPrestressMatrixBase*>
     (_property.get_property(MAST::SECTION_INTEGRATED_PRESTRESS_A_MATRIX,
                             *this).release()));
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        this->initialize_strain_operator(qp, Bmat);
        
        // get the material matrix sensitivity
        prestress_A->total(*this->sensitivity_param,
                           xyz[qp], _system.time, prestress_mat_A);
        prestress_A->convert_to_vector(prestress_mat_A, prestress_vec_A);
        
        // now calculate the force vector
        Bmat.vector_mult_transpose(tmp_vec2_n2, prestress_vec_A);
        f.add(-JxW[qp], tmp_vec2_n2);
        
        if (request_jacobian) {
            // nothing to be done here unless a nonlinear strain is included
        }
    }
    
    return request_jacobian;
}





bool
MAST::StructuralElement3D::thermal_force (bool request_jacobian,
                                          DenseRealVector& f,
                                          DenseRealMatrix& jac,
                                          MAST::BoundaryCondition& p)
{
    FEMOperatorMatrix Bmat;
    
    const std::vector<libMesh::Real>& JxW = _fe->get_JxW();
    const std::vector<Point>& xyz = _fe->get_xyz();
    const unsigned int n_phi = (unsigned int)_fe->get_phi().size();
    const unsigned int n1= 6, n2=6*n_phi;
    DenseRealMatrix material_exp_A_mat, tmp_mat1_n1n2, tmp_mat2_n2n2, stress;
    DenseRealVector  tmp_vec1_n1, tmp_vec2_n1, tmp_vec3_n2, delta_t;
    
    tmp_mat1_n1n2.resize(n1, n2); tmp_mat2_n2n2.resize(n2, n2);
    tmp_vec1_n1.resize(n1); tmp_vec2_n1.resize(n1);
    tmp_vec3_n2.resize(n2); delta_t.resize(1);
    
    
    
    Bmat.reinit(n1, _system.n_vars(), n_phi); // three stress-strain components
    
    std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > > mat
    (_property.get_property(MAST::SECTION_INTEGRATED_MATERIAL_THERMAL_EXPANSION_A_MATRIX,
                            *this).release());
    MAST::FieldFunction<libMesh::Real>& temp_func =
    dynamic_cast<MAST::FieldFunction<libMesh::Real>&>(p.function());
    MAST::FieldFunction<libMesh::Real>& ref_temp_func =
    dynamic_cast<MAST::FieldFunction<libMesh::Real>&>
    (dynamic_cast<MAST::Temperature&>(p).reference_temperature_function());
    
    libMesh::Real t, t0;
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        (*mat)(xyz[qp], _system.time, material_exp_A_mat);
        temp_func(xyz[qp], _system.time, t);
        ref_temp_func(xyz[qp], _system.time, t0);
        delta_t(0) = t-t0;
        
        material_exp_A_mat.vector_mult(tmp_vec1_n1, delta_t); // [C]{alpha (T - T0)}
        
        this->initialize_strain_operator(qp, Bmat);
        
        Bmat.vector_mult_transpose(tmp_vec3_n2, tmp_vec1_n1);
        f.add(JxW[qp], tmp_vec3_n2);
    }
    
    // Jacobian contribution from von Karman strain
    return false;
}




bool
MAST::StructuralElement3D::thermal_force_sensitivity (bool request_jacobian,
                                                      DenseRealVector& f,
                                                      DenseRealMatrix& jac,
                                                      MAST::BoundaryCondition& p)
{
    FEMOperatorMatrix Bmat;
    
    const std::vector<libMesh::Real>& JxW = _fe->get_JxW();
    const std::vector<Point>& xyz = _fe->get_xyz();
    const unsigned int n_phi = (unsigned int)_fe->get_phi().size();
    const unsigned int n1= 6, n2=6*n_phi;
    DenseRealMatrix material_exp_A_mat, material_exp_A_mat_sens,
    tmp_mat1_n1n2, tmp_mat2_n2n2, stress;
    DenseRealVector  tmp_vec1_n1, tmp_vec2_n1, tmp_vec3_n2, delta_t,
    delta_t_sens;
    
    tmp_mat1_n1n2.resize(n1, n2); tmp_mat2_n2n2.resize(n2, n2);
    tmp_vec1_n1.resize(n1); tmp_vec2_n1.resize(n1);
    tmp_vec3_n2.resize(n2); delta_t.resize(1); delta_t_sens.resize(1);
    
    
    
    Bmat.reinit(n1, _system.n_vars(), n_phi); // three stress-strain components
    
    std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > > mat
    (_property.get_property(MAST::SECTION_INTEGRATED_MATERIAL_THERMAL_EXPANSION_A_MATRIX,
                            *this).release());
    MAST::FieldFunction<libMesh::Real>& temp_func =
    dynamic_cast<MAST::FieldFunction<libMesh::Real>&>(p.function());
    MAST::FieldFunction<libMesh::Real>& ref_temp_func =
    dynamic_cast<MAST::FieldFunction<libMesh::Real>&>
    (dynamic_cast<MAST::Temperature&>(p).reference_temperature_function());
    
    libMesh::Real t, t0, t_sens;
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        (*mat)(xyz[qp], _system.time, material_exp_A_mat);
        mat->total(*this->sensitivity_param,
                   xyz[qp], _system.time, material_exp_A_mat_sens);

        temp_func(xyz[qp], _system.time, t);
        temp_func.total(*this->sensitivity_param,
                        xyz[qp], _system.time, t_sens);
        ref_temp_func(xyz[qp], _system.time, t0);
        delta_t(0) = t-t0;
        delta_t_sens(0) = t_sens;
        
        
        this->initialize_strain_operator(qp, Bmat);
        
        material_exp_A_mat.vector_mult(tmp_vec1_n1, delta_t_sens); // [C]{alpha dT/dp}
        Bmat.vector_mult_transpose(tmp_vec3_n2, tmp_vec1_n1);
        f.add(JxW[qp], tmp_vec3_n2);
        
        material_exp_A_mat_sens.vector_mult(tmp_vec1_n1, delta_t); // d([C].{alpha})/dp (T - T0)}
        Bmat.vector_mult_transpose(tmp_vec3_n2, tmp_vec1_n1);
        f.add(JxW[qp], tmp_vec3_n2);
    }
    
    // Jacobian contribution from von Karman strain
    return false;
}




void
MAST::StructuralElement3D::initialize_strain_operator(const unsigned int qp,
                                                      FEMOperatorMatrix& Bmat) {
    const std::vector<std::vector<RealVectorValue> >& dphi = _fe->get_dphi();
    
    unsigned int n_phi = (unsigned int)dphi.size();
    DenseRealVector phi; phi.resize(n_phi);
    
    // now set the shape function values
    // dN/dx
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
        phi(i_nd) = dphi[i_nd][qp](0);
    Bmat.set_shape_function(0, 0, phi); //  epsilon_xx = du/dx
    Bmat.set_shape_function(3, 1, phi); //  gamma_xy = dv/dx + ...
    Bmat.set_shape_function(5, 2, phi); //  gamma_zx = dw/dx + ...
    
    // dN/dy
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
        phi(i_nd) = dphi[i_nd][qp](1);
    Bmat.set_shape_function(1, 1, phi); //  epsilon_yy = dv/dy
    Bmat.set_shape_function(3, 0, phi); //  gamma_xy = du/dy + ...
    Bmat.set_shape_function(4, 2, phi); //  gamma_yz = dw/dy + ...
    
    // dN/dz
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
        phi(i_nd) = dphi[i_nd][qp](2);
    Bmat.set_shape_function(2, 2, phi); //  epsilon_xx = dw/dz
    Bmat.set_shape_function(4, 1, phi); //  gamma_xy = dv/dz + ...
    Bmat.set_shape_function(5, 0, phi); //  gamma_zx = du/dz + ...
    
}


