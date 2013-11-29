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
#include "ThermalElems/temperature_function.h"
#include "BoundaryConditions/boundary_condition.h"


bool
MAST::StructuralElement3D::internal_force (bool request_jacobian,
                                     DenseVector<Real>& f,
                                     DenseMatrix<Real>& jac)
{
    FEMOperatorMatrix Bmat;
    
    const std::vector<Real>& JxW = _fe->get_JxW();
    const std::vector<Point>& xyz = _fe->get_xyz();
    const unsigned int n_phi = (unsigned int)JxW.size();
    const unsigned int n1=6, n2=6*n_phi;
    DenseMatrix<Real> material_mat, tmp_mat1_n1n2, tmp_mat2_n2n2;
    DenseVector<Real>  phi, tmp_vec1_n1, tmp_vec2_n2;
    
    tmp_mat1_n1n2.resize(n1, n2); tmp_mat2_n2n2.resize(n2, n2);
    phi.resize(n_phi); tmp_vec1_n1.resize(n1); tmp_vec2_n2.resize(n2);
    
    
    Bmat.reinit(6, _system.n_vars(), _elem.n_nodes()); // six stress-strain components
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        this->initialize_strain_operator(qp, Bmat);
        
        // if temperature is specified, the initialize it to the current location
        if (_temperature)
            _temperature->initialize(xyz[qp]);
        
        // get the material matrix
        _property.calculate_matrix(_elem,
                                   MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX,
                                   material_mat);
        
        // calculate the stress
        Bmat.left_multiply(tmp_mat1_n1n2, material_mat);
        tmp_mat1_n1n2.vector_mult(tmp_vec1_n1, local_solution); // this is stress
        
        // now calculate the internal force vector
        Bmat.vector_mult_transpose(tmp_vec2_n2, tmp_vec1_n1);
        f.add(-JxW[qp], tmp_vec2_n2);
        
        // add the prestress
        
        if (request_jacobian) {
            
            Bmat.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat1_n1n2);
            jac.add(-JxW[qp], tmp_mat2_n2n2);
        }
    }
    
    return request_jacobian;
}



bool
MAST::StructuralElement3D::internal_force_sensitivity (bool request_jacobian,
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
    
    FEMOperatorMatrix Bmat;
    
    const std::vector<Real>& JxW = _fe->get_JxW();
    const std::vector<Point>& xyz = _fe->get_xyz();
    const unsigned int n_phi = (unsigned int)JxW.size();
    const unsigned int n1=6, n2=6*n_phi;
    DenseMatrix<Real> material_mat, tmp_mat1_n1n2, tmp_mat2_n2n2;
    DenseVector<Real>  phi, tmp_vec1_n1, tmp_vec2_n2;
    
    tmp_mat1_n1n2.resize(n1, n2); tmp_mat2_n2n2.resize(n2, n2);
    phi.resize(n_phi); tmp_vec1_n1.resize(n1); tmp_vec2_n2.resize(n2);
    
    
    Bmat.reinit(6, _system.n_vars(), _elem.n_nodes()); // six stress-strain components
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        this->initialize_strain_operator(qp, Bmat);
        
        // if temperature is specified, the initialize it to the current location
        if (_temperature)
            _temperature->initialize(xyz[qp]);
        
        // get the material matrix
        _property.calculate_matrix_sensitivity(_elem,
                                               MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX,
                                               material_mat,
                                               *(this->sensitivity_params));
        
        // calculate the stress
        Bmat.left_multiply(tmp_mat1_n1n2, material_mat);
        tmp_mat1_n1n2.vector_mult(tmp_vec1_n1, local_solution); // this is stress
        
        // now calculate the internal force vector
        Bmat.vector_mult_transpose(tmp_vec2_n2, tmp_vec1_n1);
        f.add(-JxW[qp], tmp_vec2_n2);
        
        // add the prestress
        
        if (request_jacobian) {
            
            Bmat.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat1_n1n2);
            jac.add(-JxW[qp], tmp_mat2_n2n2);
        }
    }
    
    return request_jacobian;
}




bool
MAST::StructuralElement3D::thermal_force (bool request_jacobian,
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
    DenseVector<Real>  phi, temperature, tmp_vec1_n1, tmp_vec2_n1, tmp_vec3_n2;
    
    phi.resize(n_phi); tmp_vec1_n1.resize(n1); tmp_vec2_n1.resize(n1);
    tmp_vec3_n2.resize(n2); temperature.resize(6);
    
    Bmat.reinit(6, _system.n_vars(), _elem.n_nodes()); // six stress-strain components
    
    Real temperature_value, ref_temperature;
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        this->initialize_strain_operator(qp, Bmat);
        
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
        f.add(JxW[qp], tmp_vec3_n2);
    }
    
    return false;
}




bool
MAST::StructuralElement3D::thermal_force_sensitivity (bool request_jacobian,
                                                DenseVector<Real>& f,
                                                DenseMatrix<Real>& jac)
{
    
    libmesh_error(); // to be implemented
    
    return false;
}




void
MAST::StructuralElement3D::initialize_strain_operator(const unsigned int qp,
                                                FEMOperatorMatrix& Bmat) {
    const std::vector<std::vector<RealVectorValue> >& dphi = _fe->get_dphi();
    
    unsigned int n_phi = (unsigned int)dphi.size();
    DenseVector<Real> phi; phi.resize(n_phi);
    
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


bool
MAST::StructuralElement3D::surface_pressure_force(bool request_jacobian,
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
    _get_side_fe_and_qrule(_elem, side, fe, qrule);
    
    const std::vector<Real> &JxW = fe->get_JxW();
    
    // Physical location of the quadrature points
    const std::vector<Point>& qpoint = fe->get_xyz();
    const std::vector<std::vector<Real> >& phi = fe->get_phi();
    const unsigned int n_phi = (unsigned int)phi.size();
    const unsigned int n1=3, n2=6*n_phi;
    
    // boundary normals
    const std::vector<Point>& face_normals = fe->get_normals();
    Real press;
    
    DenseVector<Real> phi_vec, force, tmp_vec_n2;
    phi_vec.resize(n_phi); force.resize(2*n1); tmp_vec_n2.resize(n2);
    
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
        
        f.add(-JxW[qp], tmp_vec_n2);
    }
#endif // LIBMESH_USE_COMPLEX_NUMBERS
    return (request_jacobian && follower_forces);
}



bool
MAST::StructuralElement3D::surface_pressure_force_sensitivity(bool request_jacobian,
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


