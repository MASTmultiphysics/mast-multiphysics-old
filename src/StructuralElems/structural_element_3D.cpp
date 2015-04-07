/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2014  Manav Bhatia
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

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
    FEMOperatorMatrix Bmat, Bincmat;
    
    const std::vector<Real>& JxW = _fe->get_JxW();
    const std::vector<libMesh::Point>& xyz = _fe->get_xyz();
    const unsigned int n_phi = (unsigned int)_fe->n_shape_functions();
    const unsigned int n1=6, n2=3*n_phi;
    DenseRealMatrix material_mat, mat1_n1n2, mat2_n2n2,
    m_n13(n1,3), m_33(3,3), m_3n2(3,n2), m_n23(n2,3),
    Kda(n2,3), Kad(3,n2), Kaa(3,3); // a is the incompatible modes adn d are the deformation modes
    DenseRealVector vec1_n1, vec2_n2, v1_3(3), v2_3(3);
    
    mat1_n1n2.resize(n1, n2); mat2_n2n2.resize(n2, n2);
    vec1_n1.resize(n1); vec2_n2.resize(n2);
    
    
    Bmat.reinit(n1, 3, _elem.n_nodes()); // six stress-strain components
    Bincmat.reinit(n1, 3, 1);            // six stress-strain components
    
    std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > > mat_stiff
    (_property.get_property(MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX,
                            *this).release());
        
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        this->initialize_strain_operator(qp, Bmat);
        this->initialize_incompatible_strain_operator(qp, Bincmat);
        
        // get the material matrix
        (*mat_stiff)(xyz[qp], _system.time, material_mat);
        
        // calculate the stress
        Bmat.left_multiply(mat1_n1n2, material_mat); // D B
        mat1_n1n2.vector_mult(vec1_n1, local_solution); // this is stress
        
        // calculate the contribution to incompatible mode matrix
        Bincmat.right_multiply_transpose(m_3n2, mat1_n1n2);  // Ba^T D Bd
        Kad.add(JxW[qp], m_3n2);
        
        Bincmat.left_multiply(m_n13, material_mat);
        Bmat.right_multiply_transpose(m_n23, m_n13);
        Kda.add(JxW[qp], m_n23);
        
        Bincmat.right_multiply_transpose(m_33, m_n13);
        Kaa.add(JxW[qp], m_33);
        
        // now calculate the internal force vector
        Bmat.vector_mult_transpose(vec2_n2, vec1_n1);
        f.add(-JxW[qp], vec2_n2);
        
        if (request_jacobian) {
            
            Bmat.right_multiply_transpose(mat2_n2n2, mat1_n1n2);
            jac.add(-JxW[qp], mat2_n2n2);
        }
    }
    
    // add the force and Jacobian contributions from the static condensation
    // of the Kaa matrix
    // calculate inv(Kaa)
    for (unsigned int i=0; i<3; i++) {
        v1_3.zero(); v1_3(i) = 1.;
        Kaa.lu_solve(v1_3, v2_3);
        m_33.set_column(i, v2_3);
    }

    Kad.left_multiply(m_33);
    Kad.left_multiply(Kda);

    // contribution to residual vector
    Kad.vector_mult(vec2_n2, local_solution);
    f.add(1., vec2_n2);
    
    // contribution to Jacobian
    if (request_jacobian) jac.add(1., Kad);
    
    // place small values at diagonal of the rotational terms
    if (request_jacobian) {
        for (unsigned int i=n2/2; i<n2; i++)
            jac(i,i) += 1.0e-12;
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
    
    const std::vector<Real>& JxW = _fe->get_JxW();
    const std::vector<libMesh::Point>& xyz = _fe->get_xyz();
    const unsigned int n_phi = (unsigned int)_fe->n_shape_functions();
    const unsigned int n1=6, n2=6*n_phi;
    DenseRealMatrix material_mat, mat1_n1n2, mat2_n2n2;
    DenseRealVector vec1_n1, vec2_n2;
    
    mat1_n1n2.resize(n1, n2); mat2_n2n2.resize(n2, n2);
    vec1_n1.resize(n1); vec2_n2.resize(n2);
    
    
    Bmat.reinit(n1, 3, _elem.n_nodes()); // six stress-strain components

    std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > > mat_stiff
    (_property.get_property(MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX,
                            *this).release());

    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        this->initialize_strain_operator(qp, Bmat);
        
        // get the material matrix
        mat_stiff->total(*this->sensitivity_param,
                         xyz[qp], _system.time, material_mat);

        // calculate the stress
        Bmat.left_multiply(mat1_n1n2, material_mat);
        mat1_n1n2.vector_mult(vec1_n1, local_solution); // this is stress
        
        // now calculate the internal force vector
        Bmat.vector_mult_transpose(vec2_n2, vec1_n1);
        f.add(-JxW[qp], vec2_n2);
        
        if (request_jacobian) {
            
            Bmat.right_multiply_transpose(mat2_n2n2, mat1_n1n2);
            jac.add(-JxW[qp], mat2_n2n2);
        }
    }
    
    return request_jacobian;
}




bool
MAST::StructuralElement3D::prestress_force (bool request_jacobian,
                                           DenseRealVector& f,
                                           DenseRealMatrix& jac)
{
    if (!_property.if_prestressed())
        return false;

    FEMOperatorMatrix Bmat;
    
    const std::vector<Real>& JxW = _fe->get_JxW();
    const std::vector<libMesh::Point>& xyz = _fe->get_xyz();
    const unsigned int n_phi = (unsigned int)_fe->n_shape_functions();
    const unsigned int n1=6, n2=6*n_phi;
    DenseRealMatrix prestress_mat_A, mat1_n1n2, mat2_n2n2;
    DenseRealVector  vec1_n1, vec2_n2, prestress_vec_A;
    
    mat1_n1n2.resize(n1, n2); mat2_n2n2.resize(n2, n2);
    vec1_n1.resize(n1); vec2_n2.resize(n2);
    
    Bmat.reinit(n1, 3, _elem.n_nodes()); // six stress-strain components
    
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
        Bmat.vector_mult_transpose(vec2_n2, prestress_vec_A);
        f.add(-JxW[qp], vec2_n2);
        
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
    
    const std::vector<Real>& JxW = _fe->get_JxW();
    const std::vector<libMesh::Point>& xyz = _fe->get_xyz();
    const unsigned int n_phi = (unsigned int)_fe->n_shape_functions();
    const unsigned int n1=6, n2=6*n_phi;
    DenseRealMatrix prestress_mat_A, mat1_n1n2, mat2_n2n2;
    DenseRealVector prestress_vec_A, vec1_n1, vec2_n2;
    
    mat1_n1n2.resize(n1, n2); mat2_n2n2.resize(n2, n2);
    vec1_n1.resize(n1); vec2_n2.resize(n2);
    
    
    Bmat.reinit(n1, 3, _elem.n_nodes()); // six stress-strain components
    
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
        Bmat.vector_mult_transpose(vec2_n2, prestress_vec_A);
        f.add(-JxW[qp], vec2_n2);
        
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
    
    const std::vector<Real>& JxW = _fe->get_JxW();
    const std::vector<libMesh::Point>& xyz = _fe->get_xyz();
    const unsigned int n_phi = (unsigned int)_fe->n_shape_functions();
    const unsigned int n1=6, n2=3*n_phi;
    DenseRealMatrix material_exp_A_mat, mat1_n1n2, mat2_n2n2, stress;
    DenseRealVector  vec1_n1, vec2_n1, vec3_n2, delta_t;
    
    mat1_n1n2.resize(n1, n2); mat2_n2n2.resize(n2, n2);
    vec1_n1.resize(n1); vec2_n1.resize(n1);
    vec3_n2.resize(n2); delta_t.resize(1);
    
    
    
    Bmat.reinit(n1, 3, n_phi); // three stress-strain components
    
    std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > > mat
    (_property.get_property(MAST::SECTION_INTEGRATED_MATERIAL_THERMAL_EXPANSION_A_MATRIX,
                            *this).release());
    MAST::FieldFunction<Real>& temp_func =
    dynamic_cast<MAST::FieldFunction<Real>&>(p.function());
    MAST::FieldFunction<Real>& ref_temp_func =
    dynamic_cast<MAST::FieldFunction<Real>&>
    (dynamic_cast<MAST::Temperature&>(p).reference_temperature_function());
    
    Real t, t0;
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        (*mat)(xyz[qp], _system.time, material_exp_A_mat);
        temp_func(xyz[qp], _system.time, t);
        ref_temp_func(xyz[qp], _system.time, t0);
        delta_t(0) = t-t0;
        
        material_exp_A_mat.vector_mult(vec1_n1, delta_t); // [C]{alpha (T - T0)}
        
        this->initialize_strain_operator(qp, Bmat);
        
        Bmat.vector_mult_transpose(vec3_n2, vec1_n1);
        
        f.add(JxW[qp], vec3_n2);
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
    
    const std::vector<Real>& JxW = _fe->get_JxW();
    const std::vector<libMesh::Point>& xyz = _fe->get_xyz();
    const unsigned int n_phi = (unsigned int)_fe->n_shape_functions();
    const unsigned int n1= 6, n2=6*n_phi;
    DenseRealMatrix material_exp_A_mat, material_exp_A_mat_sens,
    mat1_n1n2, mat2_n2n2, stress;
    DenseRealVector  vec1_n1, vec2_n1, vec3_n2, delta_t,
    delta_t_sens;
    
    mat1_n1n2.resize(n1, n2); mat2_n2n2.resize(n2, n2);
    vec1_n1.resize(n1); vec2_n1.resize(n1);
    vec3_n2.resize(n2); delta_t.resize(1); delta_t_sens.resize(1);
    
    
    
    Bmat.reinit(n1, 3, n_phi); // three stress-strain components
    
    std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > > mat
    (_property.get_property(MAST::SECTION_INTEGRATED_MATERIAL_THERMAL_EXPANSION_A_MATRIX,
                            *this).release());
    MAST::FieldFunction<Real>& temp_func =
    dynamic_cast<MAST::FieldFunction<Real>&>(p.function());
    MAST::FieldFunction<Real>& ref_temp_func =
    dynamic_cast<MAST::FieldFunction<Real>&>
    (dynamic_cast<MAST::Temperature&>(p).reference_temperature_function());
    
    Real t, t0, t_sens;
    
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
        
        material_exp_A_mat.vector_mult(vec1_n1, delta_t_sens); // [C]{alpha dT/dp}
        Bmat.vector_mult_transpose(vec3_n2, vec1_n1);
        f.add(JxW[qp], vec3_n2);
        
        material_exp_A_mat_sens.vector_mult(vec1_n1, delta_t); // d([C].{alpha})/dp (T - T0)}
        Bmat.vector_mult_transpose(vec3_n2, vec1_n1);
        f.add(JxW[qp], vec3_n2);
    }
    
    // Jacobian contribution from von Karman strain
    return false;
}




void
MAST::StructuralElement3D::initialize_strain_operator(const unsigned int qp,
                                                      FEMOperatorMatrix& Bmat) {
    const std::vector<std::vector<libMesh::RealVectorValue> >& dphi = _fe->get_dphi();
    
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




void
MAST::StructuralElement3D::initialize_incompatible_strain_operator(const unsigned int qp,
                                                                   FEMOperatorMatrix& Bmat) {
    
    DenseRealVector phi_vec; phi_vec.resize(3);

    // get the location of element coordinates
    const std::vector<libMesh::Point>& q_point = _qrule->get_points();
    const Real
    xi  = q_point[qp](0),
    eta = q_point[qp](1),
    phi = q_point[qp](2);
    
    // now set the shape function values
    // dN/dx
    phi_vec(0) = -xi  * dxidx(0,0);
    phi_vec(1) = -eta * dxidx(1,0);
    phi_vec(2) = -phi * dxidx(2,0);
    Bmat.set_shape_function(0, 0, phi_vec); //  epsilon_xx = du/dx
    Bmat.set_shape_function(3, 1, phi_vec); //  gamma_xy = dv/dx + ...
    Bmat.set_shape_function(5, 2, phi_vec); //  gamma_zx = dw/dx + ...
    
    // dN/dy
    phi_vec(0) = -xi  * dxidx(0,1);
    phi_vec(1) = -eta * dxidx(1,1);
    phi_vec(2) = -phi * dxidx(2,1);
    Bmat.set_shape_function(1, 1, phi_vec); //  epsilon_yy = dv/dy
    Bmat.set_shape_function(3, 0, phi_vec); //  gamma_xy = du/dy + ...
    Bmat.set_shape_function(4, 2, phi_vec); //  gamma_yz = dw/dy + ...
    
    // dN/dz
    phi_vec(0) = -xi  * dxidx(0,2);
    phi_vec(1) = -eta * dxidx(1,2);
    phi_vec(2) = -phi * dxidx(2,2);
    Bmat.set_shape_function(2, 2, phi_vec); //  epsilon_xx = dw/dz
    Bmat.set_shape_function(4, 1, phi_vec); //  gamma_xy = dv/dz + ...
    Bmat.set_shape_function(5, 0, phi_vec); //  gamma_zx = du/dz + ...
    
}




void
MAST::StructuralElement3D::_init_incompatible_fe_mapping( const libMesh::Elem& e) {
    
    unsigned int nv = _system.n_vars();
    
    libmesh_assert (nv);
    libMesh::FEType fe_type = _system.variable_type(0); // all variables are assumed to be of same type
    
    
    for (unsigned int i=1; i != nv; ++i)
        libmesh_assert(fe_type == _system.variable_type(i));
    
    // Create an adequate quadrature rule
    std::auto_ptr<libMesh::FEBase> fe(libMesh::FEBase::build(e.dim(), fe_type).release());
    const std::vector<libMesh::Point> pts(1);
    
    fe->get_dxyzdxi();
    fe->get_dxyzdeta();
    fe->get_dxyzdzeta();
    
    fe->reinit(&e, &pts);
    
    dxidx.resize(3,3);
    
    const std::vector<libMesh::RealGradient>&
    dxyzdxi  = fe->get_dxyzdxi(),
    dxyzdeta = fe->get_dxyzdeta(),
    dxyzdphi = fe->get_dxyzdzeta();
    
    dxidx(0,0) =  (dxyzdeta[0](1)*dxyzdphi[0](2)-dxyzdeta[0](2)*dxyzdphi[0](1));
    dxidx(0,1) = -(dxyzdeta[0](0)*dxyzdphi[0](2)-dxyzdeta[0](2)*dxyzdphi[0](0));
    dxidx(0,2) =  (dxyzdeta[0](0)*dxyzdphi[0](1)-dxyzdeta[0](1)*dxyzdphi[0](0));
    
    dxidx(1,0) = -(dxyzdxi[0](1)*dxyzdphi[0](2)-dxyzdxi[0](2)*dxyzdphi[0](1));
    dxidx(1,1) =  (dxyzdxi[0](0)*dxyzdphi[0](2)-dxyzdxi[0](2)*dxyzdphi[0](0));
    dxidx(1,2) = -(dxyzdxi[0](0)*dxyzdphi[0](1)-dxyzdxi[0](1)*dxyzdphi[0](0));
    
    dxidx(2,0) =  (dxyzdxi[0](1)*dxyzdeta[0](2)-dxyzdxi[0](2)*dxyzdeta[0](1));
    dxidx(2,1) = -(dxyzdxi[0](0)*dxyzdeta[0](2)-dxyzdxi[0](2)*dxyzdeta[0](0));
    dxidx(2,2) =  (dxyzdxi[0](0)*dxyzdeta[0](1)-dxyzdxi[0](1)*dxyzdeta[0](0));

}
