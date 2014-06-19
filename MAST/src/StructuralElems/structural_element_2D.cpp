//
//  structural_element_2D.cpp
//  MAST
//
//  Created by Manav Bhatia on 11/14/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//


// MAST includes
#include "StructuralElems/structural_element_2D.h"
#include "PropertyCards/element_property_card_2D.h"
#include "Numerics/fem_operator_matrix.h"
#include "BoundaryConditions/boundary_condition.h"
#include "BoundaryConditions/temperature.h"


void
MAST::Local2DElem::_create_local_elem() {
    
    libmesh_assert(_elem.dim() == 2);
    
    // first node is the origin of the new cs
    // calculate the coordinate system for the plane of the element
    libMesh::Point v1, v2, v3, p;
    v1 = *_elem.get_node(1); v1 -= *_elem.get_node(0); v1 /= v1.size(); // local x
    v2 = *_elem.get_node(2); v2 -= *_elem.get_node(0); v2 /= v2.size();
    v3 = v1.cross(v2); v3 /= v3.size();      // local z
    v2 = v3.cross(v1); v2 /= v2.size();      // local y
    
    _T_mat.resize(3,3);

    _local_elem = libMesh::Elem::build(_elem.type()).release();
    _local_nodes.resize(_elem.n_nodes());
    for (unsigned int i=0; i<_elem.n_nodes(); i++) {
        _local_nodes[i] = new libMesh::Node;
        _local_nodes[i]->set_id() = _elem.get_node(i)->id();
        _local_elem->set_node(i) = _local_nodes[i];
    }

    // now the transformation matrix from old to new cs
    //        an_i vn_i = a_i v_i
    //        an_j = a_i v_i.vn_j  = a_i t_ij = T^t a_i
    //        t_ij = v_i.vn_j
    
    for (unsigned int i=0; i<3; i++) {
        _T_mat(i,0) = v1(i);
        _T_mat(i,1) = v2(i);
        _T_mat(i,2) = v3(i);
    }
    
    // now calculate the new coordinates with respect to the origin
    for (unsigned int i=0; i<_local_nodes.size(); i++) {
        p = *_elem.get_node(i);
        p -= *_elem.get_node(0); // local wrt origin
        for (unsigned int j=0; j<3; j++)
            for (unsigned int k=0; k<3; k++)
                (*_local_nodes[i])(j) += _T_mat(k,j)*p(k);
    }
}



MAST::StructuralElement2D::StructuralElement2D(libMesh::System& sys,
                                               const libMesh::Elem& elem,
                                               const MAST::ElementPropertyCardBase& p):
MAST::BendingStructuralElem(sys, elem, p)
{ }






void
MAST::StructuralElement2D::initialize_direct_strain_operator(const unsigned int qp,
                                                             FEMOperatorMatrix& Bmat) {
    
    const std::vector<std::vector<libMesh::RealVectorValue> >& dphi = _fe->get_dphi();
    
    unsigned int n_phi = (unsigned int)dphi.size();
    DenseRealVector phi; phi.resize(n_phi);

    libmesh_assert_equal_to(Bmat.m(), 3);
    libmesh_assert_equal_to(Bmat.n(), 6*n_phi);

    // now set the shape function values
    // dN/dx
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
        phi(i_nd) = dphi[i_nd][qp](0);
    Bmat.set_shape_function(0, 0, phi); //  epsilon_xx = du/dx
    Bmat.set_shape_function(2, 1, phi); //  gamma_xy = dv/dx + ...
    
    // dN/dy
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
        phi(i_nd) = dphi[i_nd][qp](1);
    Bmat.set_shape_function(1, 1, phi); //  epsilon_yy = dv/dy
    Bmat.set_shape_function(2, 0, phi); //  gamma_xy = du/dy + ...
}



void
MAST::StructuralElement2D::initialize_von_karman_strain_operator(const unsigned int qp,
                                                                 DenseRealVector& vk_strain,
                                                                 DenseRealMatrix& vk_dwdxi_mat,
                                                                 FEMOperatorMatrix& Bmat_vk) {
    
    const std::vector<std::vector<libMesh::RealVectorValue> >& dphi = _fe->get_dphi();
    const unsigned int n_phi = (unsigned int)dphi.size();
    
    libmesh_assert_equal_to(vk_strain.size(), 3);
    libmesh_assert_equal_to(vk_dwdxi_mat.m(), 3);
    libmesh_assert_equal_to(vk_dwdxi_mat.n(), 2);
    libmesh_assert_equal_to(Bmat_vk.m(), 2);
    libmesh_assert_equal_to(Bmat_vk.n(), 6*n_phi);

    Real dw=0.;
    vk_strain.zero();
    vk_dwdxi_mat.zero();
    
    DenseRealVector phi_vec; phi_vec.resize(n_phi);
    
    dw = 0.;
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ ) {
        phi_vec(i_nd) = dphi[i_nd][qp](0);  // dphi/dx
        dw += phi_vec(i_nd)*local_solution(2*n_phi+i_nd); // dw/dx
    }
    Bmat_vk.set_shape_function(0, 2, phi_vec); // dw/dx
    vk_dwdxi_mat(0, 0) = dw;  // epsilon-xx : dw/dx
    vk_dwdxi_mat(2, 1) = dw;  // gamma-xy : dw/dx
    vk_strain(0) = 0.5*dw*dw; // 1/2 * (dw/dx)^2
    vk_strain(2) = dw;        // (dw/dx)*(dw/dy)  only dw/dx is provided here
    
    dw = 0.;
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ ) {
        phi_vec(i_nd) = dphi[i_nd][qp](1);  // dphi/dy
        dw += phi_vec(i_nd)*local_solution(2*n_phi+i_nd); // dw/dy
    }
    Bmat_vk.set_shape_function(1, 2, phi_vec); // dw/dy
    vk_dwdxi_mat(1, 1) = dw;  // epsilon-yy : dw/dy
    vk_dwdxi_mat(2, 0) = dw;  // gamma-xy : dw/dy
    vk_strain(1) = 0.5*dw*dw; // 1/2 * (dw/dy)^2
    vk_strain(2) *= dw;       // (dw/dx)*(dw/dy)
}




void
MAST::StructuralElement2D::initialize_von_karman_strain_operator_sensitivity
(const unsigned int qp,
 DenseRealMatrix &vk_dwdxi_mat_sens) {
    const std::vector<std::vector<libMesh::RealVectorValue> >& dphi = _fe->get_dphi();
    const unsigned int n_phi = (unsigned int)dphi.size();
    
    libmesh_assert_equal_to(vk_dwdxi_mat_sens.m(), 3);
    libmesh_assert_equal_to(vk_dwdxi_mat_sens.n(), 2);
    
    Real dw=0.;
    vk_dwdxi_mat_sens.zero();
    
    DenseRealVector phi_vec; phi_vec.resize(n_phi);
    
    dw = 0.;
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ ) {
        phi_vec(i_nd) = dphi[i_nd][qp](0);  // dphi/dx
        dw += phi_vec(i_nd)*local_solution_sens(2*n_phi+i_nd); // dw/dx
    }
    vk_dwdxi_mat_sens(0, 0) = dw;  // epsilon-xx : dw/dx
    vk_dwdxi_mat_sens(2, 1) = dw;  // gamma-xy : dw/dx
    
    dw = 0.;
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ ) {
        phi_vec(i_nd) = dphi[i_nd][qp](1);  // dphi/dy
        dw += phi_vec(i_nd)*local_solution_sens(2*n_phi+i_nd); // dw/dy
    }
    vk_dwdxi_mat_sens(1, 1) = dw;  // epsilon-yy : dw/dy
    vk_dwdxi_mat_sens(2, 0) = dw;  // gamma-xy : dw/dy
}




bool
MAST::StructuralElement2D::internal_force (bool request_jacobian,
                                             DenseRealVector& f,
                                             DenseRealMatrix& jac,
                                             bool if_ignore_ho_jac)
{
    FEMOperatorMatrix Bmat_mem, Bmat_bend, Bmat_vk;
    
    const std::vector<Real>& JxW = _fe->get_JxW();
    const std::vector<libMesh::Point>& xyz = _fe->get_xyz();
    const unsigned int n_phi = (unsigned int)_fe->get_phi().size();
    const unsigned int n1= this->n_direct_strain_components(), n2=6*n_phi,
    n3 = this->n_von_karman_strain_components();
    DenseRealMatrix material_A_mat, material_B_mat, material_D_mat,
    mat1_n1n2, mat2_n2n2, mat3,
    mat4_n3n2, vk_dwdxi_mat, stress, stress_l, local_jac;
    DenseRealVector  vec1_n1, vec2_n1, vec3_n2,
    vec4_n3, vec5_n3, local_f;
    
    mat1_n1n2.resize(n1, n2); mat2_n2n2.resize(n2, n2);
    mat4_n3n2.resize(n3, n2); local_jac.resize(n2, n2);
    vk_dwdxi_mat.resize(n1,n3); stress.resize(2,2); stress_l.resize(2, 2);
    local_f.resize(n2); vec1_n1.resize(n1); vec2_n1.resize(n1);
    vec3_n2.resize(n2); vec4_n3.resize(n3); vec5_n3.resize(n3);
    
    
    Bmat_mem.reinit(n1, _system.n_vars(), n_phi); // three stress-strain components
    Bmat_bend.reinit(n1, _system.n_vars(), n_phi);
    Bmat_vk.reinit(n3, _system.n_vars(), n_phi); // only dw/dx and dw/dy
    
    bool if_vk = (_property.strain_type() == MAST::VON_KARMAN_STRAIN),
    if_bending = (_property.bending_model(_elem, _fe->get_fe_type()) != MAST::NO_BENDING);
    
    std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > >
    mat_stiff_A(_property.get_property
                (MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX,
                 *this).release()),
    mat_stiff_B(_property.get_property
                (MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_B_MATRIX,
                 *this).release()),
    mat_stiff_D(_property.get_property
                (MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_D_MATRIX,
                 *this).release());
    
    
    libMesh::Point p;
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        this->global_coordinates(xyz[qp], p);
        
        // get the material matrix
        (*mat_stiff_A)(p, _system.time, material_A_mat);
        
        if (if_bending) {
            (*mat_stiff_B)(p, _system.time, material_B_mat);
            (*mat_stiff_D)(p, _system.time, material_D_mat);
        }
        
        // now calculte the quantity for these matrices
        _internal_force_operation(if_bending, if_vk, n2, qp, JxW,
                                  request_jacobian, if_ignore_ho_jac,
                                  local_f, local_jac,
                                  Bmat_mem, Bmat_bend, Bmat_vk,
                                  stress, stress_l, vk_dwdxi_mat, material_A_mat,
                                  material_B_mat, material_D_mat, vec1_n1,
                                  vec2_n1, vec3_n2, vec4_n3,
                                  vec5_n3, mat1_n1n2, mat2_n2n2,
                                  mat3, mat4_n3n2);
        
    }
    
    
    // now calculate the transverse shear contribution if appropriate for the
    // element
    if (if_bending && _bending_operator->include_transverse_shear_energy())
        _bending_operator->calculate_transverse_shear_force(request_jacobian,
                                                            local_f, local_jac,
                                                            NULL);
    
    
    // now transform to the global coorodinate system
    transform_to_global_system(local_f, vec3_n2);
    f.add(1., vec3_n2);
    if (request_jacobian) {
        // for 2D elements
        if (_elem.dim() == 2) {
            // add small values to the diagonal of the theta_z dofs
            for (unsigned int i=0; i<n_phi; i++)
                local_jac(5*n_phi+i, 5*n_phi+i) = -1.0e-8;
        }
        transform_to_global_system(local_jac, mat2_n2n2);
        jac.add(1., mat2_n2n2);
    }
    
    return request_jacobian;
}





bool
MAST::StructuralElement2D::internal_force_sensitivity (bool request_jacobian,
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
    
    FEMOperatorMatrix Bmat_mem, Bmat_bend, Bmat_vk;
    
    const std::vector<Real>& JxW = _fe->get_JxW();
    const std::vector<libMesh::Point>& xyz = _fe->get_xyz();
    const unsigned int n_phi = (unsigned int)_fe->get_phi().size();
    const unsigned int n1= this->n_direct_strain_components(), n2=6*n_phi,
    n3 = this->n_von_karman_strain_components();
    DenseRealMatrix material_A_mat, material_B_mat, material_D_mat,
    material_trans_shear_mat, mat1_n1n2, mat2_n2n2, mat3,
    mat4_n3n2, vk_dwdxi_mat, stress, stress_l, local_jac;
    DenseRealVector  vec1_n1, vec2_n1, vec3_n2,
    vec4_n3, vec5_n3, local_f;
    
    mat1_n1n2.resize(n1, n2); mat2_n2n2.resize(n2, n2);
    mat4_n3n2.resize(n3, n2); local_jac.resize(n2, n2);
    vk_dwdxi_mat.resize(n1,n3); stress.resize(2,2); stress_l.resize(2, 2);
    local_f.resize(n2);
    vec1_n1.resize(n1); vec2_n1.resize(n1);
    vec3_n2.resize(n2); vec4_n3.resize(n3); vec5_n3.resize(n3);
    
    
    Bmat_mem.reinit(n1, _system.n_vars(), n_phi); // three stress-strain components
    Bmat_bend.reinit(n1, _system.n_vars(), n_phi);
    Bmat_vk.reinit(n3, _system.n_vars(), n_phi); // only dw/dx and dw/dy
    
    bool if_vk = (_property.strain_type() == MAST::VON_KARMAN_STRAIN),
    if_bending = (_property.bending_model(_elem, _fe->get_fe_type()) != MAST::NO_BENDING);
    
    std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > >
    mat_stiff_A
    (_property.get_property(MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX,
                            *this).release()),
    mat_stiff_B
    (_property.get_property(MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_B_MATRIX,
                            *this).release()),
    mat_stiff_D
    (_property.get_property(MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_D_MATRIX,
                            *this).release());
    
    libMesh::Point p;
    
    // first calculate the sensitivity due to the parameter
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        this->global_coordinates(xyz[qp], p);
        
        // get the material matrix
        mat_stiff_A->total(*this->sensitivity_param,
                           p, _system.time, material_A_mat);
        
        if (if_bending) {
            mat_stiff_B->total(*this->sensitivity_param,
                               p, _system.time, material_B_mat);
            mat_stiff_D->total(*this->sensitivity_param,
                               p, _system.time, material_D_mat);
        }
        
        // now calculte the quantity for these matrices
        // this accounts for the sensitivity of the material property matrices
        _internal_force_operation(if_bending, if_vk, n2, qp, JxW,
                                  request_jacobian, if_ignore_ho_jac,
                                  local_f, local_jac,
                                  Bmat_mem, Bmat_bend, Bmat_vk,
                                  stress, stress_l, vk_dwdxi_mat, material_A_mat,
                                  material_B_mat, material_D_mat, vec1_n1,
                                  vec2_n1, vec3_n2, vec4_n3,
                                  vec5_n3, mat1_n1n2, mat2_n2n2,
                                  mat3, mat4_n3n2);
        
        // this accounts for the sensitivity of the linear stress as a result of
        // static solution. This is needed only for cases that require linearized
        // geometric stiffness matrix, for example in buckling or natural frequency
        // analysis
        if (if_bending && if_vk && if_ignore_ho_jac && request_jacobian) {
            (*mat_stiff_A)(p, _system.time, material_A_mat);
            (*mat_stiff_B)(p, _system.time, material_B_mat);
            
            _linearized_geometric_stiffness_sensitivity_with_static_solution
            (n2, qp, JxW, local_jac, Bmat_mem, Bmat_bend, Bmat_vk, stress_l,
             vk_dwdxi_mat, material_A_mat, material_B_mat, vec1_n1,
             vec2_n1, mat1_n1n2, mat2_n2n2, mat3);
        }
    }
    
    // now transform to the global coorodinate system
    transform_to_global_system(local_f, vec3_n2);
    f.add(1., vec3_n2);
    if (request_jacobian) {
        transform_to_global_system(local_jac, mat2_n2n2);
        jac.add(1., mat2_n2n2);
    }
    
    // now calculate the transverse shear contribution if appropriate for the
    // element
    if (if_bending && _bending_operator->include_transverse_shear_energy())
        _bending_operator->calculate_transverse_shear_force(request_jacobian,
                                                            local_f, local_jac,
                                                            this->sensitivity_param);
    
    return request_jacobian;
}


void
MAST::StructuralElement2D::_internal_force_operation
(bool if_bending,
 bool if_vk,
 const unsigned int n2,
 const unsigned int qp,
 const std::vector<Real>& JxW,
 bool request_jacobian,
 bool if_ignore_ho_jac,
 DenseRealVector& local_f,
 DenseRealMatrix& local_jac,
 FEMOperatorMatrix& Bmat_mem,
 FEMOperatorMatrix& Bmat_bend,
 FEMOperatorMatrix& Bmat_vk,
 DenseRealMatrix& stress,
 DenseRealMatrix& stress_l,
 DenseRealMatrix& vk_dwdxi_mat,
 DenseRealMatrix& material_A_mat,
 DenseRealMatrix& material_B_mat,
 DenseRealMatrix& material_D_mat,
 DenseRealVector& vec1_n1,
 DenseRealVector& vec2_n1,
 DenseRealVector& vec3_n2,
 DenseRealVector& vec4_2,
 DenseRealVector& vec5_2,
 DenseRealMatrix& mat1_n1n2,
 DenseRealMatrix& mat2_n2n2,
 DenseRealMatrix& mat3,
 DenseRealMatrix& mat4_2n2)
{
    this->initialize_direct_strain_operator(qp, Bmat_mem);
    
    // first handle constant throught the thickness stresses: membrane and vonKarman
    Bmat_mem.vector_mult(vec1_n1, local_solution);
    material_A_mat.vector_mult(vec2_n1, vec1_n1); // linear direct stress
    
    // copy the stress values to a matrix
    stress_l(0,0) = vec2_n1(0); // sigma_xx
    stress_l(0,1) = vec2_n1(2); // sigma_xy
    stress_l(1,0) = vec2_n1(2); // sigma_yx
    stress_l(1,1) = vec2_n1(1); // sigma_yy
    
    stress = stress_l;

    // get the bending strain operator
    vec2_n1.zero(); // used to store vk strain, if applicable
    if (if_bending) {
        _bending_operator->initialize_bending_strain_operator(qp, Bmat_bend);
        
        Bmat_bend.vector_mult(vec2_n1, local_solution);
        material_B_mat.vector_mult(vec1_n1, vec2_n1);
        stress_l(0,0) += vec2_n1(0); // sigma_xx
        stress_l(0,1) += vec2_n1(2); // sigma_xy
        stress_l(1,0) += vec2_n1(2); // sigma_yx
        stress_l(1,1) += vec2_n1(1); // sigma_yy
        
        stress(0,0) += vec2_n1(0); // sigma_xx
        stress(0,1) += vec2_n1(2); // sigma_xy
        stress(1,0) += vec2_n1(2); // sigma_yx
        stress(1,1) += vec2_n1(1); // sigma_yy
        
        
        if (if_vk)  { // get the vonKarman strain operator if needed
            this->initialize_von_karman_strain_operator(qp,
                                                        vec2_n1, // epsilon_vk
                                                        vk_dwdxi_mat,
                                                        Bmat_vk);

            material_A_mat.vector_mult(vec1_n1, vec2_n1); // stress
            stress(0,0) += vec1_n1(0); // sigma_xx
            stress(0,1) += vec1_n1(2); // sigma_xy
            stress(1,0) += vec1_n1(2); // sigma_yx
            stress(1,1) += vec1_n1(1); // sigma_yy

        }
    }
    
    // add the linear and nonlinear direct strains
    Bmat_mem.vector_mult(vec1_n1, local_solution);
    vec2_n1.add(1., vec1_n1);  // epsilon_mem + epsilon_vk
    
    // copy the total integrated stress to the vector
    vec1_n1(0) = stress(0,0);
    vec1_n1(1) = stress(1,1);
    vec1_n1(2) = stress(0,1);
    
    // now the internal force vector
    // this includes the membrane strain operator with all A and B material operators
    Bmat_mem.vector_mult_transpose(vec3_n2, vec1_n1);
    local_f.add(-JxW[qp], vec3_n2);
    
    if (if_bending) {
        if (if_vk) {
            // von Karman strain
            vk_dwdxi_mat.vector_mult_transpose(vec4_2, vec1_n1);
            Bmat_vk.vector_mult_transpose(vec3_n2, vec4_2);
            local_f.add(-JxW[qp], vec3_n2);
        }
        
        // now coupling with the bending strain
        // B_bend^T [B] B_mem
        material_B_mat.vector_mult(vec1_n1, vec2_n1);
        Bmat_bend.vector_mult_transpose(vec3_n2, vec1_n1);
        local_f.add(-JxW[qp], vec3_n2);
        
        // now bending stress
        Bmat_bend.vector_mult(vec2_n1, local_solution);
        material_D_mat.vector_mult(vec1_n1, vec2_n1);
        Bmat_bend.vector_mult_transpose(vec3_n2, vec1_n1);
        local_f.add(-JxW[qp], vec3_n2);
    }
    
    if (request_jacobian) {
        // membrane - membrane
        Bmat_mem.left_multiply(mat1_n1n2, material_A_mat);
        Bmat_mem.right_multiply_transpose(mat2_n2n2, mat1_n1n2);
        local_jac.add(-JxW[qp], mat2_n2n2);
        
        if (if_bending) {
            if (if_vk) {
                // membrane - vk
                mat3.resize(vk_dwdxi_mat.m(), n2);
                Bmat_vk.left_multiply(mat3, vk_dwdxi_mat);
                mat3.left_multiply(material_A_mat);
                Bmat_mem.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac.add(-JxW[qp], mat2_n2n2);
                
                // vk - membrane
                Bmat_mem.left_multiply(mat1_n1n2, material_A_mat);
                vk_dwdxi_mat.get_transpose(mat3);
                mat3.right_multiply(mat1_n1n2);
                Bmat_vk.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac.add(-JxW[qp], mat2_n2n2);
                
                // if only the first order term of the Jacobian is needed, for
                // example for linearized buckling analysis, then the linear
                // stress combined with the variation of the von Karman strain
                // is included. Otherwise, all terms are included
                if (if_ignore_ho_jac) {
                    // vk - vk: first order term
                    mat3.resize(2, n2);
                    Bmat_vk.left_multiply(mat3, stress_l);
                    Bmat_vk.right_multiply_transpose(mat2_n2n2, mat3);
                    local_jac.add(-JxW[qp], mat2_n2n2);
                }
                else {
                    // vk - vk
                    mat3.resize(2, n2);
                    Bmat_vk.left_multiply(mat3, stress);
                    Bmat_vk.right_multiply_transpose(mat2_n2n2, mat3);
                    local_jac.add(-JxW[qp], mat2_n2n2);
                    
                    mat3.resize(vk_dwdxi_mat.m(), n2);
                    Bmat_vk.left_multiply(mat3, vk_dwdxi_mat);
                    mat3.left_multiply(material_A_mat);
                    mat3.left_multiply_transpose(vk_dwdxi_mat);
                    Bmat_vk.right_multiply_transpose(mat2_n2n2, mat3);
                    local_jac.add(-JxW[qp], mat2_n2n2);
                }
                
                // bending - vk
                mat3.resize(vk_dwdxi_mat.m(), n2);
                Bmat_vk.left_multiply(mat3, vk_dwdxi_mat);
                mat3.left_multiply_transpose(material_B_mat);
                Bmat_bend.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac.add(-JxW[qp], mat2_n2n2);
                
                // vk - bending
                Bmat_bend.left_multiply(mat1_n1n2, material_B_mat);
                vk_dwdxi_mat.get_transpose(mat3);
                mat3.right_multiply(mat1_n1n2);
                Bmat_vk.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac.add(-JxW[qp], mat2_n2n2);
            }
            
            // bending - membrane
            Bmat_mem.left_multiply(mat1_n1n2, material_B_mat);
            Bmat_bend.right_multiply_transpose(mat2_n2n2, mat1_n1n2);
            local_jac.add(-JxW[qp], mat2_n2n2);
            
            // membrane - bending
            Bmat_bend.left_multiply(mat1_n1n2, material_B_mat);
            Bmat_mem.right_multiply_transpose(mat2_n2n2, mat1_n1n2);
            local_jac.add(-JxW[qp], mat2_n2n2);
            
            // bending - bending
            Bmat_bend.left_multiply(mat1_n1n2, material_D_mat);
            Bmat_bend.right_multiply_transpose(mat2_n2n2, mat1_n1n2);
            local_jac.add(-JxW[qp], mat2_n2n2);
        }
    }
}




void
MAST::StructuralElement2D::_linearized_geometric_stiffness_sensitivity_with_static_solution
(const unsigned int n2,
 const unsigned int qp,
 const std::vector<Real>& JxW,
 DenseRealMatrix& local_jac,
 FEMOperatorMatrix& Bmat_mem,
 FEMOperatorMatrix& Bmat_bend,
 FEMOperatorMatrix& Bmat_vk,
 DenseRealMatrix& stress_l,
 DenseRealMatrix& vk_dwdxi_mat,
 DenseRealMatrix& material_A_mat,
 DenseRealMatrix& material_B_mat,
 DenseRealVector& vec1_n1,
 DenseRealVector& vec2_n1,
 DenseRealMatrix& mat1_n1n2,
 DenseRealMatrix& mat2_n2n2,
 DenseRealMatrix& mat3)
{
    this->initialize_direct_strain_operator(qp, Bmat_mem);
    _bending_operator->initialize_bending_strain_operator(qp, Bmat_bend);

    // first handle constant throught the thickness stresses: membrane and vonKarman
    Bmat_mem.vector_mult(vec1_n1, local_solution_sens);
    material_A_mat.vector_mult(vec2_n1, vec1_n1); // linear direct stress
    
    // copy the stress values to a matrix
    stress_l(0,0) = vec2_n1(0); // sigma_xx
    stress_l(0,1) = vec2_n1(2); // sigma_xy
    stress_l(1,0) = vec2_n1(2); // sigma_yx
    stress_l(1,1) = vec2_n1(1); // sigma_yy

    // get the von Karman operator matrix
    this->initialize_von_karman_strain_operator(qp,
                                                vec2_n1, // epsilon_vk
                                                vk_dwdxi_mat,
                                                Bmat_vk);

    // sensitivity of the vk_dwdxi matrix due to solution sensitivity
    this->initialize_von_karman_strain_operator_sensitivity(qp,
                                                            vk_dwdxi_mat);

    // membrane - vk
    mat3.resize(vk_dwdxi_mat.m(), n2);
    Bmat_vk.left_multiply(mat3, vk_dwdxi_mat);
    mat3.left_multiply(material_A_mat);
    Bmat_mem.right_multiply_transpose(mat2_n2n2, mat3);
    local_jac.add(-JxW[qp], mat2_n2n2);
    
    // vk - membrane
    Bmat_mem.left_multiply(mat1_n1n2, material_A_mat);
    vk_dwdxi_mat.get_transpose(mat3);
    mat3.right_multiply(mat1_n1n2);
    Bmat_vk.right_multiply_transpose(mat2_n2n2, mat3);
    local_jac.add(-JxW[qp], mat2_n2n2);
    
    // vk - vk: first order term
    mat3.resize(2, n2);
    Bmat_vk.left_multiply(mat3, stress_l);
    Bmat_vk.right_multiply_transpose(mat2_n2n2, mat3);
    local_jac.add(-JxW[qp], mat2_n2n2);

    // bending - vk
    mat3.resize(vk_dwdxi_mat.m(), n2);
    Bmat_vk.left_multiply(mat3, vk_dwdxi_mat);
    mat3.left_multiply_transpose(material_B_mat);
    Bmat_bend.right_multiply_transpose(mat2_n2n2, mat3);
    local_jac.add(-JxW[qp], mat2_n2n2);
    
    // vk - bending
    Bmat_bend.left_multiply(mat1_n1n2, material_B_mat);
    vk_dwdxi_mat.get_transpose(mat3);
    mat3.right_multiply(mat1_n1n2);
    Bmat_vk.right_multiply_transpose(mat2_n2n2, mat3);
    local_jac.add(-JxW[qp], mat2_n2n2);
}




bool
MAST::StructuralElement2D::prestress_force (bool request_jacobian,
                                              DenseRealVector& f,
                                              DenseRealMatrix& jac)
{
    if (!_property.if_prestressed())
        return false;
    
    FEMOperatorMatrix Bmat_mem, Bmat_bend, Bmat_vk;
    
    const std::vector<Real>& JxW = _fe->get_JxW();
    const std::vector<libMesh::Point>& xyz = _fe->get_xyz();
    const unsigned int n_phi = (unsigned int)_fe->get_phi().size();
    const unsigned int n1= this->n_direct_strain_components(), n2=6*n_phi,
    n3 = this->n_von_karman_strain_components();
    DenseRealMatrix mat2_n2n2, mat3, vk_dwdxi_mat, local_jac,
    prestress_mat_A, prestress_mat_B;
    DenseRealVector vec2_n1, vec3_n2, vec4_n3, vec5_n3,
    local_f, prestress_vec_A, prestress_vec_B;
    
    mat2_n2n2.resize(n2, n2); local_jac.resize(n2, n2);
    vk_dwdxi_mat.resize(n1,n3); local_f.resize(n2); vec2_n1.resize(n1);
    vec3_n2.resize(n2); vec4_n3.resize(n3); vec5_n3.resize(n3);
    
    
    Bmat_mem.reinit(n1, _system.n_vars(), n_phi); // three stress-strain components
    Bmat_bend.reinit(n1, _system.n_vars(), n_phi);
    Bmat_vk.reinit(n3, _system.n_vars(), n_phi); // only dw/dx and dw/dy
    
    bool if_vk = (_property.strain_type() == MAST::VON_KARMAN_STRAIN),
    if_bending = (_property.bending_model(_elem, _fe->get_fe_type()) != MAST::NO_BENDING);
    
    std::auto_ptr<MAST::SectionIntegratedPrestressMatrixBase>
    prestress_A
    (dynamic_cast<MAST::SectionIntegratedPrestressMatrixBase*>
     (_property.get_property(MAST::SECTION_INTEGRATED_PRESTRESS_A_MATRIX,
                             *this).release())),
    prestress_B
    (dynamic_cast<MAST::SectionIntegratedPrestressMatrixBase*>
     (_property.get_property(MAST::SECTION_INTEGRATED_PRESTRESS_B_MATRIX,
                             *this).release()));
    
    libMesh::Point p;
    
    // now calculate the quantity
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        this->global_coordinates(xyz[qp], p);
        
        (*prestress_A)(p, _system.time, prestress_mat_A);
        prestress_A->convert_to_vector(prestress_mat_A, prestress_vec_A);
        (*prestress_B)(p, _system.time, prestress_mat_B);
        prestress_B->convert_to_vector(prestress_mat_B, prestress_vec_B);
        
        this->initialize_direct_strain_operator(qp, Bmat_mem);
        
        // get the bending strain operator if needed
        vec2_n1.zero(); // used to store vk strain, if applicable
        if (if_bending) {
            _bending_operator->initialize_bending_strain_operator(qp, Bmat_bend);
            
            if (if_vk)  // get the vonKarman strain operator if needed
                this->initialize_von_karman_strain_operator(qp,
                                                            vec2_n1,
                                                            vk_dwdxi_mat,
                                                            Bmat_vk);
        }
        
        // first handle constant throught the thickness stresses: membrane and vonKarman
        // multiply this with the constant through the thickness strain
        // membrane strain
        Bmat_mem.vector_mult_transpose(vec3_n2, prestress_vec_A);
        local_f.add(-JxW[qp], vec3_n2); // epsilon_mem * sigma_0
        
        if (if_bending) {
            if (if_vk) {
                // von Karman strain
                vk_dwdxi_mat.vector_mult_transpose(vec4_n3, prestress_vec_A);
                Bmat_vk.vector_mult_transpose(vec3_n2, vec4_n3);
                local_f.add(-JxW[qp], vec3_n2); // epsilon_vk * sigma_0
            }
            
            // now coupling with the bending strain
            Bmat_bend.vector_mult_transpose(vec3_n2, prestress_vec_B);
            local_f.add(-JxW[qp], vec3_n2); // epsilon_bend * sigma_0
        }
        
        if (request_jacobian) {
            if (if_bending && if_vk) {
                mat3.resize(2, n2);
                Bmat_vk.left_multiply(mat3, prestress_mat_A);
                Bmat_vk.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac.add(-JxW[qp], mat2_n2n2);
            }
        }
    }
    
    // now transform to the global coorodinate system
    transform_to_global_system(local_f, vec3_n2);
    f.add(1., vec3_n2);
    if (request_jacobian && if_vk) {
        transform_to_global_system(local_jac, mat2_n2n2);
        jac.add(1., mat2_n2n2);
    }
    
    // only the nonlinear strain returns a Jacobian for prestressing
    return (request_jacobian && if_vk);
}





bool
MAST::StructuralElement2D::prestress_force_sensitivity (bool request_jacobian,
                                                          DenseRealVector& f,
                                                          DenseRealMatrix& jac)
{
    if (!_property.if_prestressed())
        return false;
    
    FEMOperatorMatrix Bmat_mem, Bmat_bend, Bmat_vk;
    
    const std::vector<Real>& JxW = _fe->get_JxW();
    const std::vector<libMesh::Point>& xyz = _fe->get_xyz();
    const unsigned int n_phi = (unsigned int)_fe->get_phi().size();
    const unsigned int n1= this->n_direct_strain_components(), n2=6*n_phi,
    n3 = this->n_von_karman_strain_components();
    DenseRealMatrix mat2_n2n2, mat3, vk_dwdxi_mat, local_jac,
    prestress_mat_A, prestress_mat_B;
    DenseRealVector vec2_n1, vec3_n2, vec4_n3, vec5_n3,
    local_f, prestress_vec_A, prestress_vec_B;
    
    mat2_n2n2.resize(n2, n2); local_jac.resize(n2, n2);
    vk_dwdxi_mat.resize(n1,n3); local_f.resize(n2); vec2_n1.resize(n1);
    vec3_n2.resize(n2); vec4_n3.resize(n3); vec5_n3.resize(n3);
    
    
    Bmat_mem.reinit(n1, _system.n_vars(), n_phi); // three stress-strain components
    Bmat_bend.reinit(n1, _system.n_vars(), n_phi);
    Bmat_vk.reinit(n3, _system.n_vars(), n_phi); // only dw/dx and dw/dy
    
    bool if_vk = (_property.strain_type() == MAST::VON_KARMAN_STRAIN),
    if_bending = (_property.bending_model(_elem, _fe->get_fe_type()) != MAST::NO_BENDING);
    
    std::auto_ptr<MAST::SectionIntegratedPrestressMatrixBase>
    prestress_A
    (dynamic_cast<MAST::SectionIntegratedPrestressMatrixBase*>
     (_property.get_property(MAST::SECTION_INTEGRATED_PRESTRESS_A_MATRIX,
                             *this).release())),
    prestress_B
    (dynamic_cast<MAST::SectionIntegratedPrestressMatrixBase*>
     (_property.get_property(MAST::SECTION_INTEGRATED_PRESTRESS_B_MATRIX,
                             *this).release()));
    
    libMesh::Point p;
    
    // transform to the local coordinate system
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        this->global_coordinates(xyz[qp], p);
        
        prestress_A->total(*this->sensitivity_param,
                           p, _system.time, prestress_mat_A);
        prestress_A->convert_to_vector(prestress_mat_A, prestress_vec_A);
        prestress_B->total(*this->sensitivity_param,
                           p, _system.time, prestress_mat_B);
        prestress_B->convert_to_vector(prestress_mat_B, prestress_vec_B);
        
        this->initialize_direct_strain_operator(qp, Bmat_mem);
        
        // get the bending strain operator if needed
        vec2_n1.zero(); // used to store vk strain, if applicable
        if (if_bending) {
            _bending_operator->initialize_bending_strain_operator(qp, Bmat_bend);
            
            if (if_vk)  // get the vonKarman strain operator if needed
                this->initialize_von_karman_strain_operator(qp,
                                                            vec2_n1,
                                                            vk_dwdxi_mat,
                                                            Bmat_vk);
        }
        
        // first handle constant throught the thickness stresses: membrane and vonKarman
        // multiply this with the constant through the thickness strain
        // membrane strain
        Bmat_mem.vector_mult_transpose(vec3_n2, prestress_vec_A);
        local_f.add(-JxW[qp], vec3_n2); // epsilon_mem * sigma_0
        
        if (if_bending) {
            if (if_vk) {
                // von Karman strain
                vk_dwdxi_mat.vector_mult_transpose(vec4_n3, prestress_vec_A);
                Bmat_vk.vector_mult_transpose(vec3_n2, vec4_n3);
                local_f.add(-JxW[qp], vec3_n2); // epsilon_vk * sigma_0
            }
            
            // now coupling with the bending strain
            Bmat_bend.vector_mult_transpose(vec3_n2, prestress_vec_B);
            local_f.add(-JxW[qp], vec3_n2); // epsilon_bend * sigma_0
        }
        
        if (request_jacobian) {
            if (if_bending && if_vk) {
                mat3.resize(2, n2);
                Bmat_vk.left_multiply(mat3, prestress_mat_A);
                Bmat_vk.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac.add(-JxW[qp], mat2_n2n2);
            }
        }
    }
    
    // now transform to the global coorodinate system
    transform_to_global_system(local_f, vec3_n2);
    f.add(1., vec3_n2);
    if (request_jacobian && if_vk) {
        transform_to_global_system(local_jac, mat2_n2n2);
        jac.add(1., mat2_n2n2);
    }
    
    // only the nonlinear strain returns a Jacobian for prestressing
    return (request_jacobian && if_vk);
}





bool
MAST::StructuralElement2D::thermal_force (bool request_jacobian,
                                            DenseRealVector& f,
                                            DenseRealMatrix& jac,
                                            MAST::BoundaryCondition& p)
{
    FEMOperatorMatrix Bmat_mem, Bmat_bend, Bmat_vk;
    
    const std::vector<Real>& JxW = _fe->get_JxW();
    const std::vector<libMesh::Point>& xyz = _fe->get_xyz();
    const unsigned int n_phi = (unsigned int)_fe->get_phi().size();
    const unsigned int n1= this->n_direct_strain_components(), n2=6*n_phi,
    n3 = this->n_von_karman_strain_components();
    DenseRealMatrix material_exp_A_mat, material_exp_B_mat,
    mat1_n1n2, mat2_n2n2, mat3,
    mat4_n3n2, vk_dwdxi_mat, stress, local_jac;
    DenseRealVector  vec1_n1, vec2_n1, vec3_n2,
    vec4_2, vec5_n3, local_f, delta_t;
    
    mat1_n1n2.resize(n1, n2); mat2_n2n2.resize(n2, n2);
    mat4_n3n2.resize(n3, n2); local_jac.resize(n2, n2);
    vk_dwdxi_mat.resize(n1,n3); stress.resize(2,2); local_f.resize(n2);
    vec1_n1.resize(n1); vec2_n1.resize(n1);
    vec3_n2.resize(n2); vec4_2.resize(2); vec5_n3.resize(n3);
    delta_t.resize(1);
    
    
    Bmat_mem.reinit(n1, _system.n_vars(), n_phi); // three stress-strain components
    Bmat_bend.reinit(n1, _system.n_vars(), n_phi);
    Bmat_vk.reinit(n3, _system.n_vars(), n_phi); // only dw/dx and dw/dy
    
    bool if_vk = (_property.strain_type() == MAST::VON_KARMAN_STRAIN),
    if_bending = (_property.bending_model(_elem, _fe->get_fe_type()) != MAST::NO_BENDING);
    
    std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > > expansion_A
    (_property.get_property(MAST::SECTION_INTEGRATED_MATERIAL_THERMAL_EXPANSION_A_MATRIX,
                            *this).release()),
    expansion_B
    (_property.get_property(MAST::SECTION_INTEGRATED_MATERIAL_THERMAL_EXPANSION_B_MATRIX,
                            *this).release());
    
    // temperature function
    MAST::FieldFunction<Real>& temp_func =
    dynamic_cast<MAST::FieldFunction<Real>&>(p.function());
    MAST::FieldFunction<Real>& ref_temp_func =
    dynamic_cast<MAST::FieldFunction<Real>&>
    (dynamic_cast<MAST::Temperature&>(p).reference_temperature_function());
    
    Real t, t0;
    libMesh::Point pt;
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        this->global_coordinates(xyz[qp], pt);
        
        // this is moved inside the domain since
        (*expansion_A)(pt, _system.time, material_exp_A_mat);
        (*expansion_B)(pt, _system.time, material_exp_B_mat);

        // get the temperature function
        temp_func(xyz[qp], _system.time, t);
        ref_temp_func(xyz[qp], _system.time, t0);
        delta_t(0) = t-t0;

        material_exp_A_mat.vector_mult(vec1_n1, delta_t); // [C]{alpha (T - T0)} (with membrane strain)
        material_exp_B_mat.vector_mult(vec2_n1, delta_t); // [C]{alpha (T - T0)} (with bending strain)
        stress(0,0) = vec1_n1(0); // sigma_xx
        stress(0,1) = vec1_n1(2); // sigma_xy
        stress(1,0) = vec1_n1(2); // sigma_yx
        stress(1,1) = vec1_n1(1); // sigma_yy
        
        this->initialize_direct_strain_operator(qp, Bmat_mem);
        
        // membrane strain
        Bmat_mem.vector_mult_transpose(vec3_n2, vec1_n1);
        local_f.add(JxW[qp], vec3_n2);
        
        if (if_bending) {
            // bending strain
            _bending_operator->initialize_bending_strain_operator(qp, Bmat_bend);
            Bmat_bend.vector_mult_transpose(vec3_n2, vec2_n1);
            local_f.add(JxW[qp], vec3_n2);
            
            // von Karman strain
            if (if_vk) {
                // get the vonKarman strain operator if needed
                this->initialize_von_karman_strain_operator(qp,
                                                            vec2_n1, // epsilon_vk
                                                            vk_dwdxi_mat,
                                                            Bmat_vk);
                // von Karman strain
                vk_dwdxi_mat.vector_mult_transpose(vec4_2, vec1_n1);
                Bmat_vk.vector_mult_transpose(vec3_n2, vec4_2);
                local_f.add(JxW[qp], vec3_n2);
            }
            
            if (request_jacobian && if_vk) { // Jacobian only for vk strain
                // vk - vk
                mat3.resize(2, n2);
                Bmat_vk.left_multiply(mat3, stress);
                Bmat_vk.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac.add(JxW[qp], mat2_n2n2);
            }
        }
    }
    
    // now transform to the global coorodinate system
    transform_to_global_system(local_f, vec3_n2);
    f.add(1., vec3_n2);
    if (request_jacobian && if_vk) {
        transform_to_global_system(local_jac, mat2_n2n2);
        jac.add(1., mat2_n2n2);
    }
    
    // Jacobian contribution from von Karman strain
    return request_jacobian && if_vk;
}




bool
MAST::StructuralElement2D::thermal_force_sensitivity (bool request_jacobian,
                                                        DenseRealVector& f,
                                                        DenseRealMatrix& jac,
                                                        MAST::BoundaryCondition& p)
{
    FEMOperatorMatrix Bmat_mem, Bmat_bend, Bmat_vk;
    
    const std::vector<Real>& JxW = _fe->get_JxW();
    const std::vector<libMesh::Point>& xyz = _fe->get_xyz();
    const unsigned int n_phi = (unsigned int)_fe->get_phi().size();
    const unsigned int n1= this->n_direct_strain_components(), n2=6*n_phi,
    n3 = this->n_von_karman_strain_components();
    DenseRealMatrix material_exp_A_mat, material_exp_B_mat,
    material_exp_A_mat_sens, material_exp_B_mat_sens,
    mat1_n1n2, mat2_n2n2, mat3,
    mat4_n3n2, vk_dwdxi_mat, stress, local_jac;
    DenseRealVector  vec1_n1, vec2_n1, vec3_n2,
    vec4_2, vec5_n1, local_f, delta_t, delta_t_sens;
    
    mat1_n1n2.resize(n1, n2); mat2_n2n2.resize(n2, n2);
    mat4_n3n2.resize(n3, n2); local_jac.resize(n2, n2);
    vk_dwdxi_mat.resize(n1,n3); stress.resize(2,2); local_f.resize(n2);
    vec1_n1.resize(n1); vec2_n1.resize(n1);
    vec3_n2.resize(n2); vec4_2.resize(2); vec5_n1.resize(n1);
    delta_t.resize(1); delta_t_sens.resize(1);
    
    
    Bmat_mem.reinit(n1, _system.n_vars(), n_phi); // three stress-strain components
    Bmat_bend.reinit(n1, _system.n_vars(), n_phi);
    Bmat_vk.reinit(n3, _system.n_vars(), n_phi); // only dw/dx and dw/dy
    
    bool if_vk = (_property.strain_type() == MAST::VON_KARMAN_STRAIN),
    if_bending = (_property.bending_model(_elem, _fe->get_fe_type()) != MAST::NO_BENDING);
    
    std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > > expansion_A
    (_property.get_property(MAST::SECTION_INTEGRATED_MATERIAL_THERMAL_EXPANSION_A_MATRIX,
                            *this).release()),
    expansion_B
    (_property.get_property(MAST::SECTION_INTEGRATED_MATERIAL_THERMAL_EXPANSION_B_MATRIX,
                            *this).release());
    
    // temperature function
    MAST::FieldFunction<Real>& temp_func =
    dynamic_cast<MAST::FieldFunction<Real>&>(p.function());
    MAST::FieldFunction<Real>& ref_temp_func =
    dynamic_cast<MAST::FieldFunction<Real>&>
    (dynamic_cast<MAST::Temperature&>(p).reference_temperature_function());
    
    Real t, t0, t_sens;
    libMesh::Point pt;
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        this->global_coordinates(xyz[qp], pt);
        
        // this is moved inside the domain since
        (*expansion_A)(pt, _system.time, material_exp_A_mat);
        expansion_A->total(*this->sensitivity_param,
                           pt, _system.time, material_exp_A_mat_sens);
        (*expansion_B)(pt, _system.time, material_exp_B_mat);
        expansion_B->total(*this->sensitivity_param,
                           pt, _system.time, material_exp_B_mat_sens);
        
        // get the temperature function
        temp_func(xyz[qp], _system.time, t);
        ref_temp_func(xyz[qp], _system.time, t0);
        delta_t(0) = t-t0;
        
        // get the temperature function
        temp_func(xyz[qp], _system.time, t);
        temp_func.total(*this->sensitivity_param,
                        xyz[qp], _system.time, t_sens);
        ref_temp_func(xyz[qp], _system.time, t0);
        delta_t(0) = t-t0;
        delta_t_sens(0) = t_sens;

        // now prepare the membrane force sensitivity
        material_exp_A_mat.vector_mult(vec1_n1, delta_t_sens); // [C]{alpha dT/dp} (with membrane strain)
        material_exp_A_mat_sens.vector_mult(vec2_n1, delta_t); // d([C]alpha)/dp (T - T0)} (with membrane strain)
        vec1_n1.add(1., vec2_n1);
        stress(0,0) = vec1_n1(0); // sigma_xx
        stress(0,1) = vec1_n1(2); // sigma_xy
        stress(1,0) = vec1_n1(2); // sigma_yx
        stress(1,1) = vec1_n1(1); // sigma_yy
        
        material_exp_B_mat.vector_mult(vec2_n1, delta_t_sens); // [C]{alpha dT/dp} (with bending strain)
        material_exp_B_mat_sens.vector_mult(vec5_n1, delta_t); // d([C] alpha)/dp (T - T0) (with bending strain)
        vec2_n1.add(1., vec5_n1);
        
        
        this->initialize_direct_strain_operator(qp, Bmat_mem);
        
        // membrane strain
        Bmat_mem.vector_mult_transpose(vec3_n2, vec1_n1);
        local_f.add(JxW[qp], vec3_n2);
        
        if (if_bending) {
            // bending strain
            _bending_operator->initialize_bending_strain_operator(qp, Bmat_bend);
            Bmat_bend.vector_mult_transpose(vec3_n2, vec2_n1);
            local_f.add(JxW[qp], vec3_n2);
            
            // von Karman strain
            if (if_vk) {
                // get the vonKarman strain operator if needed
                this->initialize_von_karman_strain_operator(qp,
                                                            vec2_n1, // epsilon_vk
                                                            vk_dwdxi_mat,
                                                            Bmat_vk);
                // von Karman strain
                vk_dwdxi_mat.vector_mult_transpose(vec4_2, vec1_n1);
                Bmat_vk.vector_mult_transpose(vec3_n2, vec4_2);
                local_f.add(JxW[qp], vec3_n2);
            }
            
            if (request_jacobian && if_vk) { // Jacobian only for vk strain
                
                // vk - vk
                mat3.resize(2, n2);
                Bmat_vk.left_multiply(mat3, stress);
                Bmat_vk.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac.add(JxW[qp], mat2_n2n2);
            }
        }
    }
    
    
    // now transform to the global coorodinate system
    transform_to_global_system(local_f, vec3_n2);
    f.add(1., vec3_n2);
    if (request_jacobian && if_vk) {
        transform_to_global_system(local_jac, mat2_n2n2);
        jac.add(1., mat2_n2n2);
    }
    
    // Jacobian contribution from von Karman strain
    return request_jacobian && if_vk;
}








Real
MAST::StructuralElement2D::max_von_mises_stress() {
    Real max_val = 0.;
    MAST::Stress s;
    
    FEMOperatorMatrix Bmat_mem, Bmat_bend, Bmat_vk;
    
    const std::vector<Real>& JxW = _fe->get_JxW();
    const std::vector<libMesh::Point>& xyz = _fe->get_xyz();
    const unsigned int n_phi = (unsigned int)_fe->get_phi().size();
    const unsigned int n1= this->n_direct_strain_components(), n2=6*n_phi,
    n3 = this->n_von_karman_strain_components();
    DenseRealMatrix material_mat, vk_dwdxi_mat;
    DenseRealVector  vec1_n1, vec2_n1, vec3_n1, strain;
    
    vk_dwdxi_mat.resize(n1,n3);
    vec1_n1.resize(n1); vec2_n1.resize(n1);
    vec3_n1.resize(n1); strain.resize(n1);
    
    
    Bmat_mem.reinit(n1, _system.n_vars(), n_phi); // three stress-strain components
    Bmat_bend.reinit(n1, _system.n_vars(), n_phi);
    Bmat_vk.reinit(n3, _system.n_vars(), n_phi); // only dw/dx and dw/dy
    
    bool if_vk = (_property.strain_type() == MAST::VON_KARMAN_STRAIN),
    if_bending = (_property.bending_model(_elem, _fe->get_fe_type()) != MAST::NO_BENDING);
    
    MAST::BendingOperator2D& bending_2d = dynamic_cast<MAST::BendingOperator2D&>(*_bending_operator);
    
    std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > > material
    (_property.get_material().get_property(MAST::MATERIAL_STIFFNESS_MATRIX,
                                           _property,
                                           2).release());
    
    libMesh::Point p;
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        this->global_coordinates(xyz[qp], p);
        
        (*material)(p, 0., material_mat);
        
        this->initialize_direct_strain_operator(qp, Bmat_mem);
        
        // get the bending strain operator if needed
        vec2_n1.zero(); // used to store vk strain, if applicable
        if (if_bending && if_vk)  // get the vonKarman strain operator if needed
            this->initialize_von_karman_strain_operator(qp,
                                                        vec2_n1, // epsilon_vk
                                                        vk_dwdxi_mat,
                                                        Bmat_vk);
    
        
        // first handle constant throught the thickness stresses: membrane and vonKarman
        Bmat_mem.vector_mult(vec1_n1, local_solution);
        vec2_n1.add(1., vec1_n1);  // epsilon_mem + epsilon_vk
        
        if (if_bending) {
            bending_2d.initialize_bending_strain_operator_for_z(qp, 1., Bmat_bend);
            // calculate the strain at z=1, and the strain at z=-1 would be
            // negative of this
            Bmat_bend.vector_mult(vec3_n1, local_solution);
        }
        
        strain = vec2_n1;
        strain += vec3_n1; // extension and bending strain at +z
        
        // multiply this with the constant-through-the-thickness strain
        // membrane strain
        material_mat.vector_mult(vec1_n1, strain); // stress

        // copy the stress values to the stress tensor
        s(0,0) = vec1_n1(0); // sigma_xx
        s(0,1) = vec1_n1(2); // sigma_xy
        s(1,0) = vec1_n1(2); // sigma_yx
        s(1,1) = vec1_n1(1); // sigma_yy

        // store the maximum value
        max_val = std::max(s.von_mises_stress(), max_val);
        
        // now do the same for -z
        strain = vec2_n1;
        strain -= vec3_n1;
        
        // multiply this with the constant-through-the-thickness strain
        // membrane strain
        material_mat.vector_mult(vec1_n1, strain); // stress
        
        // copy the stress values to the stress tensor
        s(0,0) = vec1_n1(0); // sigma_xx
        s(0,1) = vec1_n1(2); // sigma_xy
        s(1,0) = vec1_n1(2); // sigma_yx
        s(1,1) = vec1_n1(1); // sigma_yy
        
        // store the maximum value
        max_val = std::max(s.von_mises_stress(), max_val);
    }
    
    return max_val;
}

