//
//  structural_element_1D.cpp
//  MAST
//
//  Created by Manav Bhatia on 11/14/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//


// MAST includes
#include "StructuralElems/structural_element_1D.h"
#include "PropertyCards/element_property_card_1D.h"
#include "Numerics/fem_operator_matrix.h"
#include "BoundaryConditions/boundary_condition.h"
#include "BoundaryConditions/temperature.h"


void
MAST::Local1DElem::_create_local_elem() {
    
    libmesh_assert(_elem.dim() == 1);
    
    // first node is the origin of the new cs
    // calculate the coordinate system for the plane of the element
    libMesh::Point v1, v2, v3, p;
    v1 = *_elem.get_node(1); v1 -= *_elem.get_node(0); v1 /= v1.size(); // local x
    v2 = _local_y;                           // vector in local x-y plane
    v3 = v1.cross(v2);                       // local z
    libmesh_assert_greater(v3.size(), 0.);   // 0. implies x == y
    v3 /= v3.size();
    v2 = v3.cross(v1); v2 /= v2.size();      // local y
    
    _T_mat.resize(3,3);
    
    _local_elem = Elem::build(_elem.type()).release();
    _local_nodes.resize(_elem.n_nodes());
    for (unsigned int i=0; i<_elem.n_nodes(); i++) {
        _local_nodes[i] = new Node;
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



MAST::StructuralElement1D::StructuralElement1D(libMesh::System& sys,
                                               const libMesh::Elem& elem,
                                               const MAST::ElementPropertyCardBase& p):
MAST::BendingStructuralElem(sys, elem, p)
{ }






void
MAST::StructuralElement1D::initialize_direct_strain_operator(const unsigned int qp,
                                                             FEMOperatorMatrix& Bmat) {
    
    const std::vector<std::vector<RealVectorValue> >& dphi = _fe->get_dphi();
    
    unsigned int n_phi = (unsigned int)dphi.size();
    DenseRealVector phi; phi.resize(n_phi);
    
    libmesh_assert_equal_to(Bmat.m(), 2);
    libmesh_assert_equal_to(Bmat.n(), 6*n_phi);

    // now set the shape function values
    // dN/dx
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
        phi(i_nd) = dphi[i_nd][qp](0);
    Bmat.set_shape_function(0, 0, phi); //  epsilon_xx = du/dx
    Bmat.set_shape_function(1, 3, phi); //  torsion operator = dtheta_x/dx
}



void
MAST::StructuralElement1D::initialize_von_karman_strain_operator(const unsigned int qp,
                                                                 DenseRealVector& vk_strain,
                                                                 DenseRealMatrix& vk_dvdxi_mat,
                                                                 DenseRealMatrix& vk_dwdxi_mat,
                                                                 FEMOperatorMatrix& Bmat_v_vk,
                                                                 FEMOperatorMatrix& Bmat_w_vk) {

    const std::vector<std::vector<RealVectorValue> >& dphi = _fe->get_dphi();
    const unsigned int n_phi = (unsigned int)dphi.size();

    libmesh_assert_equal_to(vk_strain.size(), 2);
    libmesh_assert_equal_to(vk_dvdxi_mat.m(), 2);
    libmesh_assert_equal_to(vk_dvdxi_mat.n(), 2);
    libmesh_assert_equal_to(Bmat_v_vk.m(), 2);
    libmesh_assert_equal_to(Bmat_v_vk.n(), 6*n_phi);
    libmesh_assert_equal_to(vk_dwdxi_mat.m(), 2);
    libmesh_assert_equal_to(vk_dwdxi_mat.n(), 2);
    libmesh_assert_equal_to(Bmat_w_vk.m(), 2);
    libmesh_assert_equal_to(Bmat_w_vk.n(), 6*n_phi);
    
    libMesh::Real dv=0., dw=0.;
    vk_strain.zero();
    vk_dvdxi_mat.zero();
    vk_dwdxi_mat.zero();
    
    DenseRealVector phi_vec; phi_vec.resize(n_phi);
    
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ ) {
        phi_vec(i_nd) = dphi[i_nd][qp](0);                // dphi/dx
        dv += phi_vec(i_nd)*local_solution(n_phi+i_nd);   // dv/dx
        dw += phi_vec(i_nd)*local_solution(2*n_phi+i_nd); // dw/dx
    }
    
    Bmat_v_vk.set_shape_function(0, 1, phi_vec); // dv/dx
    Bmat_w_vk.set_shape_function(0, 2, phi_vec); // dw/dx
    vk_dvdxi_mat(0, 0) = dv;                   // epsilon-xx : dv/dx
    vk_dwdxi_mat(0, 0) = dw;                   // epsilon-xx : dw/dx
    vk_strain(0) = 0.5*(dv*dv+dw*dw);          // 1/2 * [(dv/dx)^2 + (dw/dx)^2]
}




void
MAST::StructuralElement1D::initialize_von_karman_strain_operator_sensitivity
(const unsigned int qp,
 DenseRealMatrix& vk_dvdxi_mat_sens,
 DenseRealMatrix& vk_dwdxi_mat_sens) {
    
    const std::vector<std::vector<RealVectorValue> >& dphi = _fe->get_dphi();
    const unsigned int n_phi = (unsigned int)dphi.size();
    
    libmesh_assert_equal_to(vk_dvdxi_mat_sens.m(), 2);
    libmesh_assert_equal_to(vk_dvdxi_mat_sens.n(), 2);
    libmesh_assert_equal_to(vk_dwdxi_mat_sens.m(), 2);
    libmesh_assert_equal_to(vk_dwdxi_mat_sens.n(), 2);
    
    libMesh::Real dv=0., dw=0.;
    vk_dvdxi_mat_sens.zero();
    vk_dwdxi_mat_sens.zero();
    
    DenseRealVector phi_vec; phi_vec.resize(n_phi);
    
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ ) {
        phi_vec(i_nd) = dphi[i_nd][qp](0);                // dphi/dx
        dv += phi_vec(i_nd)*local_solution_sens(n_phi+i_nd);   // dv/dx
        dw += phi_vec(i_nd)*local_solution_sens(2*n_phi+i_nd); // dw/dx
    }
    
    vk_dvdxi_mat_sens(0, 0) = dv;                   // epsilon-xx : dv/dx
    vk_dwdxi_mat_sens(0, 0) = dw;                   // epsilon-xx : dw/dx
}



bool
MAST::StructuralElement1D::internal_force (bool request_jacobian,
                                             DenseRealVector& f,
                                             DenseRealMatrix& jac,
                                             bool if_ignore_ho_jac)
{
    FEMOperatorMatrix Bmat_mem, Bmat_bend, Bmat_v_vk, Bmat_w_vk;
    
    const std::vector<libMesh::Real>& JxW = _fe->get_JxW();
    const std::vector<Point>& xyz = _fe->get_xyz();
    const unsigned int n_phi = (unsigned int)_fe->get_phi().size();
    const unsigned int n1= this->n_direct_strain_components(), n2=6*n_phi,
    n3 = this->n_von_karman_strain_components();
    DenseRealMatrix material_A_mat, material_B_mat, material_D_mat,
    mat1_n1n2, mat2_n2n2, mat3,
    mat4_n3n2, vk_dvdxi_mat, vk_dwdxi_mat, stress, stress_l, local_jac;
    DenseRealVector  vec1_n1, vec2_n1, vec3_n2,
    vec4_n3, vec5_n3, local_f;
    
    mat1_n1n2.resize(n1, n2); mat2_n2n2.resize(n2, n2);
    mat4_n3n2.resize(n3, n2); local_jac.resize(n2, n2);
    vk_dvdxi_mat.resize(n1,n3); vk_dwdxi_mat.resize(n1,n3);
    stress.resize(2,2); stress_l.resize(2, 2);
    local_f.resize(n2); vec1_n1.resize(n1); vec2_n1.resize(n1);
    vec3_n2.resize(n2); vec4_n3.resize(n3); vec5_n3.resize(n3);
    
    
    Bmat_mem.reinit(n1, _system.n_vars(), n_phi); // three stress-strain components
    Bmat_bend.reinit(n1, _system.n_vars(), n_phi);
    Bmat_v_vk.reinit(n3, _system.n_vars(), n_phi); // only dv/dx and dv/dy
    Bmat_w_vk.reinit(n3, _system.n_vars(), n_phi); // only dw/dx and dw/dy
    
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
                                  Bmat_mem, Bmat_bend, Bmat_v_vk, Bmat_w_vk,
                                  stress, stress_l, vk_dvdxi_mat, vk_dwdxi_mat,
                                  material_A_mat,
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
        transform_to_global_system(local_jac, mat2_n2n2);
        jac.add(1., mat2_n2n2);
    }
    
    return request_jacobian;
}





bool
MAST::StructuralElement1D::internal_force_sensitivity (bool request_jacobian,
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
    
    FEMOperatorMatrix Bmat_mem, Bmat_bend, Bmat_v_vk, Bmat_w_vk;
    
    const std::vector<libMesh::Real>& JxW = _fe->get_JxW();
    const std::vector<Point>& xyz = _fe->get_xyz();
    const unsigned int n_phi = (unsigned int)_fe->get_phi().size();
    const unsigned int n1= this->n_direct_strain_components(), n2=6*n_phi,
    n3 = this->n_von_karman_strain_components();
    DenseRealMatrix material_A_mat, material_B_mat, material_D_mat,
    material_trans_shear_mat, mat1_n1n2, mat2_n2n2, mat3,
    mat4_n3n2, vk_dvdxi_mat, vk_dwdxi_mat, stress, stress_l, local_jac;
    DenseRealVector  vec1_n1, vec2_n1, vec3_n2,
    vec4_n3, vec5_n3, local_f;
    
    mat1_n1n2.resize(n1, n2); mat2_n2n2.resize(n2, n2);
    mat4_n3n2.resize(n3, n2); local_jac.resize(n2, n2);
    vk_dvdxi_mat.resize(n1,n3); vk_dwdxi_mat.resize(n1,n3);
    stress.resize(2,2); stress_l.resize(2, 2); local_f.resize(n2);
    vec1_n1.resize(n1); vec2_n1.resize(n1);
    vec3_n2.resize(n2); vec4_n3.resize(n3); vec5_n3.resize(n3);
    
    
    Bmat_mem.reinit(n1, _system.n_vars(), n_phi); // three stress-strain components
    Bmat_bend.reinit(n1, _system.n_vars(), n_phi);
    Bmat_v_vk.reinit(n3, _system.n_vars(), n_phi); // only dv/dx and dv/dy
    Bmat_w_vk.reinit(n3, _system.n_vars(), n_phi); // only dw/dx and dw/dy
    
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
        _internal_force_operation(if_bending, if_vk, n2, qp, JxW,
                                  request_jacobian, if_ignore_ho_jac,
                                  local_f, local_jac,
                                  Bmat_mem, Bmat_bend, Bmat_v_vk, Bmat_w_vk,
                                  stress, stress_l, vk_dvdxi_mat, vk_dwdxi_mat,
                                  material_A_mat,
                                  material_B_mat, material_D_mat, vec1_n1,
                                  vec2_n1, vec3_n2, vec4_n3,
                                  vec5_n3, mat1_n1n2, mat2_n2n2,
                                  mat3, mat4_n3n2);
        
        // this accounts for the sensitivity of the linear stress as a result of
        // static solution. This is needed only for cases that require linearized
        // geometric stiffness matrix, for example in buckling analysis
        if (if_bending && if_vk && if_ignore_ho_jac && request_jacobian) {
            (*mat_stiff_A)(p, _system.time, material_A_mat);
            (*mat_stiff_B)(p, _system.time, material_B_mat);
            
            _linearized_geometric_stiffness_sensitivity_with_static_solution
            (n2, qp, JxW, local_jac,
             Bmat_mem, Bmat_bend, Bmat_v_vk, Bmat_w_vk,
             stress_l, vk_dvdxi_mat, vk_dwdxi_mat,
             material_A_mat, material_B_mat, vec1_n1,
             vec2_n1, mat1_n1n2, mat2_n2n2,
             mat3);
        }

    }
    
    // now calculate the transverse shear contribution if appropriate for the
    // element
    if (if_bending && _bending_operator->include_transverse_shear_energy())
        _bending_operator->calculate_transverse_shear_force(request_jacobian,
                                                            local_f, local_jac,
                                                            this->sensitivity_param);

    // now transform to the global coorodinate system
    transform_to_global_system(local_f, vec3_n2);
    f.add(1., vec3_n2);
    if (request_jacobian) {
        transform_to_global_system(local_jac, mat2_n2n2);
        jac.add(1., mat2_n2n2);
    }
    
    return request_jacobian;
}


void
MAST::StructuralElement1D::_internal_force_operation
(bool if_bending,
 bool if_vk,
 const unsigned int n2,
 const unsigned int qp,
 const std::vector<libMesh::Real>& JxW,
 bool request_jacobian,
 bool if_ignore_ho_jac,
 DenseRealVector& local_f,
 DenseRealMatrix& local_jac,
 FEMOperatorMatrix& Bmat_mem,
 FEMOperatorMatrix& Bmat_bend,
 FEMOperatorMatrix& Bmat_v_vk,
 FEMOperatorMatrix& Bmat_w_vk,
 DenseRealMatrix& stress,
 DenseRealMatrix& stress_l,
 DenseRealMatrix& vk_dvdxi_mat,
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
    stress(0,0)   = vec2_n1(0);
    stress(1,1) = vec1_n1(0); // temporary storage of the membrane strain
    
    // get the bending strain operator
    vec2_n1.zero(); // used to store vk strain, if applicable
    if (if_bending) {
        _bending_operator->initialize_bending_strain_operator(qp, Bmat_bend);

        //  evaluate the bending stress and add that to the stress vector
        // for evaluation in the nonlinear stress term
        Bmat_bend.vector_mult(vec2_n1, local_solution);
        material_B_mat.vector_mult(vec1_n1, vec2_n1);
        stress_l(0,0) += vec1_n1(0);
        stress(0,0)   += vec1_n1(0);

        if (if_vk) {  // get the vonKarman strain operator if needed
            
            this->initialize_von_karman_strain_operator(qp,
                                                        vec2_n1, // epsilon_vk
                                                        vk_dvdxi_mat,
                                                        vk_dwdxi_mat,
                                                        Bmat_v_vk,
                                                        Bmat_w_vk);
            material_A_mat.vector_mult(vec1_n1, vec2_n1);
            stress(0,0) += vec1_n1(0); // total strain that multiplies with the membrane strain
            stress(1,1) += vec2_n1(0); // add the two strains to get the direct strain
        }
    }
    
    // copy the stress to use here.
    vec1_n1.zero();
    vec1_n1(0) = stress(0,0);

    // now the internal force vector
    // this includes the membrane strain operator with all A and B material operators
    Bmat_mem.vector_mult_transpose(vec3_n2, vec1_n1);
    local_f.add(-JxW[qp], vec3_n2);
    
    if (if_bending) {
        if (if_vk) {
            // von Karman strain: direct stress
            vk_dvdxi_mat.vector_mult_transpose(vec4_2, vec1_n1);
            Bmat_v_vk.vector_mult_transpose(vec3_n2, vec4_2);
            local_f.add(-JxW[qp], vec3_n2);
            
            // von Karman strain: direct stress
            vk_dwdxi_mat.vector_mult_transpose(vec4_2, vec1_n1);
            Bmat_w_vk.vector_mult_transpose(vec3_n2, vec4_2);
            local_f.add(-JxW[qp], vec3_n2);
        }
        
        // use the direct strain from the temprary storage
        vec2_n1(0)  = stress(1,1);
        stress(1,1) = 0.;
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
                // membrane - vk: v-displacement
                mat3.resize(vk_dvdxi_mat.m(), n2);
                Bmat_v_vk.left_multiply(mat3, vk_dvdxi_mat);
                mat3.left_multiply(material_A_mat);
                Bmat_mem.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac.add(-JxW[qp], mat2_n2n2);
                
                // membrane - vk: w-displacement
                mat3.resize(vk_dwdxi_mat.m(), n2);
                Bmat_w_vk.left_multiply(mat3, vk_dwdxi_mat);
                mat3.left_multiply(material_A_mat);
                Bmat_mem.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac.add(-JxW[qp], mat2_n2n2);
                
                // vk - membrane: v-displacement
                Bmat_mem.left_multiply(mat1_n1n2, material_A_mat);
                vk_dvdxi_mat.get_transpose(mat3);
                mat3.right_multiply(mat1_n1n2);
                Bmat_v_vk.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac.add(-JxW[qp], mat2_n2n2);
                
                // vk - membrane: w-displacement
                Bmat_mem.left_multiply(mat1_n1n2, material_A_mat);
                vk_dwdxi_mat.get_transpose(mat3);
                mat3.right_multiply(mat1_n1n2);
                Bmat_w_vk.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac.add(-JxW[qp], mat2_n2n2);
                
                // if only the first order term of the Jacobian is needed, for
                // example for linearized buckling analysis, then the linear
                // stress combined with the variation of the von Karman strain
                // is included. Otherwise, all terms are included
                if (if_ignore_ho_jac) {
                    // vk - vk: v-displacement: first order term
                    mat3.resize(2, n2);
                    Bmat_v_vk.left_multiply(mat3, stress_l);
                    Bmat_v_vk.right_multiply_transpose(mat2_n2n2, mat3);
                    local_jac.add(-JxW[qp], mat2_n2n2);

                    // vk - vk: v-displacement: first order term
                    mat3.resize(2, n2);
                    Bmat_w_vk.left_multiply(mat3, stress_l);
                    Bmat_w_vk.right_multiply_transpose(mat2_n2n2, mat3);
                    local_jac.add(-JxW[qp], mat2_n2n2);
                }
                else {
                    // vk - vk: v-displacement
                    mat3.resize(2, n2);
                    Bmat_v_vk.left_multiply(mat3, stress);
                    Bmat_v_vk.right_multiply_transpose(mat2_n2n2, mat3);
                    local_jac.add(-JxW[qp], mat2_n2n2);
                    
                    mat3.resize(vk_dvdxi_mat.m(), n2);
                    Bmat_v_vk.left_multiply(mat3, vk_dvdxi_mat);
                    mat3.left_multiply(material_A_mat);
                    mat3.left_multiply_transpose(vk_dvdxi_mat);
                    Bmat_v_vk.right_multiply_transpose(mat2_n2n2, mat3);
                    local_jac.add(-JxW[qp], mat2_n2n2);

                    // vk - vk: w-displacement
                    mat3.resize(2, n2);
                    Bmat_w_vk.left_multiply(mat3, stress);
                    Bmat_w_vk.right_multiply_transpose(mat2_n2n2, mat3);
                    local_jac.add(-JxW[qp], mat2_n2n2);
                    
                    mat3.resize(vk_dwdxi_mat.m(), n2);
                    Bmat_w_vk.left_multiply(mat3, vk_dwdxi_mat);
                    mat3.left_multiply(material_A_mat);
                    mat3.left_multiply_transpose(vk_dwdxi_mat);
                    Bmat_w_vk.right_multiply_transpose(mat2_n2n2, mat3);
                    local_jac.add(-JxW[qp], mat2_n2n2);
                    
                    // coupling of v, w-displacements
                    mat3.resize(vk_dwdxi_mat.m(), n2);
                    Bmat_w_vk.left_multiply(mat3, vk_dwdxi_mat);
                    mat3.left_multiply(material_A_mat);
                    mat3.left_multiply_transpose(vk_dvdxi_mat);
                    Bmat_v_vk.right_multiply_transpose(mat2_n2n2, mat3);
                    local_jac.add(-JxW[qp], mat2_n2n2);
                    
                    mat3.resize(vk_dvdxi_mat.m(), n2);
                    Bmat_v_vk.left_multiply(mat3, vk_dvdxi_mat);
                    mat3.left_multiply(material_A_mat);
                    mat3.left_multiply_transpose(vk_dwdxi_mat);
                    Bmat_w_vk.right_multiply_transpose(mat2_n2n2, mat3);
                    local_jac.add(-JxW[qp], mat2_n2n2);

                }
                
                // bending - vk: v-displacement
                mat3.resize(vk_dvdxi_mat.m(), n2);
                Bmat_v_vk.left_multiply(mat3, vk_dvdxi_mat);
                mat3.left_multiply_transpose(material_B_mat);
                Bmat_bend.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac.add(-JxW[qp], mat2_n2n2);

                // bending - vk: w-displacement
                mat3.resize(vk_dwdxi_mat.m(), n2);
                Bmat_w_vk.left_multiply(mat3, vk_dwdxi_mat);
                mat3.left_multiply_transpose(material_B_mat);
                Bmat_bend.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac.add(-JxW[qp], mat2_n2n2);

                // vk - bending: v-displacement
                Bmat_bend.left_multiply(mat1_n1n2, material_B_mat);
                vk_dvdxi_mat.get_transpose(mat3);
                mat3.right_multiply(mat1_n1n2);
                Bmat_v_vk.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac.add(-JxW[qp], mat2_n2n2);

                // vk - bending: w-displacement
                Bmat_bend.left_multiply(mat1_n1n2, material_B_mat);
                vk_dwdxi_mat.get_transpose(mat3);
                mat3.right_multiply(mat1_n1n2);
                Bmat_w_vk.right_multiply_transpose(mat2_n2n2, mat3);
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
MAST::StructuralElement1D::_linearized_geometric_stiffness_sensitivity_with_static_solution
(const unsigned int n2,
 const unsigned int qp,
 const std::vector<libMesh::Real>& JxW,
 DenseRealMatrix& local_jac,
 FEMOperatorMatrix& Bmat_mem,
 FEMOperatorMatrix& Bmat_bend,
 FEMOperatorMatrix& Bmat_v_vk,
 FEMOperatorMatrix& Bmat_w_vk,
 DenseRealMatrix& stress_l,
 DenseRealMatrix& vk_dvdxi_mat,
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
    Bmat_mem.vector_mult(vec1_n1, local_solution);
    material_A_mat.vector_mult(vec2_n1, vec1_n1); // linear direct stress
    
    // copy the stress values to a matrix
    stress_l(0,0) = vec2_n1(0); // sigma_xx
    
    // get the von Karman operator matrix
    this->initialize_von_karman_strain_operator(qp,
                                                vec2_n1, // epsilon_vk
                                                vk_dvdxi_mat,
                                                vk_dwdxi_mat,
                                                Bmat_v_vk,
                                                Bmat_w_vk);

    // sensitivity of the vk_dwdxi matrix due to solution sensitivity
    this->initialize_von_karman_strain_operator_sensitivity(qp,
                                                            vk_dvdxi_mat,
                                                            vk_dwdxi_mat);
    
    
    // membrane - vk: v-displacement
    mat3.resize(vk_dvdxi_mat.m(), n2);
    Bmat_v_vk.left_multiply(mat3, vk_dvdxi_mat);
    mat3.left_multiply(material_A_mat);
    Bmat_mem.right_multiply_transpose(mat2_n2n2, mat3);
    local_jac.add(-JxW[qp], mat2_n2n2);
    
    // membrane - vk: w-displacement
    mat3.resize(vk_dwdxi_mat.m(), n2);
    Bmat_w_vk.left_multiply(mat3, vk_dwdxi_mat);
    mat3.left_multiply(material_A_mat);
    Bmat_mem.right_multiply_transpose(mat2_n2n2, mat3);
    local_jac.add(-JxW[qp], mat2_n2n2);
    
    // vk - membrane: v-displacement
    Bmat_mem.left_multiply(mat1_n1n2, material_A_mat);
    vk_dvdxi_mat.get_transpose(mat3);
    mat3.right_multiply(mat1_n1n2);
    Bmat_v_vk.right_multiply_transpose(mat2_n2n2, mat3);
    local_jac.add(-JxW[qp], mat2_n2n2);
    
    // vk - membrane: w-displacement
    Bmat_mem.left_multiply(mat1_n1n2, material_A_mat);
    vk_dwdxi_mat.get_transpose(mat3);
    mat3.right_multiply(mat1_n1n2);
    Bmat_w_vk.right_multiply_transpose(mat2_n2n2, mat3);
    local_jac.add(-JxW[qp], mat2_n2n2);
    
    // vk - vk: v-displacement: first order term
    mat3.resize(2, n2);
    Bmat_v_vk.left_multiply(mat3, stress_l);
    Bmat_v_vk.right_multiply_transpose(mat2_n2n2, mat3);
    local_jac.add(-JxW[qp], mat2_n2n2);
    
    // vk - vk: v-displacement: first order term
    mat3.resize(2, n2);
    Bmat_w_vk.left_multiply(mat3, stress_l);
    Bmat_w_vk.right_multiply_transpose(mat2_n2n2, mat3);
    local_jac.add(-JxW[qp], mat2_n2n2);
    
    // bending - vk: v-displacement
    mat3.resize(vk_dvdxi_mat.m(), n2);
    Bmat_v_vk.left_multiply(mat3, vk_dvdxi_mat);
    mat3.left_multiply_transpose(material_B_mat);
    Bmat_bend.right_multiply_transpose(mat2_n2n2, mat3);
    local_jac.add(-JxW[qp], mat2_n2n2);
    
    // bending - vk: w-displacement
    mat3.resize(vk_dwdxi_mat.m(), n2);
    Bmat_w_vk.left_multiply(mat3, vk_dwdxi_mat);
    mat3.left_multiply_transpose(material_B_mat);
    Bmat_bend.right_multiply_transpose(mat2_n2n2, mat3);
    local_jac.add(-JxW[qp], mat2_n2n2);
    
    // vk - bending: v-displacement
    Bmat_bend.left_multiply(mat1_n1n2, material_B_mat);
    vk_dvdxi_mat.get_transpose(mat3);
    mat3.right_multiply(mat1_n1n2);
    Bmat_v_vk.right_multiply_transpose(mat2_n2n2, mat3);
    local_jac.add(-JxW[qp], mat2_n2n2);
    
    // vk - bending: w-displacement
    Bmat_bend.left_multiply(mat1_n1n2, material_B_mat);
    vk_dwdxi_mat.get_transpose(mat3);
    mat3.right_multiply(mat1_n1n2);
    Bmat_w_vk.right_multiply_transpose(mat2_n2n2, mat3);
    local_jac.add(-JxW[qp], mat2_n2n2);
}



bool
MAST::StructuralElement1D::prestress_force (bool request_jacobian,
                                            DenseRealVector& f,
                                            DenseRealMatrix& jac)
{
    if (!_property.if_prestressed())
        return false;
    
    FEMOperatorMatrix Bmat_mem, Bmat_bend, Bmat_v_vk, Bmat_w_vk;
    
    const std::vector<libMesh::Real>& JxW = _fe->get_JxW();
    const std::vector<Point>& xyz = _fe->get_xyz();
    const unsigned int n_phi = (unsigned int)_fe->get_phi().size();
    const unsigned int n1= this->n_direct_strain_components(), n2=6*n_phi,
    n3 = this->n_von_karman_strain_components();
    DenseRealMatrix mat2_n2n2, mat3, vk_dvdxi_mat, vk_dwdxi_mat,
    local_jac, prestress_mat_A, prestress_mat_B;
    DenseRealVector vec2_n1, vec3_n2, vec4_n3, vec5_n3,
    local_f, prestress_vec_A, prestress_vec_B;
    
    mat2_n2n2.resize(n2, n2); local_jac.resize(n2, n2); vk_dvdxi_mat.resize(n1,n3);
    vk_dwdxi_mat.resize(n1,n3); local_f.resize(n2); vec2_n1.resize(n1);
    vec3_n2.resize(n2); vec4_n3.resize(n3); vec5_n3.resize(n3);
    
    
    Bmat_mem.reinit(n1, _system.n_vars(), n_phi); // three stress-strain components
    Bmat_bend.reinit(n1, _system.n_vars(), n_phi);
    Bmat_v_vk.reinit(n3, _system.n_vars(), n_phi); // only dv/dx and dv/dy
    Bmat_w_vk.reinit(n3, _system.n_vars(), n_phi); // only dw/dx and dw/dy
    
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
                                                            vk_dvdxi_mat,
                                                            vk_dwdxi_mat,
                                                            Bmat_v_vk,
                                                            Bmat_w_vk);
        }
        
        // first handle constant throught the thickness stresses: membrane and vonKarman
        // multiply this with the constant through the thickness strain
        // membrane strain
        Bmat_mem.vector_mult_transpose(vec3_n2, prestress_vec_A);
        local_f.add(-JxW[qp], vec3_n2); // epsilon_mem * sigma_0
        
        if (if_bending) {
            if (if_vk) {
                // von Karman strain: v-displacement
                vk_dvdxi_mat.vector_mult_transpose(vec4_n3, prestress_vec_A);
                Bmat_v_vk.vector_mult_transpose(vec3_n2, vec4_n3);
                local_f.add(-JxW[qp], vec3_n2); // epsilon_vk * sigma_0

                // von Karman strain: w-displacement
                vk_dwdxi_mat.vector_mult_transpose(vec4_n3, prestress_vec_A);
                Bmat_w_vk.vector_mult_transpose(vec3_n2, vec4_n3);
                local_f.add(-JxW[qp], vec3_n2); // epsilon_vk * sigma_0
            }
            
            // now coupling with the bending strain
            Bmat_bend.vector_mult_transpose(vec3_n2, prestress_vec_B);
            local_f.add(-JxW[qp], vec3_n2); // epsilon_bend * sigma_0
        }
        
        if (request_jacobian) {
            if (if_bending && if_vk) {
                // v-displacement
                mat3.resize(2, n2);
                Bmat_v_vk.left_multiply(mat3, prestress_mat_A);
                Bmat_v_vk.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac.add(-JxW[qp], mat2_n2n2);

                // w-displacement
                mat3.resize(2, n2);
                Bmat_w_vk.left_multiply(mat3, prestress_mat_A);
                Bmat_w_vk.right_multiply_transpose(mat2_n2n2, mat3);
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
MAST::StructuralElement1D::prestress_force_sensitivity (bool request_jacobian,
                                                        DenseRealVector& f,
                                                        DenseRealMatrix& jac)
{
    if (!_property.if_prestressed())
        return false;
    
    FEMOperatorMatrix Bmat_mem, Bmat_bend, Bmat_v_vk, Bmat_w_vk;
    
    const std::vector<libMesh::Real>& JxW = _fe->get_JxW();
    const std::vector<Point>& xyz = _fe->get_xyz();
    const unsigned int n_phi = (unsigned int)_fe->get_phi().size();
    const unsigned int n1= this->n_direct_strain_components(), n2=6*n_phi,
    n3 = this->n_von_karman_strain_components();
    DenseRealMatrix mat2_n2n2, mat3, vk_dwdxi_mat, vk_dvdxi_mat,
    local_jac, prestress_mat_A, prestress_mat_B;
    DenseRealVector vec2_n1, vec3_n2, vec4_n3, vec5_n3,
    local_f, prestress_vec_A, prestress_vec_B;
    
    mat2_n2n2.resize(n2, n2); local_jac.resize(n2, n2); vk_dvdxi_mat.resize(n1,n3);
    vk_dwdxi_mat.resize(n1,n3); local_f.resize(n2); vec2_n1.resize(n1);
    vec3_n2.resize(n2); vec4_n3.resize(n3); vec5_n3.resize(n3);
    
    
    Bmat_mem.reinit(n1, _system.n_vars(), n_phi); // three stress-strain components
    Bmat_bend.reinit(n1, _system.n_vars(), n_phi);
    Bmat_v_vk.reinit(n3, _system.n_vars(), n_phi); // only dv/dx and dv/dy
    Bmat_w_vk.reinit(n3, _system.n_vars(), n_phi); // only dw/dx and dw/dy
    
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
                                                            vk_dvdxi_mat,
                                                            vk_dwdxi_mat,
                                                            Bmat_v_vk,
                                                            Bmat_w_vk);
        }
        
        // first handle constant throught the thickness stresses: membrane and vonKarman
        // multiply this with the constant through the thickness strain
        // membrane strain
        Bmat_mem.vector_mult_transpose(vec3_n2, prestress_vec_A);
        local_f.add(-JxW[qp], vec3_n2); // epsilon_mem * sigma_0
        
        if (if_bending) {
            if (if_vk) {
                // von Karman strain: v-displacement
                vk_dvdxi_mat.vector_mult_transpose(vec4_n3, prestress_vec_A);
                Bmat_v_vk.vector_mult_transpose(vec3_n2, vec4_n3);
                local_f.add(-JxW[qp], vec3_n2); // epsilon_vk * sigma_0

                // von Karman strain: w-displacement
                vk_dwdxi_mat.vector_mult_transpose(vec4_n3, prestress_vec_A);
                Bmat_w_vk.vector_mult_transpose(vec3_n2, vec4_n3);
                local_f.add(-JxW[qp], vec3_n2); // epsilon_vk * sigma_0
            }
            
            // now coupling with the bending strain
            Bmat_bend.vector_mult_transpose(vec3_n2, prestress_vec_B);
            local_f.add(-JxW[qp], vec3_n2); // epsilon_bend * sigma_0
        }
        
        if (request_jacobian) {
            if (if_bending && if_vk) {
                // v-displacement
                mat3.resize(2, n2);
                Bmat_v_vk.left_multiply(mat3, prestress_mat_A);
                Bmat_v_vk.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac.add(-JxW[qp], mat2_n2n2);

                // w-displacement
                mat3.resize(2, n2);
                Bmat_w_vk.left_multiply(mat3, prestress_mat_A);
                Bmat_w_vk.right_multiply_transpose(mat2_n2n2, mat3);
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
MAST::StructuralElement1D::thermal_force (bool request_jacobian,
                                          DenseRealVector& f,
                                          DenseRealMatrix& jac,
                                          MAST::BoundaryCondition& p)
{
    FEMOperatorMatrix Bmat_mem, Bmat_bend, Bmat_v_vk, Bmat_w_vk;
    
    const std::vector<libMesh::Real>& JxW = _fe->get_JxW();
    const std::vector<Point>& xyz = _fe->get_xyz();
    const unsigned int n_phi = (unsigned int)_fe->get_phi().size();
    const unsigned int n1= this->n_direct_strain_components(), n2=6*n_phi,
    n3 = this->n_von_karman_strain_components();
    DenseRealMatrix material_exp_A_mat, material_exp_B_mat,
    mat1_n1n2, mat2_n2n2, mat3,
    mat4_n3n2, vk_dvdxi_mat, vk_dwdxi_mat, stress, local_jac;
    DenseRealVector  vec1_n1, vec2_n1, vec3_n2,
    vec4_2, vec5_n3, local_f, delta_t;
    
    mat1_n1n2.resize(n1, n2); mat2_n2n2.resize(n2, n2);
    mat4_n3n2.resize(n3, n2); local_jac.resize(n2, n2);
    vk_dwdxi_mat.resize(n1,n3); stress.resize(2,2); local_f.resize(n2);
    vec1_n1.resize(n1); vec2_n1.resize(n1);
    vec3_n2.resize(n2); vec4_2.resize(2); vec5_n3.resize(n3);
    delta_t.resize(1); vk_dvdxi_mat.resize(2, 2); vk_dwdxi_mat.resize(2, 2);
    
    
    Bmat_mem.reinit(n1, _system.n_vars(), n_phi); // three stress-strain components
    Bmat_bend.reinit(n1, _system.n_vars(), n_phi);
    Bmat_v_vk.reinit(n3, _system.n_vars(), n_phi); // only dv/dx and dv/dy
    Bmat_w_vk.reinit(n3, _system.n_vars(), n_phi); // only dw/dx and dw/dy
    
    bool if_vk = (_property.strain_type() == MAST::VON_KARMAN_STRAIN),
    if_bending = (_property.bending_model(_elem, _fe->get_fe_type()) != MAST::NO_BENDING);
    
    std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > > expansion_A
    (_property.get_property(MAST::SECTION_INTEGRATED_MATERIAL_THERMAL_EXPANSION_A_MATRIX,
                            *this).release()),
    expansion_B
    (_property.get_property(MAST::SECTION_INTEGRATED_MATERIAL_THERMAL_EXPANSION_B_MATRIX,
                            *this).release());
    
    std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > > mat
    (_property.get_property(MAST::SECTION_INTEGRATED_MATERIAL_THERMAL_EXPANSION_A_MATRIX,
                            *this).release());
    
    // temperature function
    MAST::FieldFunction<libMesh::Real>& temp_func =
    dynamic_cast<MAST::FieldFunction<libMesh::Real>&>(p.function());
    MAST::FieldFunction<libMesh::Real>& ref_temp_func =
    dynamic_cast<MAST::FieldFunction<libMesh::Real>&>
    (dynamic_cast<MAST::Temperature&>(p).reference_temperature_function());
    
    libMesh::Real t, t0;
    libMesh::Point pt;
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        this->global_coordinates(xyz[qp], pt);
        
        // get the material property
        (*expansion_A)(pt, _system.time, material_exp_A_mat);
        (*expansion_B)(pt, _system.time, material_exp_B_mat);
        
        // get the temperature function
        temp_func(xyz[qp], _system.time, t);
        ref_temp_func(xyz[qp], _system.time, t0);
        delta_t(0) = t-t0;

        material_exp_A_mat.vector_mult(vec1_n1, delta_t); // [C]{alpha (T - T0)} (with membrane strain)
        stress(0,0) = vec1_n1(0); // sigma_xx
        material_exp_B_mat.vector_mult(vec2_n1, delta_t); // [C]{alpha (T - T0)} (with bending strain)
        
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
                                                            vk_dvdxi_mat,
                                                            vk_dwdxi_mat,
                                                            Bmat_v_vk,
                                                            Bmat_w_vk);
                // von Karman strain: v-displacement
                vk_dvdxi_mat.vector_mult_transpose(vec4_2, vec1_n1);
                Bmat_v_vk.vector_mult_transpose(vec3_n2, vec4_2);
                local_f.add(JxW[qp], vec3_n2);

                // von Karman strain: w-displacement
                vk_dwdxi_mat.vector_mult_transpose(vec4_2, vec1_n1);
                Bmat_w_vk.vector_mult_transpose(vec3_n2, vec4_2);
                local_f.add(JxW[qp], vec3_n2);
            }
            
            if (request_jacobian && if_vk) { // Jacobian only for vk strain
                
                // vk - vk: v-displacement
                mat3.resize(2, n2);
                Bmat_v_vk.left_multiply(mat3, stress);
                Bmat_v_vk.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac.add(JxW[qp], mat2_n2n2);
                
                // vk - vk: w-displacement
                mat3.resize(2, n2);
                Bmat_w_vk.left_multiply(mat3, stress);
                Bmat_w_vk.right_multiply_transpose(mat2_n2n2, mat3);
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
MAST::StructuralElement1D::thermal_force_sensitivity (bool request_jacobian,
                                                      DenseRealVector& f,
                                                      DenseRealMatrix& jac,
                                                      MAST::BoundaryCondition& p)
{
    FEMOperatorMatrix Bmat_mem, Bmat_bend, Bmat_v_vk, Bmat_w_vk;
    
    const std::vector<libMesh::Real>& JxW = _fe->get_JxW();
    const std::vector<Point>& xyz = _fe->get_xyz();
    const unsigned int n_phi = (unsigned int)_fe->get_phi().size();
    const unsigned int n1= this->n_direct_strain_components(), n2=6*n_phi,
    n3 = this->n_von_karman_strain_components();
    DenseRealMatrix material_exp_A_mat, material_exp_B_mat,
    material_exp_A_mat_sens, material_exp_B_mat_sens,
    mat1_n1n2, mat2_n2n2, mat3,
    mat4_n3n2, vk_dvdxi_mat, vk_dwdxi_mat, stress, local_jac;
    DenseRealVector  vec1_n1, vec2_n1, vec3_n2,
    vec4_2, vec5_n1, local_f, delta_t, delta_t_sens;
    
    mat1_n1n2.resize(n1, n2); mat2_n2n2.resize(n2, n2);
    mat4_n3n2.resize(n3, n2); local_jac.resize(n2, n2);
    vk_dwdxi_mat.resize(n1,n3); stress.resize(2,2); local_f.resize(n2);
    vec1_n1.resize(n1); vec2_n1.resize(n1);
    vec3_n2.resize(n2); vec4_2.resize(2); vec5_n1.resize(n1);
    delta_t.resize(1); delta_t_sens.resize(1); vk_dvdxi_mat.resize(2, 2);
    vk_dwdxi_mat.resize(2, 2);
    
    
    Bmat_mem.reinit(n1, _system.n_vars(), n_phi); // three stress-strain components
    Bmat_bend.reinit(n1, _system.n_vars(), n_phi);
    Bmat_v_vk.reinit(n3, _system.n_vars(), n_phi); // only dv/dx and dv/dy
    Bmat_w_vk.reinit(n3, _system.n_vars(), n_phi); // only dw/dx and dw/dy
    
    bool if_vk = (_property.strain_type() == MAST::VON_KARMAN_STRAIN),
    if_bending = (_property.bending_model(_elem, _fe->get_fe_type()) != MAST::NO_BENDING);
    
    std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > > expansion_A
    (_property.get_property(MAST::SECTION_INTEGRATED_MATERIAL_THERMAL_EXPANSION_A_MATRIX,
                            *this).release()),
    expansion_B
    (_property.get_property(MAST::SECTION_INTEGRATED_MATERIAL_THERMAL_EXPANSION_B_MATRIX,
                            *this).release());
    
    std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > > mat
    (_property.get_property(MAST::SECTION_INTEGRATED_MATERIAL_THERMAL_EXPANSION_A_MATRIX,
                            *this).release());
    
    // temperature function
    MAST::FieldFunction<libMesh::Real>& temp_func =
    dynamic_cast<MAST::FieldFunction<libMesh::Real>&>(p.function());
    MAST::FieldFunction<libMesh::Real>& ref_temp_func =
    dynamic_cast<MAST::FieldFunction<libMesh::Real>&>
    (dynamic_cast<MAST::Temperature&>(p).reference_temperature_function());
    
    libMesh::Real t, t0, t_sens;
    libMesh::Point pt;
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        this->global_coordinates(xyz[qp], pt);
        
        // get the material property
        (*expansion_A)(pt, _system.time, material_exp_A_mat);
        (*expansion_B)(pt, _system.time, material_exp_B_mat);
        expansion_A->total(*this->sensitivity_param,
                           pt, _system.time, material_exp_A_mat_sens);
        expansion_B->total(*this->sensitivity_param,
                           pt, _system.time, material_exp_B_mat_sens);
        
        // get the temperature function
        temp_func(xyz[qp], _system.time, t);
        temp_func.total(*this->sensitivity_param,
                        xyz[qp], _system.time, t_sens);
        ref_temp_func(xyz[qp], _system.time, t0);
        delta_t(0) = t-t0;
        delta_t_sens(0) = t_sens;
        
        // now prepare the membrane force sensitivity
        material_exp_A_mat.vector_mult(vec1_n1, delta_t_sens); // [C]{alpha (dT/dp)} (with membrane strain)
        material_exp_A_mat_sens.vector_mult(vec2_n1, delta_t); // d([C].{alpha})/dp (T - T0)} (with membrane
        vec1_n1.add(1., vec2_n1);  // sensitivity of the thermal membrane force
        stress(0,0) = vec1_n1(0); // sigma_xx
        
        // now prepare the membrane-bending coupling force sensitivity
        material_exp_B_mat.vector_mult(vec2_n1, delta_t_sens); // [C]{alpha dT/dp} (with bending strain)
        material_exp_B_mat_sens.vector_mult(vec5_n1, delta_t); // d([C].{alpha})/dp (T - T0)} (with bending strain)
        vec2_n1.add(1., vec5_n1);  // sensitivity of the thermal membrane force
        
        
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
                                                            vk_dvdxi_mat,
                                                            vk_dwdxi_mat,
                                                            Bmat_v_vk,
                                                            Bmat_w_vk);
                // von Karman strain: v-displacement
                vk_dvdxi_mat.vector_mult_transpose(vec4_2, vec1_n1);
                Bmat_v_vk.vector_mult_transpose(vec3_n2, vec4_2);
                local_f.add(JxW[qp], vec3_n2);
                
                // von Karman strain: w-displacement
                vk_dwdxi_mat.vector_mult_transpose(vec4_2, vec1_n1);
                Bmat_w_vk.vector_mult_transpose(vec3_n2, vec4_2);
                local_f.add(JxW[qp], vec3_n2);
            }
            
            if (request_jacobian && if_vk) { // Jacobian only for vk strain
                // vk - vk: v-displacement
                mat3.resize(2, n2);
                Bmat_v_vk.left_multiply(mat3, stress);
                Bmat_v_vk.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac.add(JxW[qp], mat2_n2n2);
                
                // vk - vk: w-displacement
                mat3.resize(2, n2);
                Bmat_w_vk.left_multiply(mat3, stress);
                Bmat_w_vk.right_multiply_transpose(mat2_n2n2, mat3);
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


