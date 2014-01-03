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
#include "ThermalElems/temperature_function.h"
#include "BoundaryConditions/boundary_condition.h"


void
MAST::Local2DElem::_create_local_elem() {
    
    libmesh_assert(_elem.dim() == 2);
    
    // first node is the origin of the new cs
    // calculate the coordinate system for the plane of the element
    Point v1, v2, v3, p;
    v1 = *_elem.get_node(1); v1 -= *_elem.get_node(0); v1 /= v1.size(); // local x
    v2 = *_elem.get_node(2); v2 -= *_elem.get_node(0); v2 /= v2.size();
    v3 = v1.cross(v2); v3 /= v3.size();      // local z
    v2 = v3.cross(v1); v2 /= v2.size();      // local y
    
    _T_mat.resize(3,3);

    // if element is in xy-plane, no need to create a new element
    if (v3(2) == 1.) {
        for (unsigned int i=0; i<3; i++)
            _T_mat(i,i) = 1.;
        return;
    }
    
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



MAST::StructuralElement2D::StructuralElement2D(System& sys,
                                               const Elem& elem,
                                               const MAST::ElementPropertyCardBase& p):
MAST::BendingStructuralElem(sys, elem, p)
{ }






void
MAST::StructuralElement2D::initialize_direct_strain_operator(const unsigned int qp,
                                                             FEMOperatorMatrix& Bmat) {
    
    const std::vector<std::vector<RealVectorValue> >& dphi = _fe->get_dphi();
    
    unsigned int n_phi = (unsigned int)dphi.size();
    DenseVector<Real> phi; phi.resize(n_phi);

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
                                                                 DenseVector<Real>& vk_strain,
                                                                 DenseMatrix<Real>& vk_dwdxi_mat,
                                                                 FEMOperatorMatrix& Bmat_vk) {
    
    const std::vector<std::vector<RealVectorValue> >& dphi = _fe->get_dphi();
    const unsigned int n_phi = (unsigned int)dphi.size();
    
    libmesh_assert_equal_to(vk_strain.size(), 3);
    libmesh_assert_equal_to(vk_dwdxi_mat.m(), 3);
    libmesh_assert_equal_to(vk_dwdxi_mat.n(), 2);
    libmesh_assert_equal_to(Bmat_vk.m(), 2);
    libmesh_assert_equal_to(Bmat_vk.n(), 6*n_phi);

    Real dw=0.;
    vk_dwdxi_mat.zero();
    
    DenseVector<Real> phi_vec; phi_vec.resize(n_phi);
    
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



Real
MAST::StructuralElement2D::max_von_mises_stress() {
    Real val = 0., max_val = 0.;
    DenseMatrix<Real> s; s.resize(3,3);
    
    FEMOperatorMatrix Bmat_mem, Bmat_bend, Bmat_vk;
    
    const std::vector<Real>& JxW = _fe->get_JxW();
    const std::vector<Point>& xyz = _fe->get_xyz();
    const unsigned int n_phi = (unsigned int)_fe->get_phi().size();
    const unsigned int n1= this->n_direct_strain_components(), n2=6*n_phi,
    n3 = this->n_von_karman_strain_components();
    DenseMatrix<Real> material_mat, tmp_mat1_n1n2, tmp_mat2_n2n2, tmp_mat3,
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
    
    bool if_plane_stress = dynamic_cast<const MAST::ElementPropertyCard2D&>(_property).plane_stress();
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        // if temperature is specified, the initialize it to the current location
        if (_temperature)
            _temperature->initialize(xyz[qp]);
        
        _property.get_material().calculate_2d_matrix(MAST::MATERIAL_STIFFNESS_MATRIX,
                                                     material_mat,
                                                     if_plane_stress);
        
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
        material_mat.vector_mult(tmp_vec1_n1, tmp_vec2_n1); // stress

        // copy the stress values to a matrix
        s(0,0) = tmp_vec1_n1(0); // sigma_xx
        if (_elem.dim() == 2) { // this is not needed for 1D element
            s(0,1) = tmp_vec1_n1(2); // sigma_xy
            s(1,0) = tmp_vec1_n1(2); // sigma_yx
            s(1,1) = tmp_vec1_n1(1); // sigma_yy
        }

        // bending stress is calculated at two extreme positions
        if (if_bending) {
            
        }
        
        // von Mises stress from the stress tensor
        val =
        pow(s(0,0) - s(1,1), 2) +
        pow(s(1,1) - s(2,2), 2) +
        pow(s(2,2) - s(0,0), 2) +
        6.* (pow(s(0,1), 2) +
             pow(s(1,2), 2) +
             pow(s(2,0), 2));
        
        val = sqrt(.5 * val);
        
        // store the maximum value
        max_val = std::max(val, max_val);
    }
    
    return max_val;
}

