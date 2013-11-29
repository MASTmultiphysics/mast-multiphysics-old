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
#include "ThermalElems/temperature_function.h"
#include "BoundaryConditions/boundary_condition.h"


void
MAST::Local1DElem::_create_local_elem() {
    
    libmesh_assert(_elem.dim() == 1);
    
    // first node is the origin of the new cs
    // calculate the coordinate system for the plane of the element
    Point v1, v2, v3, p;
    v1 = *_elem.get_node(1); v1 -= *_elem.get_node(0); v1 /= v1.size(); // local x
    v2(1) = 1.;
    v3 = v1.cross(v2); v3 /= v3.size();      // local z
    v2 = v3.cross(v1); v2 /= v2.size();      // local y
    
    _T_mat.resize(3,3);
    
    // if element is in xy-plane, no need to create a new element
    if ((v2(1) == 1.) && v3(2) == 1.) {
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



MAST::StructuralElement1D::StructuralElement1D(System& sys,
                                               const Elem& elem,
                                               const MAST::ElementPropertyCardBase& p):
MAST::BendingStructuralElem(sys, elem, p)
{ }






void
MAST::StructuralElement1D::initialize_direct_strain_operator(const unsigned int qp,
                                                             FEMOperatorMatrix& Bmat) {
    
    const std::vector<std::vector<RealVectorValue> >& dphi = _fe->get_dphi();
    
    unsigned int n_phi = (unsigned int)dphi.size();
    DenseVector<Real> phi; phi.resize(n_phi);
    
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
                                                                 DenseVector<Real>& vk_strain,
                                                                 DenseMatrix<Real>& vk_dwdxi_mat,
                                                                 FEMOperatorMatrix& Bmat_vk) {

    const std::vector<std::vector<RealVectorValue> >& dphi = _fe->get_dphi();
    const unsigned int n_phi = (unsigned int)dphi.size();

    libmesh_assert_equal_to(vk_strain.size(), 2);
    libmesh_assert_equal_to(vk_dwdxi_mat.m(), 2);
    libmesh_assert_equal_to(vk_dwdxi_mat.n(), 2);
    libmesh_assert_equal_to(Bmat_vk.m(), 2);
    libmesh_assert_equal_to(Bmat_vk.n(), 6*n_phi);
    
    Real dv=0., dw=0.;
    vk_dwdxi_mat.zero();
    
    DenseVector<Real> phi_vec; phi_vec.resize(n_phi);
    
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ ) {
        phi_vec(i_nd) = dphi[i_nd][qp](0);                // dphi/dx
        dv += phi_vec(i_nd)*local_solution(n_phi+i_nd);   // dv/dx
        dw += phi_vec(i_nd)*local_solution(2*n_phi+i_nd); // dw/dx
    }

    Bmat_vk.set_shape_function(0, 1, phi_vec); // dv/dx
    Bmat_vk.set_shape_function(0, 2, phi_vec); // dw/dx
    vk_dwdxi_mat(0, 0) = dv;                   // epsilon-xx : dv/dx
    vk_dwdxi_mat(1, 1) = dw;                   // epsilon-xx : dw/dx
    vk_strain(0) = 0.5*(dv*dv+dw*dw);          // 1/2 * [(dv/dx)^2 + (dw/dx)^2]
}

