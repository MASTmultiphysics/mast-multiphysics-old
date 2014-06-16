//
//  mindlin_bending_operator.h
//  MAST
//
//  Created by Manav Bhatia on 11/27/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_mindlin_bending_operator_h__
#define __MAST_mindlin_bending_operator_h__

// libMesh includes
#include "libmesh/fe.h"
#include "libmesh/quadrature.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_matrix.h"

// MAST includes
#include "StructuralElems/bending_operator.h"
#include "PropertyCards/element_property_card_base.h"
#include "Numerics/fem_operator_matrix.h"


namespace MAST {
    class MindlinBendingOperator: public MAST::BendingOperator2D {
    public:
        MindlinBendingOperator(StructuralElementBase& elem):
        MAST::BendingOperator2D(elem),
        _fe(elem.fe()),
        _shear_quadrature_reduction(1)
        { }
        
        virtual ~MindlinBendingOperator() { }
        
        /*!
         *   returns true if this bending operator supports a transverse shear component
         */
        virtual bool include_transverse_shear_energy() const {
            return true;
        }
        
        /*!
         *   initialze the bending strain operator for Mindlin element, withouth 
         *   the z-location. This is useful for use with element stiffness matrix
         *   integration where the D matrix is calculated by section integration by
         *   the ElementPropertyCard2D.
         */
        virtual void initialize_bending_strain_operator (const unsigned int qp,
                                                         FEMOperatorMatrix& Bmat);
        
        /*!
         *    initializes the bending strain operator for the specified quadrature
         * point and z-location.
         */
        void initialize_bending_strain_operator_for_z(const unsigned int qp,
                                                      const libMesh::Real z,
                                                      FEMOperatorMatrix& Bmat_bend);
        /*!
         *   calculate the transverse shear component for the element
         */
        virtual void calculate_transverse_shear_force(bool request_jacobian,
                                                      DenseRealVector& local_f,
                                                      DenseRealMatrix& local_jac,
                                                      const MAST::FieldFunctionBase* sens_params );
        
    protected:
        
        /*!
         *   reference to finite elmement object for the element
         */
        libMesh::FEBase& _fe;

        /*!
         *   reduction in quadrature for shear energy
         */
        unsigned int _shear_quadrature_reduction;
    };
}




inline void
MAST::MindlinBendingOperator::initialize_bending_strain_operator(const unsigned int qp,
                                                                 FEMOperatorMatrix& Bmat_bend) {
    this->initialize_bending_strain_operator_for_z(qp, 1., Bmat_bend);
}



inline void
MAST::MindlinBendingOperator::initialize_bending_strain_operator_for_z(const unsigned int qp,
                                                                       const libMesh::Real z,
                                                                       FEMOperatorMatrix& Bmat_bend) {
    
    const std::vector<std::vector<RealVectorValue> >& dphi = _fe.get_dphi();
    const std::vector<std::vector<libMesh::Real> >& phi = _fe.get_phi();
    
    const unsigned int n_phi = (unsigned int)phi.size();
    
    DenseRealVector phi_vec; phi_vec.resize(n_phi);
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
        phi_vec(i_nd) = dphi[i_nd][qp](0);  // dphi/dx
    
    phi_vec.scale(z);
    Bmat_bend.set_shape_function(0, 4, phi_vec); // epsilon-x: thetay
    phi_vec.scale(-1.0);
    Bmat_bend.set_shape_function(2, 3, phi_vec); // gamma-xy : thetax
    
    
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
        phi_vec(i_nd) = dphi[i_nd][qp](1);  // dphi/dy
    
    phi_vec.scale(z);
    Bmat_bend.set_shape_function(2, 4, phi_vec); // gamma-xy : thetay
                                                 //Bmat_trans.set_shape_function(1, 2, phi_vec); // gamma-yz : w
    phi_vec.scale(-1.0);
    Bmat_bend.set_shape_function(1, 3, phi_vec); // epsilon-y: thetax
    
    
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
        phi_vec(i_nd) = phi[i_nd][qp];  // phi
    
    phi_vec.scale(-1.0);
}




void
MAST::MindlinBendingOperator::calculate_transverse_shear_force
(bool request_jacobian,
 DenseRealVector& local_f,
 DenseRealMatrix& local_jac,
 const MAST::FieldFunctionBase* sens_param)
{
    const MAST::ElementPropertyCardBase& property = _structural_elem.elem_property();
    
    // make an fe and quadrature object for the requested order for integrating
    // transverse shear
    
    std::auto_ptr<libMesh::FEBase> fe;
    std::auto_ptr<libMesh::QBase> qrule;
    FEType fe_type = _fe.get_fe_type();
    
    fe.reset(libMesh::FEBase::build(_elem.dim(), fe_type).release());
    qrule.reset(fe_type.default_quadrature_rule
                (_elem.dim(),
                 property.extra_quadrature_order(_elem, fe->get_fe_type())
                 - _shear_quadrature_reduction).release());
    fe->attach_quadrature_rule(qrule.get());
    fe->get_phi();
    fe->get_JxW();
    fe->get_dphi();
    
    fe->reinit(&_elem);
    
    const std::vector<std::vector<RealVectorValue> >& dphi = fe->get_dphi();
    const std::vector<std::vector<libMesh::Real> >& phi = fe->get_phi();
    const std::vector<libMesh::Real>& JxW = fe->get_JxW();
    const std::vector<Point>& xyz = fe->get_xyz();
    
    const unsigned int n_phi = (unsigned int)phi.size(), n2 = 6*n_phi;
    DenseRealVector phi_vec; phi_vec.resize(n_phi);
    
    DenseRealVector tmp_vec3_n2, tmp_vec4_2, tmp_vec5_2;
    DenseRealMatrix material_trans_shear_mat, tmp_mat2_n2n2, tmp_mat4_2n2;
    
    tmp_vec3_n2.resize(n2); tmp_vec4_2.resize(2); tmp_vec5_2.resize(2);
    tmp_mat2_n2n2.resize(n2, n2); tmp_mat4_2n2.resize(2, n2);
    
    
    FEMOperatorMatrix Bmat_trans;
    Bmat_trans.reinit(2, 6, n_phi); // only two shear stresses
    
    std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > > mat_stiff
    (property.get_property(MAST::SECTION_INTEGRATED_MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX,
                           _structural_elem).release());
    
    libMesh::Point p;
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        _structural_elem.global_coordinates(xyz[qp], p);
        
        if (!sens_param)
            (*mat_stiff)(p, _structural_elem.system().time,
                         material_trans_shear_mat);
        else
            mat_stiff->total(*sens_param,
                             p, _structural_elem.system().time,
                             material_trans_shear_mat);
        
        
        
        // initialize the strain operator
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = dphi[i_nd][qp](0);  // dphi/dx
        
        Bmat_trans.set_shape_function(0, 2, phi_vec); // gamma-xz:  w
        
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = dphi[i_nd][qp](1);  // dphi/dy
        
        Bmat_trans.set_shape_function(1, 2, phi_vec); // gamma-yz : w
        
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];  // phi
        
        Bmat_trans.set_shape_function(0, 4, phi_vec); // gamma-xz:  thetay
        phi_vec.scale(-1.0);
        Bmat_trans.set_shape_function(1, 3, phi_vec); // gamma-yz : thetax
        
        
        // now add the transverse shear component
        Bmat_trans.vector_mult(tmp_vec4_2, _structural_elem.local_solution);
        material_trans_shear_mat.vector_mult(tmp_vec5_2, tmp_vec4_2);
        Bmat_trans.vector_mult_transpose(tmp_vec3_n2, tmp_vec5_2);
        local_f.add(-JxW[qp], tmp_vec3_n2);
        
        if (request_jacobian) {
            // now add the transverse shear component
            Bmat_trans.left_multiply(tmp_mat4_2n2, material_trans_shear_mat);
            Bmat_trans.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat4_2n2);
            local_jac.add(-JxW[qp], tmp_mat2_n2n2);
        }
    }
}


#endif  // __MAST_mindlin_bending_operator_h__
