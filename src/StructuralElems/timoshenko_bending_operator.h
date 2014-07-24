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

#ifndef __MAST_timoshenko_bending_operator_h__
#define __MAST_timoshenko_bending_operator_h__

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
    class TimoshenkoBendingOperator: public MAST::BendingOperator1D {
    public:
        TimoshenkoBendingOperator(StructuralElementBase& elem):
        MAST::BendingOperator1D(elem),
        _fe(elem.fe()),
        _shear_quadrature_reduction(2)
        { }
        
        virtual ~TimoshenkoBendingOperator() { }
        
        /*!
         *   returns true if this bending operator supports a transverse shear component
         */
        virtual bool include_transverse_shear_energy() const {
            return true;
        }
        
        /*!
         *   initialze the bending strain operator for Timoshenko beam element, withouth
         *   the y,z-location. This is useful for use with element stiffness matrix
         *   integration where the D matrix is calculated by section integration by
         *   the ElementPropertyCard1D.
         */
        virtual void initialize_bending_strain_operator (const unsigned int qp,
                                                         FEMOperatorMatrix& Bmat);
        
        /*!
         *    initializes the bending strain operator for the specified quadrature
         * point and y,z-location.
         */
        virtual void initialize_bending_strain_operator_for_yz(const unsigned int qp,
                                                               const Real y,
                                                               const Real z,
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
MAST::TimoshenkoBendingOperator::initialize_bending_strain_operator(const unsigned int qp,
                                                                    FEMOperatorMatrix& Bmat_bend) {
    this->initialize_bending_strain_operator_for_yz(qp, 1., 1., Bmat_bend);
}




inline void
MAST::TimoshenkoBendingOperator::initialize_bending_strain_operator_for_yz(const unsigned int qp,
                                                                           const Real y,
                                                                           const Real z,
                                                                           FEMOperatorMatrix& Bmat_bend) {
    
    const std::vector<std::vector<libMesh::RealVectorValue> >& dphi = _fe.get_dphi();
    const std::vector<std::vector<Real> >& phi = _fe.get_phi();
    
    const unsigned int n_phi = (unsigned int)phi.size();
    
    DenseRealVector phi_vec; phi_vec.resize(n_phi);
    
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
        phi_vec(i_nd) = dphi[i_nd][qp](0);  // dphi/dx
    phi_vec.scale(-y);
    Bmat_bend.set_shape_function(0, 5, phi_vec); // v-bending: thetaz

    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
        phi_vec(i_nd) = dphi[i_nd][qp](0);  // dphi/dx
    phi_vec.scale(z);
    Bmat_bend.set_shape_function(1, 4, phi_vec); // w-bending : thetay
}



void
MAST::TimoshenkoBendingOperator::calculate_transverse_shear_force
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
    libMesh::FEType fe_type = _fe.get_fe_type();
    
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
    
    const std::vector<std::vector<libMesh::RealVectorValue> >& dphi = fe->get_dphi();
    const std::vector<std::vector<Real> >& phi = fe->get_phi();
    const std::vector<Real>& JxW = fe->get_JxW();
    const std::vector<libMesh::Point>& xyz = fe->get_xyz();
    
    const unsigned int n_phi = (unsigned int)phi.size(), n2 = 6*n_phi;
    DenseRealVector phi_vec; phi_vec.resize(n_phi);
    
    DenseRealVector vec3_n2, vec4_2, vec5_2;
    DenseRealMatrix material_trans_shear_mat, mat2_n2n2, mat4_2n2;
    
    vec3_n2.resize(n2); vec4_2.resize(2); vec5_2.resize(2);
    mat2_n2n2.resize(n2, n2); mat4_2n2.resize(2, n2);
    
    
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
                             xyz[qp], _structural_elem.system().time,
                             material_trans_shear_mat);
        
        // initialize the strain operator
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = dphi[i_nd][qp](0);  // dphi/dx
        Bmat_trans.set_shape_function(0, 1, phi_vec); // gamma-xy:  v
        Bmat_trans.set_shape_function(1, 2, phi_vec); // gamma-xz : w
        
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];  // phi
        Bmat_trans.set_shape_function(1, 4, phi_vec); // gamma-xy:  thetay
        phi_vec.scale(-1.0);
        Bmat_trans.set_shape_function(0, 5, phi_vec); // gamma-xz : thetaz
        
        
        // now add the transverse shear component
        Bmat_trans.vector_mult(vec4_2, _structural_elem.local_solution);
        material_trans_shear_mat.vector_mult(vec5_2, vec4_2);
        Bmat_trans.vector_mult_transpose(vec3_n2, vec5_2);
        local_f.add(-JxW[qp], vec3_n2);
        
        if (request_jacobian) {
            // now add the transverse shear component
            Bmat_trans.left_multiply(mat4_2n2, material_trans_shear_mat);
            Bmat_trans.right_multiply_transpose(mat2_n2n2, mat4_2n2);
            local_jac.add(-JxW[qp], mat2_n2n2);
        }
    }
}


#endif // __MAST_timoshenko_bending_operator_h__
