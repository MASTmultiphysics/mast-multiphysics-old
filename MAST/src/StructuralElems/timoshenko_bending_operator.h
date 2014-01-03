//
//  timoshenko_bending_operator.h
//  MAST
//
//  Created by Manav Bhatia on 11/27/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

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
                                                      DenseVector<Real>& local_f,
                                                      DenseMatrix<Real>& local_jac,
                                                      const MAST::SensitivityParameters* sens_params );
        
    protected:
        
        /*!
         *   reference to finite elmement object for the element
         */
        FEBase& _fe;
        
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
    
    const std::vector<std::vector<RealVectorValue> >& dphi = _fe.get_dphi();
    const std::vector<std::vector<Real> >& phi = _fe.get_phi();
    
    const unsigned int n_phi = (unsigned int)phi.size();
    
    DenseVector<Real> phi_vec; phi_vec.resize(n_phi);
    
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
 DenseVector<Real>& local_f,
 DenseMatrix<Real>& local_jac,
 const MAST::SensitivityParameters* sens_params)
{
    const MAST::ElementPropertyCardBase& property = _structural_elem.elem_property();
    
    // make an fe and quadrature object for the requested order for integrating
    // transverse shear
    
    std::auto_ptr<FEBase> fe;
    std::auto_ptr<QBase> qrule;
    FEType fe_type = _fe.get_fe_type();
    
    fe.reset(FEBase::build(_elem.dim(), fe_type).release());
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
    const std::vector<std::vector<Real> >& phi = fe->get_phi();
    const std::vector<Real>& JxW = fe->get_JxW();
    const std::vector<Point>& xyz = fe->get_xyz();
    
    const unsigned int n_phi = (unsigned int)phi.size(), n2 = 6*n_phi;
    DenseVector<Real> phi_vec; phi_vec.resize(n_phi);
    
    DenseVector<Real> tmp_vec3_n2, tmp_vec4_2, tmp_vec5_2;
    DenseMatrix<Real> material_trans_shear_mat, tmp_mat2_n2n2, tmp_mat4_2n2;
    
    tmp_vec3_n2.resize(n2); tmp_vec4_2.resize(2); tmp_vec5_2.resize(2);
    tmp_mat2_n2n2.resize(n2, n2); tmp_mat4_2n2.resize(2, n2);
    
    
    FEMOperatorMatrix Bmat_trans;
    Bmat_trans.reinit(2, 6, n_phi); // only two shear stresses
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        // if temperature is specified, the initialize it to the current location
        //        if (_temperature)
        //            _temperature->initialize(xyz[qp]);
        
        if (!sens_params)
            property.calculate_matrix(_elem,
                                      MAST::SECTION_INTEGRATED_MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX,
                                      material_trans_shear_mat);
        else
            property.calculate_matrix_sensitivity(_elem,
                                                  MAST::SECTION_INTEGRATED_MATERIAL_THERMAL_EXPANSION_MATRIX,
                                                  material_trans_shear_mat,
                                                  *sens_params);
        
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
    
    
    // if sensitivity was requested, calculate the sensitivity wrt temperature
    // next, calculate the sensitivity due to temperature, if it is provided
    //    if (!_temperature || !sens_params)
    //        return;
    //
    //    Real dTemp_dparam = _temperature->sensitivity(*sens_params);
    //
    //    if (dTemp_dparam == 0.) // no need to do calculations if the coefficient is zero
    //        return;
    //
    //    MAST::SensitivityParameters temp_param;
    //    temp_param.add_parameter(_temperature, 1);
    //
    //    // first calculate the sensitivity due to the parameter
    //    for (unsigned int qp=0; qp<JxW.size(); qp++) {
    //
    //        // if temperature is specified, the initialize it to the current location
    //        if (_temperature)
    //            _temperature->initialize(xyz[qp]);
    //
    //        // get the material matrix
    //        property.calculate_matrix_sensitivity(_elem,
    //                                              MAST::SECTION_INTEGRATED_MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX,
    //                                              material_trans_shear_mat,
    //                                              temp_param);
    //
    //    // initialize the strain operator
    //    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
    //        phi_vec(i_nd) = dphi[i_nd][qp](0);  // dphi/dx
    //    Bmat_trans.set_shape_function(0, 1, phi_vec); // gamma-xy:  v
    //    Bmat_trans.set_shape_function(1, 2, phi_vec); // gamma-xz : w
    //
    //    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
    //        phi_vec(i_nd) = phi[i_nd][qp];  // phi
    //    Bmat_trans.set_shape_function(1, 4, phi_vec); // gamma-xy:  thetay
    //    phi_vec.scale(-1.0);
    //    Bmat_trans.set_shape_function(0, 5, phi_vec); // gamma-xz : thetaz
    //
    //
    //        // now add the transverse shear component
    //        Bmat_trans.vector_mult(tmp_vec4_2, local_solution);
    //        material_trans_shear_mat.vector_mult(tmp_vec5_2, tmp_vec4_2);
    //        Bmat_trans.vector_mult_transpose(tmp_vec3_n2, tmp_vec5_2);
    //        local_f.add(-JxW[qp]*dTemp_dparam, tmp_vec3_n2);
    //
    //        if (request_jacobian) {
    //            // now add the transverse shear component
    //            Bmat_trans.left_multiply(tmp_mat4_2n2, material_trans_shear_mat);
    //            Bmat_trans.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat4_2n2);
    //            local_jac.add(-JxW[qp]*dTemp_dparam, tmp_mat2_n2n2);
    //        }
    //    }
}


#endif // __MAST_timoshenko_bending_operator_h__
