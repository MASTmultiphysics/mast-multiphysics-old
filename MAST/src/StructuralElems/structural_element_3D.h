//
//  structural_element_3d.h
//  MAST
//
//  Created by Manav Bhatia on 10/21/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_structural_element_3d_h__
#define __MAST_structural_element_3d_h__

// MAST includes
#include "StructuralElems/structural_elem_base.h"

// Forward declerations
class FEMOperatorMatrix;

namespace MAST {

    // Forward declerations
    class DKTElem;
    class BoundaryCondition;

    
    class StructuralElement3D: public MAST::StructuralElementBase {
        
    public:
        StructuralElement3D(System& sys,
                            const Elem& elem,
                            const MAST::ElementPropertyCardBase& p):
        StructuralElementBase(sys, elem, p)
        {
            _init_fe_and_qrule(elem);
        }
        
        /*!
         *    Calculates the internal force vector and Jacobian due to
         *    strain energy
         */
        virtual bool internal_force(bool request_jacobian,
                                    DenseVector<Real>& f,
                                    DenseMatrix<Real>& jac);

        
        /*!
         *    Calculates the sensitivity of internal force vector and
         *    Jacobian due to strain energy
         */
        virtual bool internal_force_sensitivity(bool request_jacobian,
                                                DenseVector<Real>& f,
                                                DenseMatrix<Real>& jac);

        /*!
         *    Calculates the prestress force vector and Jacobian
         */
        virtual bool prestress_force (bool request_jacobian,
                                      DenseVector<Real>& f,
                                      DenseMatrix<Real>& jac)
        { libmesh_error();}

        
        /*!
         *    Calculates the sensitivity prestress force vector and Jacobian
         */
        virtual bool prestress_force_sensitivity (bool request_jacobian,
                                                  DenseVector<Real>& f,
                                                  DenseMatrix<Real>& jac)
        { libmesh_error();}

    protected:
        
        /*!
         *    Calculates the force vector and Jacobian due to thermal stresses
         */
        virtual bool thermal_force(bool request_jacobian,
                                   DenseVector<Real>& f,
                                   DenseMatrix<Real>& jac);

        /*!
         *    Calculates the sensitivity of force vector and Jacobian due to 
         *    thermal stresses
         */
        virtual bool thermal_force_sensitivity(bool request_jacobian,
                                               DenseVector<Real>& f,
                                               DenseMatrix<Real>& jac);

        /*!
         *    Calculates the force vector and Jacobian due to surface pressure.
         */
        virtual bool surface_pressure_force(bool request_jacobian,
                                            DenseVector<Real>& f,
                                            DenseMatrix<Real>& jac,
                                            const unsigned int side,
                                            MAST::BoundaryCondition& p);
        
        
        /*!
         *    Calculates the force vector and Jacobian due to surface pressure.
         */
        virtual bool surface_pressure_force_sensitivity(bool request_jacobian,
                                                        DenseVector<Real>& f,
                                                        DenseMatrix<Real>& jac,
                                                        const unsigned int side,
                                                        MAST::BoundaryCondition& p);
        
        /*!
         *   initialize strain operator matrix
         */
        void initialize_strain_operator (const unsigned int qp,
                                         FEMOperatorMatrix& Bmat);
        
        /*!
         *   matrix that transforms the global dofs to the local element coordinate
         *   system
         */
        virtual const DenseMatrix<Real>& _transformation_matrix() const {
            libmesh_error(); // should not be called for a 3D elem
        }
    };
}


#endif // __MAST_structural_element_3d_h__
