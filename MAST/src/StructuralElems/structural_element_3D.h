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
         *   returns a constant reference to the element in local coordinate system. 
         *   For a 3D element, the two are same.
         */
        virtual const MAST::LocalElemBase& local_elem() const {
            libmesh_error(); // this is not defined for a 3D element
        }

        /*!
         *    Local elements are defined for 1D and 2D elements that exist in
         *    3D space. These elements have a local coordinate system associated
         *    with the local coordinate. This method accepts the point defined
         *    the local coordinate system as the input and maps it to the
         *    global coordinate system.
         *
         *    For 3D elements, it is assumed that the local and global 
         *    coordinates are the same.
         */
        virtual void global_coordinates(const Point& local,
                                        Point& global) const {
            global = local;
        }

        /*!
         *   returns a constant reference to the element in local coordinate system.
         *   For a 3D element, the two are same.
         */
        virtual const Elem& get_elem_for_quadrature() const {
            return _elem;
        }

        /*!
         *    Calculates the internal force vector and Jacobian due to
         *    strain energy
         */
        virtual bool internal_force(bool request_jacobian,
                                    DenseVector<Real>& f,
                                    DenseMatrix<Real>& jac,
                                    bool if_ignore_ho_jac);

        
        /*!
         *    Calculates the sensitivity of internal force vector and
         *    Jacobian due to strain energy
         */
        virtual bool internal_force_sensitivity(bool request_jacobian,
                                                DenseVector<Real>& f,
                                                DenseMatrix<Real>& jac,
                                                bool if_ignore_ho_jac);

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
        
        
        /*!
         *   returns the value of maximum von Mises stress over the element
         */
        virtual Real max_von_mises_stress(){
            libmesh_error();
        }
        
        
        /*!
         *   returns the sensitivity of maximum von Mises stress over the element
         */
        virtual Real max_von_mises_stress_sensitivity(){
            libmesh_error();
        }

    protected:
        
        /*!
         *    Calculates the force vector and Jacobian due to thermal stresses
         */
        virtual bool thermal_force(bool request_jacobian,
                                   DenseVector<Real>& f,
                                   DenseMatrix<Real>& jac,
                                   MAST::BoundaryCondition& p);

        /*!
         *    Calculates the sensitivity of force vector and Jacobian due to 
         *    thermal stresses
         */
        virtual bool thermal_force_sensitivity(bool request_jacobian,
                                               DenseVector<Real>& f,
                                               DenseMatrix<Real>& jac,
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
