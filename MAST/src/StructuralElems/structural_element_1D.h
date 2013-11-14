//
//  structural_element_1D.h
//  MAST
//
//  Created by Manav Bhatia on 10/21/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_structural_element_1D_h__
#define __MAST_structural_element_1D_h__

// C++ includes
#include <memory>


// MAST includes
#include "StructuralElems/structural_elem_base.h"

// Forward declerations
class FEMOperatorMatrix;


namespace MAST {
    
    // Forward declerations
    class BernoulliElem;
    
    /*!
     *   class provides a simple mechanism to create a geometric element
     *   in the local coordinate system.
     */
    class Local1DELem {
    public:
        Local1DELem(const Elem& elem):
        _elem(elem),
        _local_elem(NULL)
        {
            this->create_local_elem();
        }
        
        
        /*!
         *   returns a constant reference to the global element.
         */
        const Elem& global_elem() const {
            return _elem;
        }
        
        
        /*!
         *   returns a constant reference to the local element.
         */
        const Elem& local_elem() const {
            return *_local_elem;
        }
        
    protected:
        
        void create_local_elem() {
            libmesh_error(); // to be implemented
        }
        
        /*!
         *   given element in global coordinate system
         */
        const Elem& _elem;
        
        /*!
         *   element created in local coordinate system
         */
        Elem* _local_elem;
        
    };
    
    
    class StructuralElement1D: public MAST::StructuralElementBase {
        
    public:
        StructuralElement1D(System& sys,
                            const Elem& elem,
                            const MAST::ElementPropertyCardBase& p);
        
        /*!
         *    Calculates the internal force vector and Jacobian due to
         *    strain energy
         */
        virtual bool internal_force(bool request_jacobian,
                                    DenseVector<Real>& f,
                                    DenseMatrix<Real>& jac);
        
        
    protected:
        
        /*!
         *    Calculates the force vector and Jacobian due to thermal stresses
         */
        virtual bool thermal_force(bool request_jacobian,
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
         *   initialize membrane strain operator matrix
         */
        void initialize_extension_strain_operator(const unsigned int qp,
                                                  FEMOperatorMatrix& Bmat);
        
        /*!
         *   initialze the bending strain operator for Mindling plate model.
         *   The bending part is included in \par Bend_par and the transverse
         *   shear part is included in \par Bmat_trans;
         */
        void initialize_timoshenko_bending_strain_operator(const unsigned int qp,
                                                           FEMOperatorMatrix& Bmat_bend,
                                                           FEMOperatorMatrix& Bmat_trans);
        
        
        /*!
         *   initialze the von Karman strain in \par vK_strain, the operator
         *   matrices needed for Jacobian calculation.
         *   vk_strain = [dw/dx 0; 0 dw/dy; dw/dy dw/dx]
         *   Bmat_vk   = [dw/dx; dw/dy]
         */
        void initialize_von_karman_strain_operator(const unsigned int qp,
                                                   DenseVector<Real>& vk_strain,
                                                   DenseMatrix<Real>& vk_dwdxi_mat,
                                                   FEMOperatorMatrix& Bmat_vk);
        
        /*!
         *    DKT element, initialized if the element needs it
         */
        std::auto_ptr<MAST::BernoulliElem> _dkt_elem;
        
        
        /*!
         *   element in local coordinate system
         */
        MAST::Local1DELem _local_elem;
    };
}


#endif // __MAST_structural_element_1d_h__
