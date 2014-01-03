//
//  structural_element_2D.h
//  MAST
//
//  Created by Manav Bhatia on 10/21/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_structural_element_2D_h__
#define __MAST_structural_element_2D_h__

// C++ includes
#include <memory>


// MAST includes
#include "StructuralElems/bending_structural_elem.h"


// Forward declerations
class FEMOperatorMatrix;


namespace MAST {

    // Forward declerations
    class BoundaryCondition;
    
    /*!
     *   class provides a simple mechanism to create a geometric element
     *   in the local coordinate system.
     */
    class Local2DElem : public MAST::LocalElemBase {
    public:
        Local2DElem(const Elem& elem):
        MAST::LocalElemBase(elem)
        {
            _create_local_elem();
        }
        
        
        virtual ~Local2DElem() {
            // the local element may not have been created
            // for cases where the original element lies in the xy-plane
            if (_local_elem) {
                delete _local_elem;
                for (unsigned int i=0; i<_local_nodes.size(); i++)
                    delete _local_nodes[i];
            }
        }
        
        
    protected:
        
        /*!
         *    creation of an element in the local coordinate system
         */
        void _create_local_elem();
    };
    
    
    
    class StructuralElement2D: public MAST::BendingStructuralElem {
        
    public:
        StructuralElement2D(System& sys,
                            const Elem& elem,
                            const MAST::ElementPropertyCardBase& p);


        /*!
         *    row dimension of the direct strain matrix, also used for the
         *    bending operator row dimension
         */
        virtual unsigned int n_direct_strain_components() {
            return 3;
        }
        
        /*!
         *    row dimension of the von Karman strain matrix
         */
        virtual unsigned int n_von_karman_strain_components() {
            return 2;
        }

        /*!
         *   returns the value of maximum von Mises stress over the element
         */
        virtual Real max_von_mises_stress();
        
        
        /*!
         *   returns the sensitivity of maximum von Mises stress over the element
         */
        virtual Real max_von_mises_stress_sensitivity() {
            libmesh_error();
        }

    protected:


        /*!
         *   initialize membrane strain operator matrix
         */
        virtual void initialize_direct_strain_operator(const unsigned int qp,
                                                       FEMOperatorMatrix& Bmat);
        
        /*!
         *   initialze the von Karman strain in \par vK_strain, the operator
         *   matrices needed for Jacobian calculation.
         *   vk_strain = [dw/dx 0; 0 dw/dy; dw/dy dw/dx]
         *   Bmat_vk   = [dw/dx; dw/dy]
         */
        virtual void initialize_von_karman_strain_operator(const unsigned int qp,
                                                           DenseVector<Real>& vk_strain,
                                                           DenseMatrix<Real>& vk_dwdxi_mat,
                                                           FEMOperatorMatrix& Bmat_vk);
    };
}



#endif  // __MAST_structural_element_2D_h__
