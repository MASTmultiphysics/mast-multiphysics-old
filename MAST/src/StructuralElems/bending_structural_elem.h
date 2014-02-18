//
//  bending_structural_elem_base.h
//  MAST
//
//  Created by Manav Bhatia on 11/27/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_bending_structural_elem_base_h__
#define __MAST_bending_structural_elem_base_h__

// C++ includes
#include <memory>


// MAST includes
#include "StructuralElems/structural_elem_base.h"
#include "Numerics/function_base.h"
#include "Numerics/constant_function.h"

// Forward declerations
class FEMOperatorMatrix;


namespace MAST {
    
    // Forward declerations
    class BendingOperator;
    class BoundaryCondition;
    
    /*!
     *   class provides a simple mechanism to create a geometric element
     *   in the local coordinate system.
     */
    class LocalElemBase {
    public:
        LocalElemBase(const Elem& elem):
        _elem(elem),
        _local_elem(NULL)
        { }
        
        
        virtual ~LocalElemBase() { }
        
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
            if (!_local_elem) // original element lies in the xy-plane
                return _elem;
            else
                return *_local_elem;
        }
        
        /*!
         *    returns the transformation matrix for this element. This is used
         *    to map the coordinates from local to global coordinate system
         */
        const DenseMatrix<Real>& T_matrix() const {
            return _T_mat;
        }

        
        /*!
         *    returns the transformation matrix for this element. This is used
         *    to map the coordinates from local to global coordinate system
         */
        std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>> T_matrix_function() const {
            return std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>>
            (new MAST::ConstantFunction<DenseMatrix<Real> >("T_mat", _T_mat));
        }

        
        /*!
         *   maps the local coordinates to the global coordinates
         */
        void global_coordinates(const Point& global,
                                Point& local) const {
            local = 0.;
            
            // now calculate the global coordinates with respect to the origin
            for (unsigned int j=0; j<3; j++)
                for (unsigned int k=0; k<3; k++)
                    local(j) += _T_mat(j,k)*global(k);
            
            // shift to the global coordinate
            local += (*_elem.get_node(0));
        }
        

        protected:

        /*!
         *    creation of an element in the local coordinate system
         */
        void _create_local_elem();
        
        /*!
         *   given element in global coordinate system
         */
        const Elem& _elem;
        
        /*!
         *   element created in local coordinate system
         */
        Elem* _local_elem;
        
        /*!
         *   nodes for local element
         */
        std::vector<Node*> _local_nodes;
        
        /*!
         *    Transformation matrix defines T_ij = V_i^t . Vn_j, where
         *    V_i are the unit vectors of the global cs, and Vn_j are the
         *    unit vectors of the local cs. To transform a vector from global to
         *    local cs,    an_j = T^t a_i, and the reverse transformation is
         *    obtained as  a_j  = T  an_i
         */
        DenseMatrix<Real> _T_mat;
    };
    
    
    /*!
     *    this method creates a local element object and returns it in a 
     *    smart pointer
     */
    std::auto_ptr<MAST::LocalElemBase>
    build_local_elem(StructuralElementBase& elem);
    
    
    
    class BendingStructuralElem: public MAST::StructuralElementBase {
        
    public:
        BendingStructuralElem(System& sys,
                            const Elem& elem,
                            const MAST::ElementPropertyCardBase& p);
        
        /*!
         *   returns a constant reference to the element in local coordinate system
         */
        virtual const MAST::LocalElemBase& local_elem() const {
            return *_local_elem;
        }

        /*!
         *    Local elements are defined for 1D and 2D elements that exist in
         *    3D space. These elements have a local coordinate system associated
         *    with the local coordinate. This method accepts the point defined
         *    the local coordinate system as the input and maps it to the
         *    global coordinate system.
         *
         *    For 1D and 2D elements the local element is asked to perform the 
         *    mapping
         */
        virtual void global_coordinates(const Point& local,
                                        Point& global) const {
            _local_elem->global_coordinates(local, global);
        }

        
        /*!
         *   returns a constant reference to the element in local coordinate system
         */
        virtual const Elem& get_elem_for_quadrature() const {
            return _local_elem->local_elem();
        }

        /*!
         *    row dimension of the direct strain matrix, also used for the 
         *    bending operator row dimension
         */
        virtual unsigned int n_direct_strain_components() = 0;
        
        /*!
         *    row dimension of the von Karman strain matrix
         */
        virtual unsigned int n_von_karman_strain_components() = 0;
        
    protected:
        
        
        /*!
         *   initialize membrane strain operator matrix
         */
        virtual void initialize_direct_strain_operator(const unsigned int qp,
                                                       FEMOperatorMatrix& Bmat) = 0;
        
        
        /*!
         *   matrix that transforms the global dofs to the local element coordinate
         *   system
         */
        virtual const DenseMatrix<Real>& _transformation_matrix() const {
            return _local_elem->T_matrix();
        }
        
        /*!
         *    bending operator used for this elmeent
         */
        std::auto_ptr<MAST::BendingOperator> _bending_operator;
        
        
        /*!
         *   element in local coordinate system
         */
        std::auto_ptr<MAST::LocalElemBase> _local_elem;
        
    };
}


#endif  //__MAST_bending_structural_elem_base_h__
