//
//  element_property_card_2D.h
//  MAST
//
//  Created by Manav Bhatia on 10/23/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_element_property_card_2D_h__
#define __MAST_element_property_card_2D_h__

// MAST includes
#include "PropertyCards/element_property_card_base.h"
#include "PropertyCards/material_property_card_base.h"


namespace MAST
{
    
    class ElementPropertyCard2D: public MAST::ElementPropertyCardBase {
        
    public:
        ElementPropertyCard2D(unsigned int pid):
        MAST::ElementPropertyCardBase(pid),
        _bending_model(MAST::DEFAULT_BENDING),
        _if_plane_stress(true)
        { }
        
        /*!
         *   virtual destructor
         */
        virtual ~ElementPropertyCard2D() { }
        
        /*!
         *   returns the bending model to be used for the 2D element
         */
        void set_bending_model(MAST::BendingOperatorType b)  {
            _bending_model = b;
        }
        
        
        /*!
         *   returns the bending model to be used for the 2D element.
         */
        MAST::BendingOperatorType bending_model(const Elem& elem,
                                         const FEType& fe) const;
        
        
        /*!
         *    returns the extra quadrature order (on top of the system) that
         *    this element should use. This is elevated by two orders for a DKT
         *    element
         */
        virtual unsigned int extra_quadrature_order(const Elem& elem,
                                                    const FEType& fe) const {
            if (this->bending_model(elem, fe) == MAST::DKT)
                return 2;
            else
                return 0;
        }
        
        /*!
         *   sets the flag for plane stress.
         */
        void set_plane_stress(bool val) {
            _if_plane_stress = val;
        }
        
        /*!
         *   returns the flag for plane stress.
         */
        bool plane_stress() const {
            return _if_plane_stress;
        }
        
        
    protected:
        
        /*!
         *    initializes the vector to the prestress in the element
         */
        virtual void _prestress_vector(const DenseMatrix<Real>& T,
                                       DenseVector<Real>& v) const {
            v.resize(3); // zero, if the stress has not been defined
            if (_prestress.m() != 0) {

                DenseMatrix<Real> mat;
                _prestress_matrix(T, mat);
                v(0) = mat(0,0); // sigma_xx
                v(1) = mat(1,1); // sigma_yy
                v(2) = mat(0,1); // sigma_xy
            }
        }
        
        
        /*!
         *    initializes the matrix to the prestress in the element
         */
        virtual void _prestress_matrix(const DenseMatrix<Real>& T,
                                       DenseMatrix<Real>& m) const {
            m.resize(2, 2);
            if (_prestress.m() != 0) {
                DenseMatrix<Real> mat; mat = _prestress;
                mat.right_multiply(T);
                mat.left_multiply_transpose(T);
                
                for (unsigned int i=0; i<2; i++)
                    for (unsigned int j=0; j<2; j++)
                        m(i,j) = mat(i,j);
            }
        }
        

        
        /*!
         *   material property card. By default this chooses DKT for 3 noded
         *   triangles and Mindling for all other elements
         */
        MAST::BendingOperatorType _bending_model;
        
        /*!
         *   if the analysis is plne stress, otherwise it is plane strain.
         *   Note that this is true by default
         */
        bool _if_plane_stress;
        
    };
    
}

inline
MAST::BendingOperatorType
MAST::ElementPropertyCard2D::bending_model(const Elem& elem,
                                           const FEType& fe) const {
    // for a TRI3 element, default bending is DKT. For all other elements
    // the default is Mindlin. Otherwise it returns the model set for
    // this card.
    switch (elem.type()) {
        case TRI3:
            if ((fe.family == LAGRANGE) &&
                (fe.order  == FIRST) &&
                (_bending_model == MAST::DEFAULT_BENDING))
                return MAST::DKT;
            else
                return MAST::MINDLIN;
            break;
            
        default:
            if (_bending_model == MAST::DEFAULT_BENDING)
                return MAST::MINDLIN;
            else
                return _bending_model;
            break;
    }
}



#endif // __MAST_element_property_card_2D_h__
