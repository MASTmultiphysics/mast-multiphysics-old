//
//  element_property_card_1D.h
//  MAST
//
//  Created by Manav Bhatia on 10/25/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_element_property_card_1D_h__
#define __MAST_element_property_card_1D_h__


// MAST includes
#include "PropertyCards/element_property_card_base.h"
#include "PropertyCards/material_property_card_base.h"


namespace MAST
{
    class ElementPropertyCard1D: public MAST::ElementPropertyCardBase {
        
    public:
        ElementPropertyCard1D(unsigned int pid):
        MAST::ElementPropertyCardBase(pid),
        _bending_model(MAST::DEFAULT_BENDING)
        { }
        
        /*!
         *   virtual destructor
         */
        virtual ~ElementPropertyCard1D() { }
        
        
        /*!
         *   dimension of the element for which this property is defined
         */
        virtual unsigned int dim() const {
            return 1;
        }

        
        /*!
         *   returns the bending model to be used for the 2D element
         */
        void set_bending_model(MAST::BendingOperatorType b)  {
            _bending_model = b;
        }
        
        
        /*!
         *   returns the bending model to be used for the 2D element.
         */
        virtual MAST::BendingOperatorType bending_model(const Elem& elem,
                                                 const FEType& fe) const;
        
        
        /*!
         *    returns the extra quadrature order (on top of the system) that
         *    this element should use. This is elevated by two orders for a DKT
         *    element
         */
        virtual unsigned int extra_quadrature_order(const Elem& elem,
                                                    const FEType& fe) const {
            if (this->bending_model(elem, fe) == MAST::BERNOULLI)
                return 2;
            else
                return 0;
        }
        
        
        /*!
         *   returns value of the property \par val. The string values for 
         *   \par val are IYY, IZZ, IYZ
         */
        virtual Real value(const std::string& val) const = 0;
        
        /*!
         *   vector in the x-y plane of the element. This should not be the same
         *   as the element x-axis.
         */
        Point& y_vector() {
            return _local_y;
        }
        

        /*!
         *   constant reference to vector in the x-y plane of the element. 
         *   This should not be the same as the element x-axis.
         */
        const Point& y_vector() const {
            return _local_y;
        }
        
        
    protected:
        
        /*!
         *    initializes the vector to the prestress in the element
         */
        virtual void _prestress_vector(const DenseMatrix<Real>& T,
                                       DenseVector<Real>& v) const {
            v.resize(2); // zero, if the stress has not been defined

            if (_prestress.m() != 0) {
                
                DenseMatrix<Real> mat;
                _prestress_matrix(T, mat);
                v(0) = mat(0,0); // sigma_xx
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
                
                m(0,0) = mat(0,0);
            }
        }

        /*!
         *   material property card. By default this chooses DKT for 3 noded
         *   triangles and Mindling for all other elements
         */
        MAST::BendingOperatorType _bending_model;
        
        /*!
         *   vector in the x-y plane.
         */
        Point _local_y;
        
    };
    
    
}



inline
MAST::BendingOperatorType
MAST::ElementPropertyCard1D::bending_model(const Elem& elem,
                                           const FEType& fe) const {
    // for an EDGE2 element, default bending is Bernoulli. For all other elements
    // the default is Timoshenko. Otherwise it returns the model set for
    // this card.
    switch (elem.type()) {
        case EDGE2:
            if ((fe.family == LAGRANGE) &&
                (fe.order  == FIRST) &&
                (_bending_model == MAST::DEFAULT_BENDING))
                return MAST::BERNOULLI;
            else
                return MAST::TIMOSHENKO;
            break;
            
        default:
            if (_bending_model == MAST::DEFAULT_BENDING)
                return MAST::TIMOSHENKO;
            else
                return _bending_model;
            break;
    }
}


#endif  // __MAST_element_property_card_1D_h__
