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
    enum BendingModel2D {
        DKT,
        MINDLIN,
        NO_BENDING_2D,
        DEFAULT_BENDING_2D
    };
    
    
    
    class ElementPropertyCard2D: public MAST::ElementPropertyCardBase {
        
    public:
        ElementPropertyCard2D():
        MAST::ElementPropertyCardBase(),
        _bending_model(MAST::DEFAULT_BENDING_2D),
        _if_plane_stress(true)
        { }
        
        /*!
         *   virtual destructor
         */
        ~ElementPropertyCard2D() { }
        
        /*!
         *   returns the bending model to be used for the 2D element
         */
        void set_bending_model(MAST::BendingModel2D b)  {
            _bending_model = b;
        }
        
        
        /*!
         *   returns the bending model to be used for the 2D element.
         */
        MAST::BendingModel2D bending_model(const Elem& elem,
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
        virtual void _prestress_vector(DenseVector<Real>& v) const {
            if (_prestress.size() == 0)
                v.resize(3); // zero, if the stress has not been defined
            else {
                v.resize(3);
                v(0) = _prestress(0);
                v(1) = _prestress(1);
                v(2) = _prestress(3);
            }
        }
        
        
        /*!
         *    initializes the matrix to the prestress in the element
         */
        virtual void _prestress_matrix(DenseMatrix<Real>& m) const {
            m.resize(2, 2);
            if (_prestress.size() == 6) {
                for (unsigned int i=0; i<2; i++)
                    m(i,i) = _prestress(i);
                m(0,1) = _prestress(3);  // tau_xy
                m(1,0) = _prestress(3);  // tau_xy
            }
        }
        

        
        /*!
         *   material property card. By default this chooses DKT for 3 noded
         *   triangles and Mindling for all other elements
         */
        MAST::BendingModel2D _bending_model;
        
        /*!
         *   if the analysis is plne stress, otherwise it is plane strain.
         *   Note that this is true by default
         */
        bool _if_plane_stress;
        
    };
    
    
    
    
    class Solid2DSectionElementPropertyCard : public MAST::ElementPropertyCard2D {
    public:
        Solid2DSectionElementPropertyCard():
        MAST::ElementPropertyCard2D(),
        _material(NULL)
        { }
        
        /*!
         *   virtual destructor
         */
        ~Solid2DSectionElementPropertyCard() { }
        
        /*!
         *    sets the material card
         */
        void set_material(MAST::MaterialPropertyCardBase& mat) {
            _material = &mat;
        }
        
        
        /*!
         *    returns a reference to the material
         */
        const MAST::MaterialPropertyCardBase& get_material() const {
            libmesh_assert(_material); // make sure it has already been set
            return *_material;
        }
        

        /*!
         *   calculates the material matrix in \par m of type \par t.
         */
        virtual void calculate_matrix(const Elem& elem,
                                      MAST::ElemenetPropertyMatrixType t,
                                      DenseMatrix<Real>& m) const;

        /*!
         *    initializes the vector to the prestress in the element
         */
        virtual void prestress_vector(MAST::ElemenetPropertyMatrixType t,
                                      DenseVector<Real>& v) const;
        
        
        /*!
         *    initializes the matrix to the prestress in the element
         */
        virtual void prestress_matrix(MAST::ElemenetPropertyMatrixType t,
                                      DenseMatrix<Real>& m) const;
        

    protected:
        
        /*!
         *   material property card
         */
        MAST::MaterialPropertyCardBase *_material;
    };
}



inline
MAST::BendingModel2D
MAST::ElementPropertyCard2D::bending_model(const Elem& elem,
                                           const FEType& fe) const {
    // for a TRI3 element, default bending is DKT. For all other elements
    // the default is Mindlin. Otherwise it returns the model set for
    // this card.
    switch (elem.type()) {
        case TRI3:
            if ((fe.family == LAGRANGE) &&
                (fe.order  == FIRST) &&
                (_bending_model == MAST::DEFAULT_BENDING_2D))
                return MAST::DKT;
            else
                return _bending_model;
            break;
            
        default:
            if (_bending_model == MAST::DEFAULT_BENDING_2D)
                return MAST::MINDLIN;
            else
                return _bending_model;
            break;
    }
}



inline void
MAST::Solid2DSectionElementPropertyCard::calculate_matrix(const libMesh::Elem &elem,
                                                          MAST::ElemenetPropertyMatrixType t,
                                                          DenseMatrix<Real>& m) const
{
    libmesh_assert(_material); // should have been set
    
    switch (elem.dim()) {

        case 2: {
            double h = this->get<Real>("h")();
            switch (t) {
                case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX:
                    _material->calculate_2d_matrix(MAST::MATERIAL_STIFFNESS_MATRIX,
                                                   m, _if_plane_stress);
                    m.scale(h);
                    break;
                    
                case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_B_MATRIX:
                    // for solid sections with isotropic material this is zero
                    m.resize(3,3);
                    break;
                    
                case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_D_MATRIX:
                    _material->calculate_2d_matrix(MAST::MATERIAL_STIFFNESS_MATRIX,
                                                   m, _if_plane_stress);
                    m.scale(pow(h,3)/12.);
                    break;

                case MAST::SECTION_INTEGRATED_MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX:
                    _material->calculate_2d_matrix(MAST::MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX,
                                                   m, _if_plane_stress);
                    m.scale(h);
                    break;
                    
                case MAST::SECTION_INTEGRATED_MATERIAL_INERTIA_MATRIX:
                    _material->calculate_2d_matrix(MAST::MATERIAL_INERTIA_MATRIX,
                                                   m, _if_plane_stress);
                    m.scale(h);
                    // now scale the rotation dofs with small factors
                    for (unsigned int i=0; i<3; i++) {
                        m(i+3, i+3) *= 1.0e-6;
                    }
                    break;

                default:
                    libmesh_error();
                    break;
            }
        }
            break;
        
        case 1:
        case 3:
        default:
            libmesh_error();
            break;
    }
}






inline void
MAST::Solid2DSectionElementPropertyCard::prestress_vector(MAST::ElemenetPropertyMatrixType t,
                                                          DenseVector<Real>& v) const {
    double h = this->get<Real>("h")();
    switch (t) {
        case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX:
            _prestress_vector(v);
            v.scale(h);
            break;
                    
        case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_B_MATRIX:
            // for solid sections with isotropic material this is zero
            v.resize(3);
            break;
                    
        default:
            libmesh_error();
            break;
    }
}





inline void
MAST::Solid2DSectionElementPropertyCard::prestress_matrix(MAST::ElemenetPropertyMatrixType t,
                                                          DenseMatrix<Real>& m) const {
    double h = this->get<Real>("h")();
    switch (t) {
        case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX:
            _prestress_matrix(m);
            m.scale(h);
            break;
            
        case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_B_MATRIX:
            // for solid sections with isotropic material this is zero
            m.resize(2,2);
            break;
            
        default:
            libmesh_error();
            break;
    }
}



#endif // __MAST_element_property_card_2D_h__
