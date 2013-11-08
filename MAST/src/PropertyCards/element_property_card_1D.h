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
    enum BendingModel1D {
        BERNOULLI,
        TIMOSHENKO,
        NO_BENDING_1D,
        DEFAULT_BENDING_1D
    };
    
    
    
    class ElementPropertyCard1D: public MAST::ElementPropertyCardBase {
        
    public:
        ElementPropertyCard1D():
        MAST::ElementPropertyCardBase(),
        _bending_model(MAST::DEFAULT_BENDING_1D)
        { }
        
        /*!
         *   virtual destructor
         */
        ~ElementPropertyCard1D() { }
        
        /*!
         *   returns the bending model to be used for the 2D element
         */
        void set_bending_model(MAST::BendingModel1D b)  {
            _bending_model = b;
        }
        
        
        /*!
         *   returns the bending model to be used for the 2D element.
         */
        MAST::BendingModel1D bending_model(const Elem& elem,
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
        
    protected:
        
        /*!
         *   material property card. By default this chooses DKT for 3 noded
         *   triangles and Mindling for all other elements
         */
        MAST::BendingModel1D _bending_model;
        
    };
    
    
    class Solid1DSectionElementPropertyCard : public MAST::ElementPropertyCard1D {
        
        Solid1DSectionElementPropertyCard():
        MAST::ElementPropertyCard1D(),
        _material(NULL)
        { }
        
        /*!
         *   virtual destructor
         */
        ~Solid1DSectionElementPropertyCard() { }
        
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
         *   calculates the matrix in \par m of type \par t.
         */
        virtual void calculate_matrix(const Elem& elem,
                                      MAST::ElemenetPropertyMatrixType t,
                                      DenseMatrix<Real>& m) const;

        /*!
         *   calculates the sensitivity of matrix in \par m of type \par t.
         */
        virtual void calculate_matrix_sensitivity(const Elem& elem,
                                                  MAST::ElemenetPropertyMatrixType t,
                                                  DenseMatrix<Real>& m,
                                                  const MAST::SensitivityParameters& p) const;

    protected:
        
        /*!
         *   material property card
         */
        MAST::MaterialPropertyCardBase *_material;
    };
}



inline
MAST::BendingModel1D
MAST::ElementPropertyCard1D::bending_model(const Elem& elem,
                                           const FEType& fe) const {
    // for a TRI3 element, default bending is DKT. For all other elements
    // the default is Mindlin. Otherwise it returns the model set for
    // this card.
    switch (elem.type()) {
        case EDGE2:
            if ((fe.family == LAGRANGE) &&
                (fe.order  == FIRST) &&
                (_bending_model == MAST::DEFAULT_BENDING_1D))
                return MAST::BERNOULLI;
            else
                return _bending_model;
            break;
            
        default:
            if (_bending_model == MAST::DEFAULT_BENDING_1D)
                return MAST::TIMOSHENKO;
            else
                return _bending_model;
            break;
    }
}



inline void
MAST::Solid1DSectionElementPropertyCard::calculate_matrix(const libMesh::Elem &elem,
                                                          MAST::ElemenetPropertyMatrixType t,
                                                          DenseMatrix<Real>& m) const
{
    libmesh_assert(_material); // should have been set
    
    switch (elem.dim()) {
            
        case 1: {
            double h = this->get<Real>("h")(), // section height
            b = this->get<Real>("b")(),        // section width
            Area = b*h, Itrans = b*pow(h,3)/12., Ichord = h*pow(b,3)/12.;
            switch (t) {
                case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX:
                    _material->calculate_1d_matrix(MAST::MATERIAL_STIFFNESS_MATRIX,
                                                   m);
                    m.scale(Area);
                    break;
                    
                case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_B_MATRIX_1D_TRANSVERSE:
                case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_B_MATRIX_1D_CHORDWISE:
                    // for solid sections with isotropic material this is zero
                    break;
                    
                case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_D_MATRIX_1D_TRANSVERSE:
                    _material->calculate_1d_matrix(MAST::MATERIAL_STIFFNESS_MATRIX,
                                                   m);
                    m.scale(Itrans);
                    break;

                case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_D_MATRIX_1D_CHORDWISE:
                    _material->calculate_1d_matrix(MAST::MATERIAL_STIFFNESS_MATRIX,
                                                   m);
                    m.scale(Ichord);
                    break;

                    
                case MAST::SECTION_INTEGRATED_MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX:
                    _material->calculate_1d_matrix(MAST::MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX,
                                                   m);
                    m.scale(Area);
                    break;
                    
                default:
                    libmesh_error();
                    break;
            }
        }
            break;
            
        case 2:
        case 3:
        default:
            libmesh_error();
            break;
    }
}



inline void
MAST::Solid1DSectionElementPropertyCard::calculate_matrix_sensitivity(const libMesh::Elem &elem,
                                                                      MAST::ElemenetPropertyMatrixType t,
                                                                      DenseMatrix<Real>& m,
                                                                      const MAST::SensitivityParameters& p) const
{
    libmesh_assert(_material); // should have been set

    // only first order sensitivities are calculated at this point
    libmesh_assert_equal_to(p.total_order(), 1);
    
    const MAST::SensitivityParameters::ParameterMap& p_map = p.get_map();
    MAST::SensitivityParameters::ParameterMap::const_iterator it, end;
    it = p_map.begin(); end = p_map.end();
    
    const MAST::FunctionBase& f = *(it->first);

    DenseMatrix<Real> dm;
    
    switch (elem.dim()) {
            
        case 1: {
            double h = this->get<Real>("h")(), // section height
            b = this->get<Real>("b")(),        // section width
            Area = b*h, Itrans = b*pow(h,3)/12., Ichord = h*pow(b,3)/12.,
            dAreadb = h, dAreadh = b,
            dItransdb = pow(h,3)/12., dItransdh = b*pow(h,2)/4.,
            dIchorddb = h*pow(b,2)/4., dIchorddh =  pow(b,3)/12.;
            switch (t) {
                case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX: {
                    _material->calculate_1d_matrix(MAST::MATERIAL_STIFFNESS_MATRIX,
                                                   m);
                    _material->calculate_1d_matrix_sensitivity(MAST::MATERIAL_STIFFNESS_MATRIX,
                                                               dm, p);
                    if (f.name() == "h")
                        m.scale(dAreadh);
                    else if (f.name() == "b")
                        m.scale(dAreadb);
                    else
                        m.scale(0.);
                    
                    m.add(Area, dm);
                }
                    break;
                    
                case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_B_MATRIX_1D_TRANSVERSE:
                case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_B_MATRIX_1D_CHORDWISE:
                    // for solid sections with isotropic material this is zero
                    break;
                    
                case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_D_MATRIX_1D_TRANSVERSE: {
                    _material->calculate_1d_matrix(MAST::MATERIAL_STIFFNESS_MATRIX,
                                                   m);
                    _material->calculate_1d_matrix_sensitivity(MAST::MATERIAL_STIFFNESS_MATRIX,
                                                               dm, p);
                    
                    if (f.name() == "h")
                        m.scale(dItransdh);
                    else if (f.name() == "b")
                        m.scale(dItransdb);
                    else
                        m.scale(0.);
                    
                    m.add(Itrans, dm);
                }
                    break;
                    
                case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_D_MATRIX_1D_CHORDWISE: {
                    _material->calculate_1d_matrix(MAST::MATERIAL_STIFFNESS_MATRIX,
                                                   m);
                    _material->calculate_1d_matrix_sensitivity(MAST::MATERIAL_STIFFNESS_MATRIX,
                                                               dm, p);
                    
                    if (f.name() == "h")
                        m.scale(dIchorddh);
                    else if (f.name() == "b")
                        m.scale(dIchorddb);
                    else
                        m.scale(0.);
                    
                    m.add(Ichord, dm);
                }
                    break;
                    
                    
                case MAST::SECTION_INTEGRATED_MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX: {
                    _material->calculate_1d_matrix(MAST::MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX,
                                                   m);
                    _material->calculate_1d_matrix_sensitivity(MAST::MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX,
                                                               dm, p);
                    if (f.name() == "h")
                        m.scale(dAreadh);
                    else if (f.name() == "b")
                        m.scale(dAreadb);
                    else
                        m.scale(0.);
                    
                    m.add(Area, dm);
                }
                    break;
                    
                default:
                    libmesh_error();
                    break;
            }
        }
            break;
            
        case 2:
        case 3:
        default:
            libmesh_error();
            break;
    }
}


#endif  // __MAST_element_property_card_1D_h__
