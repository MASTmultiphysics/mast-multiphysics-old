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
    
    class Solid2DSectionElementPropertyCard : public MAST::ElementPropertyCardBase {
        
        Solid2DSectionElementPropertyCard():
        MAST::ElementPropertyCardBase(),
        _material(NULL),
        _if_plane_stress(false)
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
         *   material property card
         */
        MAST::MaterialPropertyCardBase *_material;
        
        /*!
         *   if the analysis is plne stress, otherwise it is plane strain.
         *   Note that this is true by default
         */
        bool _if_plane_stress;
        
    };
    
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
                                                   m,
                                                   _if_plane_stress);
                    m.scale(h);
                    break;
                    
                case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_B_MATRIX:
                    // for solid sections with isotropic material this is zero
                    break;
                    
                case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_D_MATRIX:
                    _material->calculate_2d_matrix(MAST::MATERIAL_STIFFNESS_MATRIX,
                                                   m,
                                                   _if_plane_stress);
                    m.scale(pow(h,3)/12.);
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





#endif // __MAST_element_property_card_2D_h__
