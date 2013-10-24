//
//  element_property_card_3D.h
//  MAST
//
//  Created by Manav Bhatia on 10/24/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_element_property_card_3D_h__
#define __MAST_element_property_card_3D_h__

// MAST includes
#include "PropertyCards/element_property_card_base.h"
#include "PropertyCards/material_property_card_base.h"


namespace MAST
{
    class ElementPropertyCard3D: public MAST::ElementPropertyCardBase {
        
    public:
    public:
        ElementPropertyCard3D():
        MAST::ElementPropertyCardBase(),
        _material(NULL)
        { }
        
        /*!
         *   virtual destructor
         */
        virtual ~ElementPropertyCard3D() { }
        
        /*!
         *   calculates the material matrix in \par m of type \par t.
         */
        virtual void calculate_matrix(const Elem& elem,
                                      MAST::ElemenetPropertyMatrixType t,
                                      DenseMatrix<Real>& m) const;
        
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
        
    protected:
        
        /*!
         *    pointer to the material property card
         */
        MAST::MaterialPropertyCardBase* _material;
    };
}


inline void
MAST::ElementPropertyCard3D::calculate_matrix(const libMesh::Elem &elem,
                                              MAST::ElemenetPropertyMatrixType t,
                                              DenseMatrix<Real>& m) const
{
    libmesh_assert(_material); // should have been set
    
    switch (elem.dim()) {
        case 3:
            switch (t) {
                case MAST::SECTION_INTEGRATED_MATERIAL_A_MATRIX:
                    _material->calculate_3d_matrix(MAST::MATERIAL_STIFFNESS_MATRIX, m);
                    break;
                    
                case MAST::SECTION_INTEGRATED_MATERIAL_INERTIA_MATRIX:
                    _material->calculate_3d_matrix(MAST::MATERIAL_INERTIA_MATRIX, m);
                    break;

                case MAST::SECTION_INTEGRATED_MATERIAL_THERMAL_EXPANSION:
                    _material->calculate_3d_matrix(MAST::MATERIAL_THERMAL_EXPANSION_MATRIX, m);
                    break;

            }
            break;
            
        case 1:
        case 2:
        default:
            libmesh_error();
            break;
    }
}


#endif // __MAST_element_property_card_3D_h__
