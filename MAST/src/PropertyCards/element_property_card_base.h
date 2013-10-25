//
//  element_property_card_base.h
//  MAST
//
//  Created by Manav Bhatia on 10/15/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_element_property_card_base_h__
#define __MAST_element_property_card_base_h__

// MAST includes
#include "PropertyCards/property_card_base.h"

// libMesh includes
#include "libmesh/elem.h"


namespace MAST
{
    enum StrainType {
        LINEAR_STRAIN,
        VON_KARMAN_STRAIN
    };
    
    enum ElemenetPropertyMatrixType {
        SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX,
        SECTION_INTEGRATED_MATERIAL_STIFFNESS_B_MATRIX_1D_TRANSVERSE,
        SECTION_INTEGRATED_MATERIAL_STIFFNESS_D_MATRIX_1D_TRANSVERSE,
        SECTION_INTEGRATED_MATERIAL_STIFFNESS_B_MATRIX_1D_CHORDWISE,
        SECTION_INTEGRATED_MATERIAL_STIFFNESS_D_MATRIX_1D_CHORDWISE,
        SECTION_INTEGRATED_MATERIAL_STIFFNESS_B_MATRIX,
        SECTION_INTEGRATED_MATERIAL_STIFFNESS_D_MATRIX,
        SECTION_INTEGRATED_MATERIAL_DAMPING_MATRIX,
        SECTION_INTEGRATED_MATERIAL_INERTIA_MATRIX,
        SECTION_INTEGRATED_MATERIAL_THERMAL_EXPANSION_MATRIX,
        SECTION_INTEGRATED_MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX
    };
    
    
    class ElementPropertyCardBase: public MAST::PropertyCardBase {
        
    public:
        ElementPropertyCardBase():
        MAST::PropertyCardBase(),
        _strain_type(MAST::LINEAR_STRAIN)
        { }
        
        /*!
         *   virtual destructor
         */
        virtual ~ElementPropertyCardBase() { }
        
        /*!
         *    returns the extra quadrature order (on top of the system) that 
         *    this element should use. By default this is zero, and can be
         *    changed by the derived classes
         */
        virtual unsigned int extra_quadrature_order(const Elem& elem) const {
            return 0;
        }
        
        /*!
         *   calculates the material matrix in \par m of type \par t.
         */
        virtual void calculate_matrix(const Elem& elem,
                                      MAST::ElemenetPropertyMatrixType t,
                                      DenseMatrix<Real>& m) const = 0;
        
        /*!
         *    sets the type of strain to be used, which is LINEAR_STRAIN by
         *    default
         */
        void set_strain(MAST::StrainType strain) {
            _strain_type = strain;
        }
        
        
        /*!
         *    returns the type of strain to be used for this element
         */
        const MAST::StrainType strain_type() const {
            return _strain_type;
        }
        
        
    protected:
        
        /*!
         *    type of nonlinear strain to be used for analysis
         */
        MAST::StrainType _strain_type;
    };
    
    
}



#endif // __MAST_element_property_card_base_h__
