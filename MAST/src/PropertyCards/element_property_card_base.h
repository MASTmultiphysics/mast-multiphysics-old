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
    enum ElemenetPropertyMatrixType {
        STRESS_MATERIAL_MATRIX,
        DAMPING_MATERIAL_MATRIX,
        INERTIA_MATERIAL_MATRIX,
        SECTION_INTEGRATED_STRESS_MATERIAL_MATRIX,
        SECTION_INTEGRATED_DAMPING_MATERIAL_MATRIX,
        SECTION_INTEGRATED_INERTIA_MATERIAL_MATRIX,
        SOLID_MATERIAL_MATRIX
    };
    
    
    class ElementPropertyCardBase: public MAST::PropertyCardBase {
        
    public:
        ElementPropertyCardBase():
        MAST::PropertyCardBase()
        { }
        

        /*!
         *   calculates the material matrix in \par m of type \par t
         */
        virtual void calculate_matrix(const Elem& elem,
                                      MAST::ElemenetPropertyMatrixType t,
                                      DenseMatrix<Real>& m) const;
        
        
    protected:

        /*!
         *    pointer to the material property card
         */
        MAST::MaterialPropertyCardBase* _material;
        
    };
}



#endif // __MAST_element_property_card_base_h__
