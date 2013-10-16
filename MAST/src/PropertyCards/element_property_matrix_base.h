//
//  element_property_matrix_base.h
//  MAST
//
//  Created by Manav Bhatia on 10/15/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_element_property_matrix_base_h__
#define __MAST_element_property_matrix_base_h__

// MAST includes
#include "PropertyCards/property_card_base.h"


namespace MAST
{
    enum ElemenetPropertyMatrixType {
        PLANE_STRESS_MATERIAL_MATRIX,
        PLANE_STRAIN_MATERIAL_MATRIX,
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
        virtual calculate_matrix(MAST::ElemenetPropertyMatrixType t,
                                 DenseMatrix<Real>& m);
        
        
    protected:

        /*!
         *    pointer to the material property card
         */
        MAST::MaterialPropertyCardBase* _material;
        
    };
}



#endif // __MAST_element_property_matrix_base_h__
