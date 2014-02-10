//
//  material_property_card_base.h
//  MAST
//
//  Created by Manav Bhatia on 10/15/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_material_property_card_base_h__
#define __MAST_material_property_card_base_h__

// MAST includes
#include "PropertyCards/property_card_base.h"

// libMesh includes
#include "libmesh/dense_matrix.h"

namespace MAST
{
    // Forward decleration
    template <typename ValType> class FieldFunction;
    class ElementPropertyCardBase;
    
    enum MaterialPropertyMatrixType {
        MATERIAL_STIFFNESS_MATRIX,
        MATERIAL_DAMPING_MATRIX,
        MATERIAL_INERTIA_MATRIX,
        MATERIAL_THERMAL_EXPANSION_MATRIX,
        MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX
    };

    
    class  MaterialPropertyCardBase: public MAST::PropertyCardBase {
        
    public:
        
        MaterialPropertyCardBase(unsigned int pid):
        MAST::PropertyCardBase (),
        _pid(pid)
        { }
        
        /*!
         *    returns the id for this card
         */
        unsigned int id() const {
            return _pid;
        }

        virtual std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>>
        get_property(MAST::MaterialPropertyMatrixType t,
                     const MAST::ElementPropertyCardBase& p,
                     const unsigned int dim) const = 0;
        
    protected:
        
        /*!
         *   property card id
         */
        unsigned int _pid;
        
    };
    
    
}


#endif // __MAST_material_property_card_base_h__
