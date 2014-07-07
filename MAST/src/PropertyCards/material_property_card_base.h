/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2014  Manav Bhatia
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

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

        virtual std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > >
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
