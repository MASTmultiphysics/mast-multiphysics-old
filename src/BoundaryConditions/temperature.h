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

#ifndef __MAST_temperature_h__
#define __MAST_temperature_h__


// MAST includes
#include "BoundaryConditions/boundary_condition.h"


namespace MAST {
    
    class Temperature: public MAST::BoundaryCondition {
    public:
        Temperature():
        MAST::BoundaryCondition(MAST::TEMPERATURE)
        { }
        
        virtual ~Temperature()
        { }
        
        /*!
         *   sets pointer to the function that provides the reference
         *   temperature
         */
        void set_reference_temperature_function(MAST::FieldFunction<Real>& f) {
            _reference_temperature_function = &f;
        }
        
        /*!
         *   @returns reference to the function that provides the reference
         *   temperature
         */
        MAST::FieldFunction<Real>& reference_temperature_function() {
            libmesh_assert(_function);
            return *_reference_temperature_function;
        }
        
        /*!
         *   sets pointer to the function that provides the reference
         *   temperature
         */
        const MAST::FieldFunction<Real>& reference_temperature_function() const {
            libmesh_assert(_function);
            return *_reference_temperature_function;
        }
        
        
    protected:
        
        /*!
         *   reference temperature function
         */
        MAST::FieldFunction<Real>* _reference_temperature_function;
    };
}


#endif // __MAST_temperature_h__
