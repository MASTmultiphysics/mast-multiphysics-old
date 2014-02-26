//
//  temperature.h
//  MAST
//
//  Created by Manav Bhatia on 1/7/14.
//  Copyright (c) 2014 Manav Bhatia. All rights reserved.
//

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
        void set_reference_temperature_function(MAST::FieldFunction<Number>& f) {
            _reference_temperature_function = &f;
        }
        
        /*!
         *   @returns reference to the function that provides the reference
         *   temperature
         */
        MAST::FieldFunction<Number>& reference_temperature_function() {
            libmesh_assert(_function);
            return *_reference_temperature_function;
        }
        
        /*!
         *   sets pointer to the function that provides the reference
         *   temperature
         */
        const MAST::FieldFunction<Number>& reference_temperature_function() const {
            libmesh_assert(_function);
            return *_reference_temperature_function;
        }
        
        
    protected:
        
        /*!
         *   reference temperature function
         */
        MAST::FieldFunction<Number>* _reference_temperature_function;
    };
}


#endif // __MAST_temperature_h__
