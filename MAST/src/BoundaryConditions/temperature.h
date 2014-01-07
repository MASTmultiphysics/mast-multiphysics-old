//
//  temperature.h
//  MAST
//
//  Created by Manav Bhatia on 1/7/14.
//  Copyright (c) 2014 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_temperature_h__
#define __MAST_temperature_h__

// libMesh includes
#include "libmesh/mesh_function.h"
#include "libmesh/system.h"
#include "libmesh/dof_map.h"
#include "libmesh/equation_systems.h"
#include "libmesh/mesh_serializer.h"
#include "libmesh/numeric_vector.h"


// MAST includes
#include "BoundaryConditions/boundary_condition.h"
#include "Flight/flight_condition.h"
#include "FluidElems/frequency_domain_linearized_fluid_system.h"
#include "FluidElems/fluid_elem_base.h"


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
        void set_reference_temperature_function(libMesh::FunctionBase<Number>& f) {
            _reference_temperature_function = &f;
        }
        
        /*!
         *   @returns reference to the function that provides the reference
         *   temperature
         */
        libMesh::FunctionBase<Number>& reference_temperature_function() {
            libmesh_assert(_function);
            return *_reference_temperature_function;
        }
        
        /*!
         *   sets pointer to the function that provides the reference
         *   temperature
         */
        const libMesh::FunctionBase<Number>& reference_temperature_function() const {
            libmesh_assert(_function);
            return *_reference_temperature_function;
        }
        
        
    protected:
        
        /*!
         *   reference temperature function
         */
        libMesh::FunctionBase<Number>* _reference_temperature_function;
    };
}


#endif // __MAST_temperature_h__
