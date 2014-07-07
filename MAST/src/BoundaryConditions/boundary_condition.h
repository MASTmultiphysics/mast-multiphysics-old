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

#ifndef __MAST_boundary_condition_h__
#define __MAST_boundary_condition_h__

// C++ includes
#include <memory>

// libMesh includes
#include "libmesh/function_base.h"

// MAST includes
#include "Base/MAST_data_types.h"



namespace MAST {
    
    // Forward declerations
    class FieldFunctionBase;
    
    enum BoundaryConditionType {
        SURFACE_PRESSURE,
        SMALL_DISTURBANCE_PRESSURE, // provides pressure perturbations about steady values
        SMALL_DISTURBANCE_MOTION, // provides pressure and motion
                                  //perturbations about steady state values
        PISTON_THEORY,
        DISPLACEMENT_DIRICHLET,
        SMALL_DISTURBANCE_DISPLACEMENT,
        TEMPERATURE
    };
    
    
    class BoundaryCondition {
        
    public:
        BoundaryCondition(MAST::BoundaryConditionType t):
        _bc_type(t),
        _function(NULL)
        {}
        
        virtual ~BoundaryCondition() { }
        
        
        MAST::BoundaryConditionType type() const {
            return _bc_type;
        }
        
        
        void set_function(MAST::FieldFunctionBase& f) {
            _function = &f;
        }
        
        MAST::FieldFunctionBase& function() {
            libmesh_assert(_function);
            return *_function;
        }

        const MAST::FieldFunctionBase& function() const {
            libmesh_assert(_function);
            return *_function;
        }
        
    protected:
        
        MAST::BoundaryConditionType _bc_type;
        
        
        MAST::FieldFunctionBase* _function;
    };
    
    
    /*!
     *    builds a boundary condition object and returns it in a smart-pointer
     *    object
     */
    std::auto_ptr<MAST::BoundaryCondition>
    build_boundary_condition(MAST::BoundaryConditionType t);
}


#endif  // __MAST_boundary_condition_h__
