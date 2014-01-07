//
//  boundary_condition.h
//  MAST
//
//  Created by Manav Bhatia on 11/12/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_boundary_condition_h__
#define __MAST_boundary_condition_h__

// C++ includes
#include <memory>

// libMesh includes
#include "libmesh/function_base.h"

// MAST includes
#include "Base/MAST_data_types.h"



namespace MAST {
    
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
        
        
        void set_function(libMesh::FunctionBase<Number>& f) {
            _function = &f;
        }
        
        libMesh::FunctionBase<Number>& function() {
            libmesh_assert(_function);
            return *_function;
        }

        const libMesh::FunctionBase<Number>& function() const {
            libmesh_assert(_function);
            return *_function;
        }
        
    protected:
        
        MAST::BoundaryConditionType _bc_type;
        
        
        libMesh::FunctionBase<Number>* _function;
    };
    
    
    /*!
     *    builds a boundary condition object and returns it in a smart-pointer
     *    object
     */
    std::auto_ptr<MAST::BoundaryCondition>
    build_boundary_condition(MAST::BoundaryConditionType t);
}


#endif  // __MAST_boundary_condition_h__
