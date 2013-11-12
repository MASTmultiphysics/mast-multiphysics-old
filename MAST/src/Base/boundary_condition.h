//
//  boundary_condition.h
//  MAST
//
//  Created by Manav Bhatia on 11/12/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_boundary_condition_h__
#define __MAST_boundary_condition_h__


// MAST includes
#include "Numerics/function_base.h"


namespace MAST {
    
    enum BoundaryConditionType {
        SURFACE_PRESSURE,
        DISPLACEMENT,
        INVALID_BOUNDARY_CONDITION
    };
    
    
    class BoundaryCondition {
        
    public:
        BoundaryCondition():
        _bc_type(MAST::INVALID_BOUNDARY_CONDITION),
        _function(NULL)
        {}
        
        virtual ~BoundaryCondition() { }
        
        void set_type(MAST::BoundaryConditionType t) {
            _bc_type = t;
        }
        
        MAST::BoundaryConditionType type() const {
            return _bc_type;
        }
        
        
        void set_function(MAST::FunctionBase& f) {
            _function = &f;
        }
        
        const MAST::FunctionBase& function() const {
            libmesh_assert(_function);
            return _function;
        }
        
    protected:
        
        MAST::BoundaryConditionType _bc_type;
        
        
        MAST::FunctionBase* _function;
    };
}


#endif  // __MAST_boundary_condition_h__
