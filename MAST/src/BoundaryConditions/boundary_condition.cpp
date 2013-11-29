//
//  boundary_condition.cpp
//  MAST
//
//  Created by Manav Bhatia on 11/29/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//


// MAST includes
#include "BoundaryConditions/boundary_condition.h"
#include "BoundaryConditions/surface_pressure.h"
#include "BoundaryConditions/small_disturbance_motion.h"


std::auto_ptr<MAST::BoundaryCondition>
MAST::build_boundary_condition(MAST::BoundaryConditionType t) {
    std::auto_ptr<MAST::BoundaryCondition> rval;
    
    switch (t) {
        case MAST::SURFACE_PRESSURE:
            rval.reset(new MAST::SurfacePressure);
            break;
            
        case MAST::SMALL_DISTURBANCE_MOTION:
            rval.reset(new MAST::SmallDisturbanceMotion);
            break;

        default:
            libmesh_error(); // should not get here
            break;
    }
    
    return rval;
}

