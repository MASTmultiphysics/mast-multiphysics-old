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

