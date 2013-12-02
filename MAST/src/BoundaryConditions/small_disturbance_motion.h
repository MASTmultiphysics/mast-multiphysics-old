//
//  small_disturbance_motion.h
//  MAST
//
//  Created by Manav Bhatia on 11/29/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_small_disturbance_motion_h__
#define __MAST_small_disturbance_motion_h__

// MAST includes
#include "BoundaryConditions/boundary_condition.h"
#include "BoundaryConditions/boundary_surface_motion.h"
#include "BoundaryConditions/surface_pressure.h"


namespace MAST {
    class SmallDisturbanceMotion : public MAST::BoundaryCondition {
        
    public:
        SmallDisturbanceMotion():
        MAST::BoundaryCondition(MAST::SMALL_DISTURBANCE_MOTION),
        _deformation(NULL),
        _pressure(NULL)
        { }
        
        
        virtual ~SmallDisturbanceMotion()
        { }

        
        /*!
         *   sets the deformation object
         */
        void set_deformation(MAST::SurfaceMotionBase& deform) {
            libmesh_assert(!_deformation);
            _deformation = &deform;
        }

        
        /*!
         *   sets the deformation object
         */
        void set_pressure(MAST::SmallDisturbanceSurfacePressure& press) {
            libmesh_assert(!_pressure);
            _pressure = &press;
        }

        
        /*!
         *   @returns a reference to the deformation object
         */
        MAST::SurfaceMotionBase& get_deformation() {
            libmesh_assert(_deformation);
            return *_deformation;
        }
        
        
        /*!
         *   @returns a reference to the surface pressure object
         */
        MAST::SmallDisturbanceSurfacePressure& get_pressure() {
            libmesh_assert(_pressure);
            return *_pressure;
        }

        
    protected:
        
        /*!
         *   Pointer to the small disturbance deformation for this boundary
         *   condition
         */
        MAST::SurfaceMotionBase* _deformation;
        
        
        /*!
         *   Pointer to the small disturbance pressure for this boundary condition
         */
        MAST::SmallDisturbanceSurfacePressure* _pressure;
        
    };
}


#endif  //__MAST_small_disturbance_motion_h__

