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
         *   libMesh::Pointer to the small disturbance deformation for this boundary
         *   condition
         */
        MAST::SurfaceMotionBase* _deformation;
        
        
        /*!
         *   libMesh::Pointer to the small disturbance pressure for this boundary condition
         */
        MAST::SmallDisturbanceSurfacePressure* _pressure;
        
    };
}


#endif  //__MAST_small_disturbance_motion_h__

