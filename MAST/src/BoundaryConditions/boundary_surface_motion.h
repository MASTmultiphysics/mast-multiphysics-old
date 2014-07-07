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

#ifndef __MAST__surface_motion_base__
#define __MAST__surface_motion_base__

// libmesh includes
#include "libmesh/mesh.h"
#include "libmesh/point.h"
#include "libmesh/dense_vector.h"

// MAST includes
#include "BoundaryConditions/boundary_condition.h"





namespace MAST {
    
    class SurfaceMotionBase: public MAST::BoundaryCondition
    {
    public:
        SurfaceMotionBase():
        MAST::BoundaryCondition(MAST::SMALL_DISTURBANCE_DISPLACEMENT),
        frequency(0.),
        phase_offset(0.)
        { }
        
        virtual ~SurfaceMotionBase()
        { }
        
        virtual void zero()
        {
            frequency    = 0.;
            phase_offset = 0.;
        }
        
        /*!
         *   frequency of oscillation. If the reference chord, \p b_ref and
         *   reference velocity, \p V_ref are left to unity, then this is the
         *   dimensional frequency in rad/sec, otherwise this is the reduced
         *   frequency
         */
        Real frequency;
        
        
        /*!
         *    All transient motion data is based on a sine function that
         *    multiplies the motion amplitudes, so that the motion starts with
         *    zero amplitude at t=0. This constant can be set to pi/2
         *    to use a consine multiplier.
         */
        Real phase_offset;
        
        /*!
         *   calculation of surface velocity in frequency domain. \p u_trans is
         *   the pure translation velocity component, while \p dn_rot defines the
         *   surface normal perturbation
         */
        virtual void surface_velocity(const Real t,
                                      const libMesh::Point& p,
                                      const libMesh::Point& n,
                                      DenseComplexVector& w_trans,
                                      DenseComplexVector& u_trans,
                                      DenseComplexVector& dn_rot) = 0;
        
        /*!
         *   calculation of surface velocity in time domain. \p u_trans is
         *   the pure translation velocity component, while \p dn_rot defines the
         *   surface normal perturbation
         */
        virtual void surface_velocity(const Real t,
                                      const libMesh::Point& p,
                                      const libMesh::Point& n,
                                      DenseRealVector& w_trans,
                                      DenseRealVector& u_trans,
                                      DenseRealVector& dn_rot) = 0;
        
    protected:
        
        /*!
         *   initialization function for this object
         */
        virtual void init(Real freq, Real phase)
        {
            frequency = freq;
            phase_offset = phase;
        }
        
    };
}

#endif /* defined(__MAST__surface_motion__) */
