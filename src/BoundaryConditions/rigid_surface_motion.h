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

#ifndef MAST_rigid_surface_motion_h
#define MAST_rigid_surface_motion_h

// MAST includes
#include "BoundaryConditions/boundary_surface_motion.h"

namespace MAST {
    
    
    class RigidSurfaceMotion: public MAST::SurfaceMotionBase
    {
    public:
        RigidSurfaceMotion();
        
        virtual ~RigidSurfaceMotion();
        
        libMesh::Point plunge_vector;
        
        libMesh::Point pitch_axis;
        
        libMesh::Point hinge_location;
        
        Real plunge_amplitude;
        
        Real pitch_amplitude;
        
        Real pitch_phase;
        
        virtual void zero();
        
        virtual void init(const Real freq, const Real vel, const Real phase);
        
        /*!
         *   calculation of surface velocity in frequency domain. 
         *   \p u_trans is the pure translation velocity component,
         *   while \p dn_rot defines the surface normal perturbation.
         */
        virtual void surface_velocity(const Real t,
                                      const libMesh::Point& p,
                                      const libMesh::Point& n,
                                      DenseComplexVector& w_trans,
                                      DenseComplexVector& u_trans,
                                      DenseComplexVector& dn_rot);
        
        /*!
         *   calculation of sensitivity of surface velocity components wrt the
         *   frequency. \p u_trans is
         *   the pure translation velocity component, while \p dn_rot defines the
         *   surface normal perturbation
         */
        virtual void surface_velocity_k_sens(const Real t,
                                             const libMesh::Point& p,
                                             const libMesh::Point& n,
                                             DenseComplexVector& w_trans,
                                             DenseComplexVector& u_trans,
                                             DenseComplexVector& dn_rot) {
            // to be implemented
            libmesh_assert(false);
        }
        
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
                                      DenseRealVector& dn_rot);
    };
}

#endif
