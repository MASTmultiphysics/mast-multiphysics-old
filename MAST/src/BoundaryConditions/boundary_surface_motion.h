//
//  surface_motion.h
//  MAST
//
//  Created by Manav Bhatia on 3/20/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

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
