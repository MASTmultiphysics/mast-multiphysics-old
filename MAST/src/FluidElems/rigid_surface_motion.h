//
//  rigid_surface_motion.h
//  MAST
//
//  Created by Manav Bhatia on 7/31/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef MAST_rigid_surface_motion_h
#define MAST_rigid_surface_motion_h

// MAST includes
#include "FluidElems/surface_motion.h"

class RigidSurfaceMotion: public SurfaceMotionBase
{
public:
    RigidSurfaceMotion();
    
    virtual ~RigidSurfaceMotion();
    
    Point plunge_vector;
    
    Point pitch_axis;
    
    Point hinge_location;
    
    Real plunge_amplitude;
    
    Real pitch_amplitude;
    
    Real pitch_phase;
    
    virtual void zero();
    
    virtual void init(Real freq, Real phase);
    
    /*!
     *   calculation of surface velocity in frequency domain. \p u_trans is
     *   the pure translation velocity component, while \p dn_rot defines the
     *   surface normal perturbation
     */
    virtual void surface_velocity_frequency_domain(const Point& p,
                                                   const Point& n,
                                                   DenseVector<Complex>& u_trans,
                                                   DenseVector<Complex>& dn_rot);
    
    /*!
     *   calculation of surface velocity in time domain. \p u_trans is
     *   the pure translation velocity component, while \p dn_rot defines the
     *   surface normal perturbation
     */
    virtual void surface_velocity_time_domain(const Real t,
                                              const Point& p,
                                              const Point& n,
                                              DenseVector<Real>& u_trans,
                                              DenseVector<Real>& dn_rot);
};

#endif
