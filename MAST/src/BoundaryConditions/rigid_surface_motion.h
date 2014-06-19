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
        
        virtual void init(Real freq, Real phase);
        
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
