//
//  surface_motion.h
//  FESystem
//
//  Created by Manav Bhatia on 3/20/13.
//
//

#ifndef __FESystem__surface_motion_base__
#define __FESystem__surface_motion_base__

// libmesh includes
#include "libmesh/mesh.h"
#include "libmesh/point.h"
#include "libmesh/dense_vector.h"


using namespace libMesh;


class SurfaceMotionBase
{
public:
    SurfaceMotionBase():
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
    virtual void surface_velocity_frequency_domain(const Point& p,
                                                   const Point& n,
                                                   DenseVector<Complex>& u_trans,
                                                   DenseVector<Complex>& dn_rot) = 0;

    /*!
     *   calculation of surface velocity in time domain. \p u_trans is
     *   the pure translation velocity component, while \p dn_rot defines the
     *   surface normal perturbation
     */
    virtual void surface_velocity_time_domain(const Real t,
                                              const Point& p,
                                              const Point& n,
                                              DenseVector<Number>& u_trans,
                                              DenseVector<Number>& dn_rot) = 0;

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


#endif /* defined(__FESystem__surface_motion__) */
