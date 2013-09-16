//
//  surface_motion_function.h
//  MAST
//
//  Created by Manav Bhatia on 9/11/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef MAST_surface_motion_function_h
#define MAST_surface_motion_function_h

// C++ includes
#include <functional>

// MAST includes
#include "FluidElems/surface_motion.h"

class SurfaceMotionFunction: public SurfaceMotionBase
{
public:
    SurfaceMotionFunction();
    
    virtual ~SurfaceMotionFunction();
        
    /*!
     *   The functors that provide the velocity and displacement derivative with
     *   respect to each spatial coordinate.
     */
//    virtual void init(std::function<Point(const Point &p, const Real t)> wdot,
//                      std::function<Point(const Point &p, const Real t)> dwdx,
//                      std::function<Point(const Point &p, const Real t)> dwdy,
//                      std::function<Point(const Point &p, const Real t)> dwdz){
//        _wdot = wdot;
//        _dwdx = dwdx;
//        _dwdy = dwdy;
//        _dwdz = dwdz;
//    }
    
    
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
                                              DenseVector<Number>& u_trans,
                                              DenseVector<Number>& dn_rot);
    
protected:
    
//    /*!
//     *   surface velocity functor
//     */
//    std::function<Point(const Point &p, const Real t)> _wdot;
//
//    /*!
//     *   surface displacement gradient with respect to x-coordinate
//     */
//    std::function<Point(const Point &p, const Real t)> _dwdx;
//
//    /*!
//     *   surface displacement gradient with respect to y-coordinate
//     */
//    std::function<Point(const Point &p, const Real t)> _dwdy;
//
//    /*!
//     *   surface displacement gradient with respect to z-coordinate
//     */
//    std::function<Point(const Point &p, const Real t)> _dwdz;

};


SurfaceMotionFunction::SurfaceMotionFunction():
SurfaceMotionBase()
{
    
}


SurfaceMotionFunction::~SurfaceMotionFunction()
{
    
}





void
SurfaceMotionFunction::surface_velocity_frequency_domain(const Point& p,
                                                      const Point& n,
                                                      DenseVector<Complex>& u_trans,
                                                      DenseVector<Complex>& dn_rot)
{
    libmesh_assert(false); // to be implemented
}



void
SurfaceMotionFunction::surface_velocity_time_domain(const Real t,
                                                 const Point& p,
                                                 const Point& n,
                                                 DenseVector<Number>& u_trans,
                                                 DenseVector<Number>& dn_rot)
{
    u_trans.zero();
    dn_rot.zero();
    
    // translation is obtained by direct interpolation of the u,v,w vars
    Point v, dwdx, dwdy, dwdz;
    const Real x = p(0) - 0.2, tc = 0.05, c = 0.2, h = tc*c/2., pi = acos(-1.);
    if ( (x>=0.) && (x <= 0.2)){
        v(1)     = h*sin(pi*x/c);
        dwdx(1)  = h*pi/c*cos(pi*x/c);
    }
    
    // now copy the values to u_trans
    for (unsigned int i=0; i<3; i++)
        u_trans(i) = 0.*v(i);
    
    // now prepare the rotation vector
    DenseVector<Number> rot;
    rot.resize(3);
    rot(0) = dwdy(2) - dwdz(1); // dwz/dy - dwy/dz
    rot(1) = dwdz(0) - dwdx(2); // dwx/dz - dwz/dx
    rot(2) = dwdx(1) - dwdy(0); // dwy/dx - dwx/dy
    
    // now do the cross-products
    dn_rot(0) =    rot(1) * n(2) - rot(2) * n(1);
    dn_rot(1) =  -(rot(0) * n(2) - rot(2) * n(0));
    dn_rot(2) =    rot(0) * n(1) - rot(1) * n(0);
}



#endif
