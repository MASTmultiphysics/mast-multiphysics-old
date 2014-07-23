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

#ifndef MAST_surface_motion_function_h
#define MAST_surface_motion_function_h

// C++ includes
#include <functional>

// MAST includes
#include "BoundaryConditions/boundary_surface_motion.h"

namespace MAST {
    
    class SurfaceMotionFunction: public MAST::SurfaceMotionBase
    {
    public:
        SurfaceMotionFunction();
        
        virtual ~SurfaceMotionFunction();
        
        /*!
         *   The functors that provide the velocity and displacement derivative with
         *   respect to each spatial coordinate.
         */
        //    virtual void init(std::function<Point(const libMesh::Point &p, const Real t)> wdot,
        //                      std::function<Point(const libMesh::Point &p, const Real t)> dwdx,
        //                      std::function<Point(const libMesh::Point &p, const Real t)> dwdy,
        //                      std::function<Point(const libMesh::Point &p, const Real t)> dwdz){
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
        
    protected:
        
        //    /*!
        //     *   surface velocity functor
        //     */
        //    std::function<Point(const libMesh::Point &p, const Real t)> _wdot;
        //
        //    /*!
        //     *   surface displacement gradient with respect to x-coordinate
        //     */
        //    std::function<Point(const libMesh::Point &p, const Real t)> _dwdx;
        //
        //    /*!
        //     *   surface displacement gradient with respect to y-coordinate
        //     */
        //    std::function<Point(const libMesh::Point &p, const Real t)> _dwdy;
        //
        //    /*!
        //     *   surface displacement gradient with respect to z-coordinate
        //     */
        //    std::function<Point(const libMesh::Point &p, const Real t)> _dwdz;
        
    };
}


inline
MAST::SurfaceMotionFunction::SurfaceMotionFunction():
MAST::SurfaceMotionBase()
{
    
}


inline
MAST::SurfaceMotionFunction::~SurfaceMotionFunction()
{
    
}





inline void
MAST::SurfaceMotionFunction::surface_velocity(const Real t,
                                              const libMesh::Point& p,
                                              const libMesh::Point& n,
                                              DenseComplexVector& w_trans,
                                              DenseComplexVector& u_trans,
                                              DenseComplexVector& dn_rot)
{
    // use only on the upper surface
    if (n(2) < 0. &&
        (p(0) >= 14.0e3 && p(0) <= 14.2e3) &&
        (p(1) >= -2.0e3 && p(1) <= -1.8e3)) {
        
        u_trans.zero();
        dn_rot.zero();
        
        // translation is obtained by direct interpolation of the u,v,w vars
        libMesh::Point v, dwdx, dwdy, dwdz;
        const Real x = p(0) - 2., tc = 0.001, c = 2., h = tc*c/2., pi = acos(-1.),
        omega=124.874;
        
        const Complex iota(0., 1.);
        
        if ( (x>=0.) && (x <= c)) {
            v(1)     = h*(1.-cos(pi*2.*x/c));
            dwdx(1)  = h*2.*pi/c*sin(pi*2.*x/c);
        }
        
        // now copy the values to u_trans
        for (unsigned int i=0; i<3; i++) {
            w_trans(i) = v(i);
            u_trans(i) = iota*omega*v(i);
        }
        
        // now prepare the rotation vector
        DenseRealVector rot;
        rot.resize(3);
        rot(0) = dwdy(2) - dwdz(1); // dwz/dy - dwy/dz
        rot(1) = dwdz(0) - dwdx(2); // dwx/dz - dwz/dx
        rot(2) = dwdx(1) - dwdy(0); // dwy/dx - dwx/dy
        
        // now do the cross-products
        dn_rot(0) =    rot(1) * n(2) - rot(2) * n(1);
        dn_rot(1) =  -(rot(0) * n(2) - rot(2) * n(0));
        dn_rot(2) =    rot(0) * n(1) - rot(1) * n(0);
    }
    else {
        w_trans.resize(3);
        u_trans.resize(3);
        dn_rot.resize(3);
    }
}





inline void
MAST::SurfaceMotionFunction::surface_velocity(const Real t,
                                              const libMesh::Point& p,
                                              const libMesh::Point& n,
                                              DenseRealVector& w_trans,
                                              DenseRealVector& u_trans,
                                              DenseRealVector& dn_rot)
{
    u_trans.zero();
    dn_rot.zero();
    
    // translation is obtained by direct interpolation of the u,v,w vars
    libMesh::Point v, dwdx, dwdy, dwdz;
    const Real x = p(0) - 2., tc = 0.001, c = 2., h = tc*c/2., pi = acos(-1.),
    omega=124.874;
    
    if ( (x>=0.) && (x <= c)){
        v(1)     = h*(1.-cos(pi*2.*x/c));
        dwdx(1)  = h*2.*pi/c*sin(pi*2.*x/c)*sin(omega*t);
    }
    
    // now copy the values to u_trans
    for (unsigned int i=0; i<3; i++) {
        w_trans(i) = sin(omega*t)*v(i);
        u_trans(i) = omega*cos(omega*t)*v(i);
    }
    
    // now prepare the rotation vector
    DenseRealVector rot;
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
