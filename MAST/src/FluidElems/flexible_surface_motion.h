//
//  flexible_surface_motion.h
//  MAST
//
//  Created by Manav Bhatia on 7/31/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef MAST_flexible_surface_motion_h
#define MAST_flexible_surface_motion_h

// MAST includes
#include "FluidElems/surface_motion.h"

// libMesh includes
#include "libmesh/mesh_function.h"
#include "libmesh/dof_map.h"


class FlexibleSurfaceMotion: public SurfaceMotionBase
{
public:
    FlexibleSurfaceMotion(System& sys):
    SurfaceMotionBase(),
    system(sys)
    { }
    
    virtual ~FlexibleSurfaceMotion()
    { }

    /*!
     *   system associated with the mesh and solution vector
     */
    System& system;

    virtual void zero();
    
    virtual void init(Real freq, Real phase,
                      NumericVector<Number>& sol);
    
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
    
protected:

    /*!
     *   mesh function that interpolates the solution
     */
    std::auto_ptr<MeshFunction> _function;

    /*!
     *    numeric vector that stores the solution
     */
    std::auto_ptr<NumericVector<Number> > _sol;
    
};



void
FlexibleSurfaceMotion::zero()
{
    libmesh_assert(false);
}


void
FlexibleSurfaceMotion::init(Real freq, Real phase,
                            NumericVector<Number>& sol)
{
    // first initialize the solution to the given vector
    if (!_sol.get())
    {
        _sol.reset(NumericVector<Number>::build(system.comm()).release());
        _sol->init(sol.size(), true, SERIAL);
    }
    
    // now localize the give solution to this objects's vector
    sol.localize(*_sol, system.get_dof_map().get_send_list());
    
    // if the mesh function has not been created so far, initialize it
    if (!_function.get())
    {
        std::vector<unsigned int> vars;
        system.get_all_variable_numbers(vars);
        _function.reset(new MeshFunction( system.get_equation_systems(),
                                         *_sol, system.get_dof_map(), vars));
        _function->init();
    }
    
    SurfaceMotionBase::init(freq, phase);
}



void
FlexibleSurfaceMotion::surface_velocity_frequency_domain(const Point& p,
                                                         const Point& n,
                                                         DenseVector<Complex>& u_trans,
                                                         DenseVector<Complex>& dn_rot)
{
    u_trans.zero();
    dn_rot.zero();
    
    libmesh_assert(_function.get()); // should be initialized before this call
    
    // translation is obtained by direct interpolation of the u,v,w vars
    DenseVector<Complex> v;
    (*_function)(p, 0., v);

    // now copy the values to u_trans
    Complex iota(0., 1.);
    for (unsigned int i=0; i<p.size(); i++)
        u_trans(i) = v(i) * iota * frequency;

    // perturbation of the normal requires calculation of the curl of
    // displacement at the given point
    std::vector<Gradient> gradients;
    _function->gradient(p, 0., gradients);

    // TODO: these need to be mapped from local 2D to 3D space
    
    // now prepare the rotation vector
    DenseVector<Complex> rot;
    rot.resize(3); 
    rot(0) = gradients[2](1) - gradients[1](2); // dwz/dy - dwy/dz
    rot(1) = gradients[0](2) - gradients[2](0); // dwx/dz - dwz/dx
    rot(2) = gradients[1](0) - gradients[0](1); // dwy/dx - dwx/dy
    
    // now do the cross-products
    dn_rot(0) =   rot(1) * n(2) - rot(2) * n(1);
    dn_rot(1) = -(rot(0) * n(2) - rot(2) * n(0));
    dn_rot(2) =   rot(0) * n(1) - rot(1) * n(0);
}



void
FlexibleSurfaceMotion::surface_velocity_time_domain(const Real t,
                                                 const Point& p,
                                                 const Point& n,
                                                 DenseVector<Real>& u_trans,
                                                 DenseVector<Real>& dn_rot)
{
    u_trans.zero();
    dn_rot.zero();
    
    libmesh_assert(_function.get()); // should be initialized before this call
    
    // translation is obtained by direct interpolation of the u,v,w vars
    DenseVector<Number> v;
    (*_function)(p, 0., v);
    
    // now copy the values to u_trans
    for (unsigned int i=0; i<p.size(); i++)
        u_trans(i) = std::real(v(i)) * frequency * cos(frequency*t +
                                                       phase_offset);
    
    // perturbation of the normal requires calculation of the curl of
    // displacement at the given point
    std::vector<Gradient> gradients;
    _function->gradient(p, 0., gradients);
    
    // TODO: these need to be mapped from local 2D to 3D space
    
    // now prepare the rotation vector
    DenseVector<Number> rot;
    rot(0) = gradients[2](1) - gradients[1](2); // dwz/dy - dwy/dz
    rot(1) = gradients[0](2) - gradients[2](0); // dwx/dz - dwz/dx
    rot(2) = gradients[1](0) - gradients[0](1); // dwy/dx - dwx/dy
    
    // now do the cross-products
    dn_rot(0) =   std::real(rot(1) * n(2) - rot(2) * n(1));
    dn_rot(1) =  -std::real((rot(0) * n(2) - rot(2) * n(0)));
    dn_rot(2) =   std::real(rot(0) * n(1) - rot(1) * n(0));
    
    dn_rot.scale(sin(frequency*t + phase_offset));
}



#endif
