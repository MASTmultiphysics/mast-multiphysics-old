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
#include "BoundaryConditions/boundary_surface_motion.h"
#include "Numerics/utility.h"

// libMesh includes
#include "libmesh/mesh_function.h"
#include "libmesh/dof_map.h"
#include "libmesh/mesh_serializer.h"


namespace MAST {
    
    class FlexibleSurfaceMotion: public MAST::SurfaceMotionBase
    {
    public:
        FlexibleSurfaceMotion(libMesh::System& sys):
        SurfaceMotionBase(),
        system(sys)
        {
            MeshBase& mesh = sys.get_mesh();
            _mesh_serializer.reset(new MeshSerializer(mesh, true));
        }
        
        virtual ~FlexibleSurfaceMotion()
        { }
        
        /*!
         *   system associated with the mesh and solution vector
         */
        libMesh::System& system;
        
        virtual void zero();
        
        virtual void init(libMesh::Real freq, libMesh::Real phase,
                          libMesh::NumericVector<libMesh::Real>& sol);
        
        /*!
         *   calculation of surface velocity in frequency domain. \p u_trans is
         *   the pure translation velocity component, while \p dn_rot defines the
         *   surface normal perturbation
         */
        virtual void surface_velocity(const libMesh::Real t,
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
        virtual void surface_velocity(const libMesh::Real t,
                                      const libMesh::Point& p,
                                      const libMesh::Point& n,
                                      DenseRealVector& w_trans,
                                      DenseRealVector& u_trans,
                                      DenseRealVector& dn_rot);
        
    protected:
        
        /*!
         *   mesh function that interpolates the solution
         */
        std::auto_ptr<MeshFunction> _function;
        
        /*!
         *    numeric vector that stores the solution
         */
        std::auto_ptr<libMesh::NumericVector<libMesh::Real> > _sol;
        
        
        /*!
         *   this serializes the mesh for use in interpolation
         */
        std::auto_ptr<MeshSerializer> _mesh_serializer;
    };
}


inline
void
MAST::FlexibleSurfaceMotion::zero()
{
    libmesh_assert(false);
}


inline
void
MAST::FlexibleSurfaceMotion::init(libMesh::Real freq, libMesh::Real phase,
                                  libMesh::NumericVector<libMesh::Real>& sol)
{
    // first initialize the solution to the given vector
    if (!_sol.get())
    {
        _sol.reset(libMesh::NumericVector<libMesh::Real>::build(system.comm()).release());
        _sol->init(sol.size(), true, SERIAL);
    }
    
    // now localize the give solution to this objects's vector
    sol.localize(*_sol);
    
    // if the mesh function has not been created so far, initialize it
    if (!_function.get())
    {
        std::vector<unsigned int> vars(3);
        vars[0] = system.variable_number("ux");
        vars[1] = system.variable_number("uy");
        vars[2] = system.variable_number("uz");
        _function.reset(new MeshFunction( system.get_equation_systems(),
                                         *_sol, system.get_dof_map(), vars));
        _function->init();
    }
    
    MAST::SurfaceMotionBase::init(freq, phase);
}



inline
void
MAST::FlexibleSurfaceMotion::surface_velocity(const libMesh::Real t,
                                              const libMesh::Point& p,
                                              const libMesh::Point& n,
                                              DenseComplexVector& w_trans,
                                              DenseComplexVector& u_trans,
                                              DenseComplexVector& dn_rot)
{
//#ifdef LIBMESH_USE_COMPLEX_NUMBERS
    w_trans.zero();
    u_trans.zero();
    dn_rot.zero();
    
    libmesh_assert(_function.get()); // should be initialized before this call
    
    // translation is obtained by direct interpolation of the u,v,w vars
    DenseRealVector v_real;
    DenseComplexVector v;
    (*_function)(p, t, v_real);
    v = v_real;
    
    // now copy the values to u_trans
    libMesh::Complex iota(0., 1.);
    for (unsigned int i=0; i<3; i++) {
        w_trans(i) = v(i);
        u_trans(i) = v(i) * iota * frequency;
    }
    
    // perturbation of the normal requires calculation of the curl of
    // displacement at the given point
    std::vector<Gradient> gradients;
    _function->gradient(p, 0., gradients);
    
    // TODO: these need to be mapped from local 2D to 3D space
    
    // now prepare the rotation vector
    DenseComplexVector rot;
    rot.resize(3);
    rot(0) = gradients[2](1) - gradients[1](2); // dwz/dy - dwy/dz
    rot(1) = gradients[0](2) - gradients[2](0); // dwx/dz - dwz/dx
    rot(2) = gradients[1](0) - gradients[0](1); // dwy/dx - dwx/dy
    
    // now do the cross-products
    dn_rot(0) =   rot(1) * n(2) - rot(2) * n(1);
    dn_rot(1) = -(rot(0) * n(2) - rot(2) * n(0));
    dn_rot(2) =   rot(0) * n(1) - rot(1) * n(0);

}



inline
void
MAST::FlexibleSurfaceMotion::surface_velocity(const libMesh::Real t,
                                              const libMesh::Point& p,
                                              const libMesh::Point& n,
                                              DenseRealVector& w_trans,
                                              DenseRealVector& u_trans,
                                              DenseRealVector& dn_rot)
{
    u_trans.zero();
    dn_rot.zero();
    
    libmesh_assert(_function.get()); // should be initialized before this call
    
    // translation is obtained by direct interpolation of the u,v,w vars
    DenseRealVector v;
    (*_function)(p, t, v);
    
    // now copy the values to u_trans
    for (unsigned int i=0; i<3; i++) {
        w_trans(i) = v(i);
        u_trans(i) = v(i) * frequency * cos(frequency*t + phase_offset);
    }
    
    // perturbation of the normal requires calculation of the curl of
    // displacement at the given point
    std::vector<Gradient> gradients;
    _function->gradient(p, 0., gradients);
    
    // TODO: these need to be mapped from local 2D to 3D space
    
    // now prepare the rotation vector
    DenseRealVector rot;
    rot.resize(3);
    rot(0) = gradients[2](1) - gradients[1](2); // dwz/dy - dwy/dz
    rot(1) = gradients[0](2) - gradients[2](0); // dwx/dz - dwz/dx
    rot(2) = gradients[1](0) - gradients[0](1); // dwy/dx - dwx/dy
    
    // now do the cross-products
    dn_rot(0) =    rot(1) * n(2) - rot(2) * n(1);
    dn_rot(1) =  -(rot(0) * n(2) - rot(2) * n(0));
    dn_rot(2) =    rot(0) * n(1) - rot(1) * n(0);
    
    dn_rot.scale(sin(frequency*t + phase_offset));
}



#endif
