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
            libMesh::MeshBase& mesh = sys.get_mesh();
            _mesh_serializer.reset(new libMesh::MeshSerializer(mesh, true));
        }
        
        virtual ~FlexibleSurfaceMotion()
        { }
        
        /*!
         *   system associated with the mesh and solution vector
         */
        libMesh::System& system;
        
        virtual void zero();
        
        virtual void init(Real freq, Real phase,
                          libMesh::NumericVector<Real>& sol);
        
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
        
    protected:
        
        /*!
         *   mesh function that interpolates the solution
         */
        std::auto_ptr<libMesh::MeshFunction> _function;
        
        /*!
         *    numeric vector that stores the solution
         */
        std::auto_ptr<libMesh::NumericVector<Real> > _sol;
        
        
        /*!
         *   this serializes the mesh for use in interpolation
         */
        std::auto_ptr<libMesh::MeshSerializer> _mesh_serializer;
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
MAST::FlexibleSurfaceMotion::init(Real freq, Real phase,
                                  libMesh::NumericVector<Real>& sol)
{
    // first initialize the solution to the given vector
    if (!_sol.get())
    {
        _sol.reset(libMesh::NumericVector<Real>::build(system.comm()).release());
        _sol->init(sol.size(), true, libMesh::SERIAL);
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
        _function.reset(new libMesh::MeshFunction( system.get_equation_systems(),
                                         *_sol, system.get_dof_map(), vars));
        _function->init();
    }
    
    MAST::SurfaceMotionBase::init(freq, phase);
    
    //sol.print();
}



inline
void
MAST::FlexibleSurfaceMotion::surface_velocity(const Real t,
                                              const libMesh::Point& p,
                                              const libMesh::Point& n,
                                              DenseComplexVector& w_trans,
                                              DenseComplexVector& u_trans,
                                              DenseComplexVector& dn_rot)
{
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
    Complex iota(0., 1.);
    for (unsigned int i=0; i<3; i++) {
        w_trans(i) = v(i);
        u_trans(i) = v(i) * iota * frequency;
    }
    
    // perturbation of the normal requires calculation of the curl of
    // displacement at the given point
    std::vector<libMesh::Gradient> gradients;
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
MAST::FlexibleSurfaceMotion::surface_velocity_k_sens(const Real t,
                                                     const libMesh::Point& p,
                                                     const libMesh::Point& n,
                                                     DenseComplexVector& w_trans,
                                                     DenseComplexVector& u_trans,
                                                     DenseComplexVector& dn_rot)
{
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
    Complex iota(0., 1.);
    for (unsigned int i=0; i<3; i++) {
        w_trans(i) = v(i);
        u_trans(i) = v(i) * iota;
    }
    
    // perturbation of the normal requires calculation of the curl of
    // displacement at the given point
    std::vector<libMesh::Gradient> gradients;
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

    // zero the components that are independent of the frequency
    w_trans.zero();
    dn_rot.zero();
}



inline
void
MAST::FlexibleSurfaceMotion::surface_velocity(const Real t,
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
    std::vector<libMesh::Gradient> gradients;
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
    
    //dn_rot.scale(sin(frequency*t + phase_offset));
}



#endif
