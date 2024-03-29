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

#ifndef __MAST_gaussian_bump_mesh_h__
#define __MAST_gaussian_bump_mesh_h__


// MAST includes
#include "Mesh/mesh_initializer.h"
#include "BoundaryConditions/boundary_surface_motion.h"


class GaussianBumpMesh2D: public MeshInitializer
{
public:
    GaussianBumpMesh2D():
    _height(0.),
    _x0(0.), _x1(0.),
    _y0(0.), _y1(0.)
    { }
    
    /*!
     *   initializes the object with the division for each dimension.
     */
    virtual void init( const Real height,
                       const std::vector<MeshInitializer::CoordinateDivisions*>& divs,
                       libMesh::UnstructuredMesh& mesh, libMesh::ElemType t);
    
    Real x0() const {return _x0;}

    Real x1() const {return _x1;}

    Real y0() const {return _y0;}

    Real y1() const {return _y1;}
    
    Real h() const {return _height;}

protected:
    
    /*!
     *   move the mesh points to account for a bump, and apply boudnary condition
     */
    virtual void process_mesh( );

    /*!
     *   Height of the bump
     */
    Real _height;
    
    /*!
     *   Extent of computational domain
     */
    Real _x0, _x1;
    Real _y0, _y1;
};



class GaussianBumpMesh3D: public MeshInitializer
{
public:
    GaussianBumpMesh3D():
    _height(0.),
    _x0(0.), _x1(0.),
    _y0(0.), _y1(0.),
    _z0(0.), _z1(0.)
    { }
    
    /*!
     *   initializes the object with the division for each dimension.
     */
    virtual void init(const Real height,
                      const std::vector<MeshInitializer::CoordinateDivisions*>& divs,
                      libMesh::UnstructuredMesh& mesh, libMesh::ElemType t);
    
protected:
    
    /*!
     *   move the mesh points to account for a bump, and apply boudnary condition
     */
    virtual void process_mesh( );
    
    /*!
     *   Height of the bump
     */
    Real _height;
    
    /*!
     *   Extent of computational domain
     */
    Real _x0, _x1;
    Real _y0, _y1;
    Real _z0, _z1;
};




class GaussianBumpSurfaceNormalCorrection2D: public MAST::SurfaceMotionBase
{
public:
    GaussianBumpSurfaceNormalCorrection2D(Real x0, Real x1, Real h):
    MAST::SurfaceMotionBase(),
    _x0(x0),
    _x1(x1),
    _h(h)
    { }
    
    virtual ~GaussianBumpSurfaceNormalCorrection2D() {}
    
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
                                  DenseComplexVector& dn_rot)
    { libmesh_error();}
    
    virtual void surface_velocity_k_sens(const Real t,
                                         const libMesh::Point& p,
                                         const libMesh::Point& n,
                                         DenseComplexVector& w_trans,
                                         DenseComplexVector& u_trans,
                                         DenseComplexVector& dn_rot) {
        libmesh_error();
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
    
    Real _x0, _x1, _h;
};



inline void
GaussianBumpSurfaceNormalCorrection2D::surface_velocity(const Real t,
                                                        const libMesh::Point& p,
                                                        const libMesh::Point& n,
                                                        DenseRealVector& w_trans,
                                                        DenseRealVector& u_trans,
                                                        DenseRealVector& dn_rot)
{
    w_trans.zero();
    u_trans.zero();
    dn_rot.zero();
    
    // for the point p, add the correction of surface normal n to dn_rot
    if (p(1) <= _h) {
        Real x = p(0),
        f = -25. * pow((x-_x0)-(_x1-_x0)*0.5, 2.0),
        dfdx = -50.*((x-_x0)-(_x1-_x0)*0.5),
        dydx = _h*exp(f)*dfdx;
        
        dn_rot(0) = dydx; dn_rot(1) = -1.;
        dn_rot.scale(1./dn_rot.l2_norm());  // expected surface normal
        dn_rot(0) -= n(0);                  // error in surface normal
        dn_rot(1) -= n(1);
    }
}


inline void
GaussianBumpMesh2D::init (const Real height,
                          const std::vector<MeshInitializer::CoordinateDivisions*>& divs,
                          libMesh::UnstructuredMesh& mesh, libMesh::ElemType t)
{
    libmesh_assert(divs.size() == 2);
    libmesh_assert(divs[0]->n_divs() >= 1);
    libmesh_assert(divs[1]->n_divs() >= 1);
    
    _height = height;
    _x0 = divs[0]->div_location(0);
    _x1 = divs[0]->div_location( divs[0]->n_divs() );
    _y0 = divs[1]->div_location(0);
    _y1 = divs[1]->div_location( divs[1]->n_divs() );
    MeshInitializer::init(divs, mesh, t);
}



inline void
GaussianBumpMesh2D::process_mesh( )
{
    // now move the mesh points
    libMesh::MeshBase::node_iterator   n_it  = _mesh->nodes_begin();
    const libMesh::Mesh::node_iterator n_end = _mesh->nodes_end();
    
    Real x, y;
    
    for (; n_it != n_end; n_it++)
    {
        libMesh::Node& n =  **n_it;
        
        x = n(0);
        y = n(1);
        
        n(1) += (1.0 - (y-_y0)/(_y1-_y0)) * _height *
        exp(-25. * pow( (x-_x0)-(_x1-_x0)*0.5, 2.0));
    }
}




inline void
GaussianBumpMesh3D::init (const Real height,
                          const std::vector<MeshInitializer::CoordinateDivisions*>& divs,
                          libMesh::UnstructuredMesh& mesh, libMesh::ElemType t)
{
    libmesh_assert(divs.size() == 3);
    libmesh_assert(divs[0]->n_divs() >= 1);
    libmesh_assert(divs[1]->n_divs() >= 1);
    libmesh_assert(divs[2]->n_divs() >= 1);
    
    _height = height;
    _x0 = divs[0]->div_location(0);
    _x1 = divs[0]->div_location( divs[0]->n_divs() );
    _y0 = divs[1]->div_location(0);
    _y1 = divs[1]->div_location( divs[1]->n_divs() );
    _z0 = divs[2]->div_location(0);
    _z1 = divs[2]->div_location( divs[2]->n_divs() );
    MeshInitializer::init(divs, mesh, t);
}


inline void
GaussianBumpMesh3D::process_mesh( )
{
    // now move the mesh points
    libMesh::MeshBase::node_iterator   n_it  = _mesh->nodes_begin();
    const libMesh::Mesh::node_iterator n_end = _mesh->nodes_end();
    
    Real x, y, z;
    
    for (; n_it != n_end; n_it++)
    {
        libMesh::Node& n =  **n_it;
        
        x = n(0);
        y = n(1);
        z = n(2);
        
        n(2) += (1.0 - (z-_z0)/(_z1-_z0)) * _height *
        exp(-25. * pow( (x-_x0)-(_x1-_x0)*0.5, 2.0)) *
        exp(-25. * pow( (y-_y0)-(_y1-_y0)*0.5, 2.0));
    }
}



#endif // __MAST_gaussian_bump_mesh_h__
