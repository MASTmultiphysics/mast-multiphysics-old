//
//  gaussian_bump_mesh.h
//  MAST
//
//  Created by Manav Bhatia on 9/13/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

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
    virtual void init( const libMesh::Real height,
                       const std::vector<MeshInitializer::CoordinateDivisions*>& divs,
                       UnstructuredMesh& mesh, ElemType t);
    
    libMesh::Real x0() const {return _x0;}

    libMesh::Real x1() const {return _x1;}

    libMesh::Real y0() const {return _y0;}

    libMesh::Real y1() const {return _y1;}
    
    libMesh::Real h() const {return _height;}

protected:
    
    /*!
     *   move the mesh points to account for a bump, and apply boudnary condition
     */
    virtual void process_mesh( );

    /*!
     *   Height of the bump
     */
    libMesh::Real _height;
    
    /*!
     *   Extent of computational domain
     */
    libMesh::Real _x0, _x1;
    libMesh::Real _y0, _y1;
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
    virtual void init(const libMesh::Real height,
                      const std::vector<MeshInitializer::CoordinateDivisions*>& divs,
                      UnstructuredMesh& mesh, ElemType t);
    
protected:
    
    /*!
     *   move the mesh points to account for a bump, and apply boudnary condition
     */
    virtual void process_mesh( );
    
    /*!
     *   Height of the bump
     */
    libMesh::Real _height;
    
    /*!
     *   Extent of computational domain
     */
    libMesh::Real _x0, _x1;
    libMesh::Real _y0, _y1;
    libMesh::Real _z0, _z1;
};




class GaussianBumpSurfaceNormalCorrection2D: public MAST::SurfaceMotionBase
{
public:
    GaussianBumpSurfaceNormalCorrection2D(libMesh::Real x0, libMesh::Real x1, libMesh::Real h):
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
    virtual void surface_velocity_frequency_domain(const libMesh::Point& p,
                                                   const libMesh::Point& n,
                                                   libMesh::DenseVector<libMesh::Complex>& u_trans,
                                                   libMesh::DenseVector<libMesh::Complex>& dn_rot)
    { libmesh_error();}
    
    /*!
     *   calculation of surface velocity in time domain. \p u_trans is
     *   the pure translation velocity component, while \p dn_rot defines the
     *   surface normal perturbation
     */
    virtual void surface_velocity_time_domain(const libMesh::Real t,
                                              const libMesh::Point& p,
                                              const libMesh::Point& n,
                                              libMesh::DenseVector<libMesh::Number>& u_trans,
                                              libMesh::DenseVector<libMesh::Number>& dn_rot);
    
    protected:
    
    libMesh::Real _x0, _x1, _h;
};



inline void
GaussianBumpSurfaceNormalCorrection2D::surface_velocity_time_domain(const libMesh::Real t,
                                                                    const libMesh::Point& p,
                                                                    const libMesh::Point& n,
                                                                    libMesh::DenseVector<libMesh::Number>& u_trans,
                                                                    libMesh::DenseVector<libMesh::Number>& dn_rot)
{
    // for the point p, add the correction of surface normal n to dn_rot
    libMesh::Real x = p(0),
    f = -25. * pow((x-_x0)-(_x1-_x0)*0.5, 2.0),
    dfdx = -50.*((x-_x0)-(_x1-_x0)*0.5),
    dydx = _h*exp(f)*dfdx;

    dn_rot(0) = dydx; dn_rot(1) = -1.;
    dn_rot.scale(1./dn_rot.l2_norm());  // expected surface normal
    dn_rot(0) -= n(0);                  // error in surface normal
    dn_rot(1) -= n(1);
}


inline void
GaussianBumpMesh2D::init (const libMesh::Real height,
                          const std::vector<MeshInitializer::CoordinateDivisions*>& divs,
                          UnstructuredMesh& mesh, ElemType t)
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
    MeshBase::node_iterator   n_it  = _mesh->nodes_begin();
    const Mesh::node_iterator n_end = _mesh->nodes_end();
    
    libMesh::Real x, y;
    
    for (; n_it != n_end; n_it++)
    {
        Node& n =  **n_it;
        
        x = n(0);
        y = n(1);
        
        n(1) += (1.0 - (y-_y0)/(_y1-_y0)) * _height *
        exp(-25. * pow( (x-_x0)-(_x1-_x0)*0.5, 2.0));
    }
}




inline void
GaussianBumpMesh3D::init (const libMesh::Real height,
                          const std::vector<MeshInitializer::CoordinateDivisions*>& divs,
                          UnstructuredMesh& mesh, ElemType t)
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
    MeshBase::node_iterator   n_it  = _mesh->nodes_begin();
    const Mesh::node_iterator n_end = _mesh->nodes_end();
    
    libMesh::Real x, y, z;
    
    for (; n_it != n_end; n_it++)
    {
        Node& n =  **n_it;
        
        x = n(0);
        y = n(1);
        z = n(2);
        
        n(2) += (1.0 - (z-_z0)/(_z1-_z0)) * _height *
        exp(-25. * pow( (x-_x0)-(_x1-_x0)*0.5, 2.0)) *
        exp(-25. * pow( (y-_y0)-(_y1-_y0)*0.5, 2.0));
    }
}



#endif // __MAST_gaussian_bump_mesh_h__
