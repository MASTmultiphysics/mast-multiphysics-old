//
//  Header.h
//  MAST
//
//  Created by Manav Bhatia on 9/13/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_Header_h__
#define __MAST_Header_h__

// C++ includes
#include <memory>


// MAST includes
#include "Mesh/mesh_initializer.h"


class RinglebMesh: public MeshInitializer
{
public:
    RinglebMesh():
    _tc_ratio(0.),
    _x0(0.), _x1(0.),
    _y0(0.), _y1(0.),
    _n_maxima(0),
    _panel_bc_id(0),
    _symmetry_bc_id(0),
    _cos_profile(false)
    { }
    
    /*!
     *   initializes the object with the division for each dimension.
     */
    virtual void init (unsigned int n_x, unsigned int n_y,
                       UnstructuredMesh& mesh, ElemType t);
    
protected:
    
    /*!
     *   move the mesh points to account for a bump, and apply boudnary condition
     */
    virtual void process_mesh( );
    
    /*!
     *   t/c ratio of the panel
     */
    libMesh::Real _tc_ratio;
    
    libMesh::Real _x0, _x1;
    
    libMesh::Real _y0, _y1;
    
    unsigned int _n_maxima;
    
    unsigned int _panel_bc_id, _symmetry_bc_id;
    
    bool _cos_profile;
};


class RinglebSurfaceNormalCorrection: public MAST::SurfaceMotionBase
{
public:
    RinglebSurfaceNormalCorrection(libMesh::Real x0, libMesh::Real x1, libMesh::Real h):
    MAST::SurfaceMotionBase(),
    _x0(x0),
    _x1(x1),
    _h(h)
    { }
    
    virtual ~RinglebSurfaceNormalCorrection() {}
    
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
RinglebSurfaceNormalCorrection::surface_velocity_time_domain(const libMesh::Real t,
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
RinglebMesh::init (unsigned int n_x, unsigned int n_y,
                   UnstructuredMesh& mesh, ElemType t)
{
    std::auto_ptr<MeshInitializer::CoordinateDivisions>
    x_divs (new MeshInitializer::CoordinateDivisions),
    y_divs (new MeshInitializer::CoordinateDivisions);
    std::vector<MeshInitializer::CoordinateDivisions*> divs(2);
    divs[0] = x_divs.get();
    divs[1] = y_divs.get();

    std::vector<libMesh::Real> div_loc(2), dx_relative(2);
    std::vector<unsigned int> n_divs(1);
    
    div_loc[0] = 0.; div_loc[1] = 1.;
    dx_relative[0] = 1.; dx_relative[1] = 1.;

    n_divs[0] = n_x;
    x_divs->init(1, div_loc, dx_relative, n_divs);

    n_divs[0] = n_y;
    y_divs->init(1, div_loc, dx_relative, n_divs);
    
    MeshInitializer::init(divs, mesh, t);
}



inline void
RinglebMesh::process_mesh( )
{
    // now move the mesh points
    MeshBase::node_iterator   n_it  = _mesh->nodes_begin();
    const Mesh::node_iterator n_end = _mesh->nodes_end();
    
    const libMesh::Real pi = acos(-1.);
    libMesh::Real x_val, y_val;
    
    for (; n_it != n_end; n_it++)
    {
        Node& n =  **n_it;
        
        x_val = n(0);
        if (x_val <= 0.5)
            x_val *= 2.;
        else
            x_val = (0.5-x_val)*2.+1.;
        
        libMesh::Real exp_val = 0.2, x_break=0.7;
        if (x_val <= x_break)
            x_val = x_val/x_break*pow(x_break, exp_val);
        else
            x_val = pow(x_val,exp_val);
        
        y_val = n(1);
        libMesh::Real kval, qval, aval, gammaval =1.4, rhoval, pval, Jval;
        kval = y_val * 0.8 + 0.7; // linear variation from 0.7 to 1.5
        qval = x_val * (kval-0.5) + 0.5; // linear variation from 0.5 to kval
        aval = sqrt(1.-0.5*(gammaval-1.)*qval*qval);
        rhoval = pow(aval, 2./(gammaval-1.));
        pval = pow(aval, 2.*gammaval/(gammaval-1.))/gammaval;
        Jval = 1./aval + pow(aval,-3.)/3. + pow(aval,-5.)/5. - 0.5*log((1.+aval)/(1.-aval));
        
        n(1) = sqrt(1.-pow(qval/kval,2.))/kval/rhoval/qval;
        if (n(0) > 0.5)
            n(1) *= -1.;
        n(0) = (2./kval/kval - 1./qval/qval)/2./rhoval - Jval/2.;
    }
    
    
    
}


#endif
