//
//  surface_motion.cpp
//  FESystem
//
//  Created by Manav Bhatia on 3/20/13.
//
//

// FESystem includes
#include "euler/surface_motion.h"

#ifdef LIBMESH_USE_COMPLEX_NUMBERS


SurfaceMotion::SurfaceMotion(MeshBase& m):
mesh(m),
plunge_amplitude(Number(0.,0.)),
pitch_amplitude(Number(0.,0.)),
frequency(0.)
{
    
}

SurfaceMotion::~SurfaceMotion()
{
    
}


void SurfaceMotion::zero()
{
    plunge_vector.zero();
    pitch_axis.zero();
    hinge_location.zero();
    
    plunge_amplitude = Number(0.,0.);
    pitch_amplitude = Number(0.,0.);
    frequency = 0.;
}


void SurfaceMotion::init()
{
    // make unit vectors
    if (pitch_axis.size() > 0.)
        pitch_axis /= pitch_axis.size();
    if (plunge_vector.size() > 0.)
        plunge_vector /= plunge_vector.size();
}



void SurfaceMotion::surface_velocity(const Point& p, DenseVector<Number>& vel)
{
    vel.zero();
    const Number iota(0., 1.);
    
    if (abs(pitch_amplitude) > 0.)
    {
        // normal distance from pitching axis to the given point
        Point r(p); r -= hinge_location;
        Point tmp(pitch_axis.cross(r)); // omega x r

        for (unsigned int i=0; i<vel.size(); i++)
            vel(i) += tmp(i) * pitch_amplitude;
    }
    
    if (abs(plunge_amplitude) > 0.)
    {
        for (unsigned int i=0; i<vel.size(); i++)
            vel(i) += plunge_vector(i) * plunge_amplitude;
    }
    
    vel.scale(iota*frequency);
}

#endif


