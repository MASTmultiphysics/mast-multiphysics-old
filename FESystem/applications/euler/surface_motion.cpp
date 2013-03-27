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
plunge_amplitude(0.),
pitch_amplitude(0.),
pitch_phase(0.),
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
    
    plunge_amplitude = 0.;
    pitch_amplitude = 0.;
    pitch_phase = 0.;
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
    
    if (fabs(pitch_amplitude) > 0.)
    {
        // normal distance from pitching axis to the given point
        Number pitch_scale;
        pitch_scale.real() = cos(pitch_phase);
        pitch_scale.imag() = sin(pitch_phase);
        pitch_scale *= pitch_amplitude;
        
        Point r(p); r -= hinge_location;
        Point tmp(pitch_axis.cross(r)); // omega x r

        for (unsigned int i=0; i<vel.size(); i++)
            vel(i) = tmp(i) * pitch_scale;
    }
    
    if (fabs(plunge_amplitude) > 0.)
    {
        for (unsigned int i=0; i<vel.size(); i++)
            vel(i) += plunge_vector(i) * plunge_amplitude;
    }
    
    vel.scale(iota*frequency);
}



void SurfaceMotion::surface_normal_perturbation(const Point& n, DenseVector<Number>& dnormal)
{
    dnormal.zero();
    const Number iota(0., 1.);
    
    if (fabs(pitch_amplitude) > 0.)
    {
        // normal distance from pitching axis to the given point
        Number pitch_scale;
        pitch_scale.real() = cos(pitch_phase);
        pitch_scale.imag() = sin(pitch_phase);
        
        if (dnormal.size() == 2)
        {
            // for small angles, this is the same as dnormal = pitch_amplitude * (pitch_axis x n)
            dnormal(0) = (cos(pitch_amplitude)-1.) * n(0) - sin(pitch_amplitude) * n(1);
            dnormal(1) = sin(pitch_amplitude) * n(0) + (cos(pitch_amplitude)-1.) * n(1);
        }
        
        dnormal.scale(pitch_scale); // will not be multiplied by omega
    }
}


#endif


