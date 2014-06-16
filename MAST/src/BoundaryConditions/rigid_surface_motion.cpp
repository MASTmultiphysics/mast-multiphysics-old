//
//  surface_motion.cpp
//  MAST
//
//  Created by Manav Bhatia on 3/20/13.
//
//

// MAST includes
#include "BoundaryConditions/rigid_surface_motion.h"


MAST::RigidSurfaceMotion::RigidSurfaceMotion():
SurfaceMotionBase(),
plunge_amplitude(0.),
pitch_amplitude(0.),
pitch_phase(0.)
{
    
}

MAST::RigidSurfaceMotion::~RigidSurfaceMotion()
{
    
}


void
MAST::RigidSurfaceMotion::zero()
{
    plunge_vector.zero();
    pitch_axis.zero();
    hinge_location.zero();
    
    plunge_amplitude = 0.;
    pitch_amplitude = 0.;
    pitch_phase = 0.;
    
    SurfaceMotionBase::zero();
}


void
MAST::RigidSurfaceMotion::init(libMesh::Real freq, libMesh::Real phase)
{
    // make unit vectors
    if (pitch_axis.size() > 0.)
        pitch_axis /= pitch_axis.size();
    if (plunge_vector.size() > 0.)
        plunge_vector /= plunge_vector.size();
    
    SurfaceMotionBase::init(freq, phase);
}



void
MAST::RigidSurfaceMotion::surface_velocity(const libMesh::Real t,
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
    const libMesh::Complex iota(0., 1.);
    
    if (fabs(pitch_amplitude) > 0.)
    {
        // normal distance from pitching axis to the given point
        libMesh::Complex pitch_scale;
        pitch_scale = Complex(cos(pitch_phase), sin(pitch_phase));
        pitch_scale *= pitch_amplitude;
        
        libMesh::Point r(p); r -= hinge_location;
        libMesh::Point r_rot(pitch_axis.cross(r)); // omega x r
        for (unsigned int i=0; i<3; i++)
            u_trans(i) = r_rot(i) * pitch_scale;
        

        libMesh::Point n_rotvec(pitch_axis.cross(n));
        for (unsigned int i=0; i<3; i++)
            dn_rot(i) = n_rotvec(i) * pitch_scale;
    }
    
    if (fabs(plunge_amplitude) > 0.)
        for (unsigned int i=0; i<3; i++)
            u_trans(i) += plunge_vector(i) * plunge_amplitude;

    // u_trans is in phase with velocity
    w_trans = u_trans;
    u_trans.scale(iota*frequency);
}



void
MAST::RigidSurfaceMotion::surface_velocity(const libMesh::Real t,
                                           const libMesh::Point& p,
                                           const libMesh::Point& n,
                                           DenseRealVector& w_trans,
                                           DenseRealVector& u_trans,
                                           DenseRealVector& dn_rot)
{
    u_trans.zero();
    dn_rot.zero();

    
    if (fabs(pitch_amplitude) > 0.)
    {
        // normal distance from pitching axis to the given point
        libMesh::Point r(p); r -= hinge_location;
        libMesh::Point r_rot(pitch_axis.cross(r)); // omega x r
        for (unsigned int i=0; i<3; i++)
            u_trans(i) = r_rot(i) * pitch_amplitude *
            cos(frequency*t + pitch_phase + phase_offset); // cosine for velocity
        
        
        libMesh::Point n_rotvec(pitch_axis.cross(n));
        for (unsigned int i=0; i<3; i++)
            dn_rot(i) = n_rotvec(i) * pitch_amplitude *
            sin(frequency*t + pitch_phase + phase_offset); // sine for position
    }
    
    if (fabs(plunge_amplitude) > 0.)
        for (unsigned int i=0; i<3; i++)
            u_trans(i) += plunge_vector(i) * plunge_amplitude *
            cos(frequency*t + phase_offset); // cosine for velocity
    
    // u_trans is in phase with velocity
    u_trans.scale(frequency);
}



