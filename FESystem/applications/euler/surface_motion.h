//
//  surface_motion.h
//  FESystem
//
//  Created by Manav Bhatia on 3/20/13.
//
//

#ifndef __FESystem__surface_motion__
#define __FESystem__surface_motion__

// libmesh includes
#include "libmesh/mesh.h"
#include "libmesh/point.h"
#include "libmesh/dense_vector.h"


using namespace libMesh;


#ifdef LIBMESH_USE_COMPLEX_NUMBERS

class SurfaceMotion
{
public:
    SurfaceMotion(MeshBase& m);
    
    ~SurfaceMotion();
    
    MeshBase& mesh;

    Point plunge_vector;
    
    Point pitch_axis;
    
    Point hinge_location;
    
    Real plunge_amplitude;
    
    Real pitch_amplitude;
    
    Real pitch_phase;
    
    Real frequency; // rad/sec
    
    void zero();
    
    void init();

    // calculation in frequency domain
    void surface_velocity(const Point& p, DenseVector<Number>& vel);
    
    void surface_normal_perturbation(const Point& n, DenseVector<Number>& dnormal);
    
private:
    
    
};

#endif // LIBMESH_USE_COMPLEX_NUMBERS

#endif /* defined(__FESystem__surface_motion__) */
