//
//  frequency_domain_linearized_euler.h
//  MAST
//
//  Created by Manav Bhatia on 3/8/13.
//
//

#ifndef __MAST__frequency_domain_linearized_euler__
#define __MAST__frequency_domain_linearized_euler__

// libMesh includes
#include "libmesh/libmesh_config.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fem_system.h"
#include "libmesh/auto_ptr.h"
#include "libmesh/mesh_function.h"

// MAST includes
#include "FluidElems/fluid_elem_base.h"




class SurfaceMotionBase;
class FlightCondition;

class FrequencyDomainLinearizedFluidSystem: public libMesh::FEMSystem, public FluidElemBase
{
public:
    FrequencyDomainLinearizedFluidSystem(libMesh::EquationSystems& es,
                                         const std::string& name_in,
                                         const unsigned int number_in):
    libMesh::FEMSystem(es, name_in, number_in),
    FluidElemBase(*es.parameters.get<GetPot*>("input_file")),
    perturbed_surface_motion(NULL),
    _if_localized_sol(false)
    { }
    
    void init_data();
    
    virtual void init_context(libMesh::DiffContext &context);
    
    
    void localize_fluid_solution();
    
    
    virtual bool element_time_derivative (bool request_jacobian,
                                          libMesh::DiffContext &context);
    
    virtual bool side_time_derivative (bool request_jacobian,
                                       libMesh::DiffContext &context);
    
    std::vector<unsigned int> vars;

    /*!
     *   this defines the small disturbance surface motion on top of the 
     *   steady surface motion that the body might have seen.
     */
    MAST::SurfaceMotionBase* perturbed_surface_motion;
    
protected:
    
    bool _if_localized_sol;
    
    libMesh::AutoPtr<libMesh::NumericVector<Real> > _local_fluid_solution;
    
};




class FrequencyDomainFluidPostProcessSystem : public libMesh::System
{
public:
    // Constructor
    FrequencyDomainFluidPostProcessSystem(libMesh::EquationSystems& es,
                                          const std::string& name_in,
                                          const unsigned int number_in)
    : libMesh::System(es, name_in, number_in)
    {}
    
    virtual void init_data();
    
    virtual void postprocess();
    
    unsigned int u, v, w, T, s, p, cp, a, M, // real part of the variables
                                             // imaginary part of the variables
    rho_im, u_im, v_im, w_im, T_im, s_im, p_im, cp_im, a_im, M_im;

    FlightCondition* flight_condition;
};




#endif /* defined(__MAST__frequency_domain_linearized_euler__) */
