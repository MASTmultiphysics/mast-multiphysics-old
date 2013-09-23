//
//  frequency_domain_linearized_euler.h
//  FESystem
//
//  Created by Manav Bhatia on 3/8/13.
//
//

#ifndef __FESystem__frequency_domain_linearized_euler__
#define __FESystem__frequency_domain_linearized_euler__

// libmesh config includes
#include "libmesh/libmesh_config.h"
#include "libmesh/equation_systems.h"

#ifdef LIBMESH_USE_COMPLEX_NUMBERS

// DiffSystem framework files
#include "libmesh/fem_system.h"
#include "libmesh/auto_ptr.h"
#include "libmesh/mesh_function.h"

// FEsystem includes
#include "FluidElems/fluid_elem_base.h"


using namespace libMesh;

class SurfaceMotionBase;
class FlightCondition;

class FrequencyDomainLinearizedFluidSystem: public FEMSystem, public FluidElemBase
{
public:
    FrequencyDomainLinearizedFluidSystem(EquationSystems& es,
                                         const std::string& name_in,
                                         const unsigned int number_in):
    FEMSystem(es, name_in, number_in),
    FluidElemBase(*es.parameters.get<GetPot*>("input_file")),
    perturbed_surface_motion(NULL),
    _if_localized_sol(false)
    { }
    
    void init_data();
    
    virtual void init_context(DiffContext &context);
    
    
    void localize_fluid_solution();
    
    
    virtual bool element_time_derivative (bool request_jacobian,
                                          DiffContext &context);
    
    virtual bool side_time_derivative (bool request_jacobian,
                                       DiffContext &context);
    
    std::vector<unsigned int> vars;

    /*!
     *   this defines the small disturbance surface motion on top of the 
     *   steady surface motion that the body might have seen.
     */
    SurfaceMotionBase* perturbed_surface_motion;
    
protected:
    
    bool _if_localized_sol;
    
    AutoPtr<NumericVector<Number> > _local_fluid_solution;
    
};




class FrequencyDomainFluidPostProcessSystem : public System
{
public:
    // Constructor
    FrequencyDomainFluidPostProcessSystem(EquationSystems& es,
                                          const std::string& name_in,
                                          const unsigned int number_in)
    : System(es, name_in, number_in)
    {}
    
    virtual void init_data();
    
    virtual void postprocess();
    
    unsigned int u, v, w, T, s, p, cp, a, M, // real part of the variables
                                             // imaginary part of the variables
    rho_im, u_im, v_im, w_im, T_im, s_im, p_im, cp_im, a_im, M_im;

    FlightCondition* flight_condition;
};


#endif // LIBMESH_USE_COMPLEX_NUMBERS

#endif /* defined(__FESystem__frequency_domain_linearized_euler__) */
