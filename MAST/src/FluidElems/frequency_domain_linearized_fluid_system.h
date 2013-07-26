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

#ifdef LIBMESH_USE_COMPLEX_NUMBERS

// DiffSystem framework files
#include "libmesh/fem_system.h"
#include "libmesh/auto_ptr.h"
#include "libmesh/mesh_function.h"

// FEsystem includes
#include "FluidElems/fluid_elem_base.h"


using namespace libMesh;

class SurfaceMotion;


class FrequencyDomainLinearizedEuler: public FEMSystem, public EulerElemBase
{
public:
    FrequencyDomainLinearizedEuler(EquationSystems& es,
                                   const std::string& name_in,
                                   const unsigned int number_in):
    FEMSystem(es, name_in, number_in),
    EulerElemBase(),
    _if_localized_sol(false)
    { }
    
    void init_data();
    
    virtual void init_context(DiffContext &context);

    
    void localize_fluid_solution();
    
    
    virtual bool element_time_derivative (bool request_jacobian,
                                          DiffContext &context);
    
    virtual bool side_time_derivative (bool request_jacobian,
                                       DiffContext &context);
    
    AutoPtr<SurfaceMotion> surface_motion;
    
    std::vector<unsigned int> vars;
    
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
};


#endif // LIBMESH_USE_COMPLEX_NUMBERS

#endif /* defined(__FESystem__frequency_domain_linearized_euler__) */
