//
//  assembleEuler.h
//  FESystem
//
//  Created by Manav Bhatia on 2/21/13.
//
//

#ifndef __FESystem__assembleEuler__
#define __FESystem__assembleEuler__

// libmesh includes
#include "libmesh/libmesh_config.h"

#ifndef LIBMESH_USE_COMPLEX_NUMBERS

// DiffSystem framework files
#include "libmesh/fem_system.h"

// FESystem includes
#include "FluidElems/fluid_elem_base.h"

using namespace libMesh;


Real euler_solution_value(const Point& p,
                          const Parameters& parameters,
                          const std::string& sys_name,
                          const std::string& var_name);



void init_euler_variables(EquationSystems& es,
                          const std::string& system_name);



class EulerSystem : public FEMSystem, public EulerElemBase
{
public:
    // Constructor
    EulerSystem(EquationSystems& es,
                const std::string& name_in,
                const unsigned int number_in)
    : FEMSystem(es, name_in, number_in),
    EulerElemBase(),
    _rho_norm_old(1.),
    _rho_norm_curr(1.),
    dc_recalculate_tolerance(1.0e-8),
    if_use_stored_dc_coeff(false)
    {}
    
    void init_data();
    
    virtual void init_context(DiffContext &context);
    
    virtual bool element_time_derivative (bool request_jacobian,
                                          DiffContext &context);
    
    virtual bool side_time_derivative (bool request_jacobian,
                                       DiffContext &context);
    
    virtual bool mass_residual (bool request_jacobian,
                                DiffContext& context);
    
    virtual void postprocess();
    

    void evaluate_recalculate_dc_flag();
    
    /*!
     *    tolerance threshold for turning on/off the stored dc coeff flag
     */
    Real dc_recalculate_tolerance;
    
    /*!
     *     flag to tell the system to use the stored discontinuity capturing terms
     */
    bool if_use_stored_dc_coeff;

    std::vector<unsigned int> vars;

protected:
    
    /*!
     *    Current and old norms of density in the flow-field
     */
    Real _rho_norm_old, _rho_norm_curr;
};





class FluidPostProcessSystem : public System
{
public:
    // Constructor
    FluidPostProcessSystem(EquationSystems& es,
                           const std::string& name_in,
                           const unsigned int number_in)
    : System(es, name_in, number_in)
    {}
    
    virtual void init_data();
    
    virtual void postprocess();
    
    unsigned int u, v, w, T, s, p, cp, a, M;
    
    FlightCondition* flight_condition;
};

#endif // LIBMESH_USE_COMPLEX_NUMBERS

#endif /* defined(__FESystem__assembleEuler__) */
