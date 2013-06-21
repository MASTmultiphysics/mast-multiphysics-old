//
//  assembleEuler.h
//  FESystem
//
//  Created by Manav Bhatia on 2/21/13.
//
//

#ifndef __FESystem__conduction_system__
#define __FESystem__conduction_system__


// DiffSystem framework files
#include "libmesh/fem_system.h"

using namespace libMesh;

enum ConductionBoundaryConditionType {
    TEMPERATURE, ADIABATIC,
    HEAT_FLUX,
    CONVECTION_HEAT_FLUX,
    EXTERNAL_RADIATION,
    INTERNAL_RADIATION };


Real temperature_solution_value(const Point& p,
                                const Parameters& parameters,
                                const std::string& sys_name,
                                const std::string& var_name);



void init_temperature_variables(EquationSystems& es,
                                const std::string& system_name);



class ConductionSystem : public FEMSystem
{
public:
    // Constructor
    ConductionSystem(EquationSystems& es,
                     const std::string& name_in,
                     const unsigned int number_in)
    : FEMSystem(es, name_in, number_in),
    _k_thermal(0.),
    _cp_thermal(0.)
    {}
    
    void init_data();
    
    virtual void init_context(DiffContext &context);
    
    virtual bool element_time_derivative (bool request_jacobian,
                                          DiffContext &context);
    
    virtual bool side_time_derivative (bool request_jacobian,
                                       DiffContext &context);
    
    virtual bool mass_residual (bool request_jacobian,
                                DiffContext& context);
    
    unsigned int T_var;

protected:
    
    Real _k_thermal, _cp_thermal;
    
    std::multimap<unsigned int, ConductionBoundaryConditionType> _boundary_condition;
};





#endif /* defined(__FESystem_conduction_system__) */
