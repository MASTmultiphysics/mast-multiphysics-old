//
//  conduction_system.h
//  FESystem
//
//  Created by Manav Bhatia on 6/24/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __FESystem__conduction_system__
#define __FESystem__conduction_system__

#ifndef LIBMESH_USE_COMPLEX_NUMBERS


// DiffSystem framework files
#include "libmesh/fem_system.h"


using namespace libMesh;


Real conduction_solution_value(const Point& p,
                               const Parameters& parameters,
                               const std::string& sys_name,
                               const std::string& var_name);



void init_conduction_variables(EquationSystems& es,
                               const std::string& system_name);


enum ConductionBoundaryConditionType {TEMPERATURE_DIRICHLET, HEAT_FLUX, SURFACE_CONVECTION_FLUX, SURFACE_RADIATION};


class ConductionSystem : public FEMSystem
{
public:
    // Constructor
    ConductionSystem(EquationSystems& es,
                     const std::string& name_in,
                     const unsigned int number_in)
    : FEMSystem(es, name_in, number_in),
    Tvar(0),
    _k_thermal(0.), _cp_thermal(0.)
    {}
    
    void init_data();
    
    virtual void init_context(DiffContext &context);
    
    virtual bool element_time_derivative (bool request_jacobian,
                                          DiffContext &context);
    
    virtual bool side_time_derivative (bool request_jacobian,
                                       DiffContext &context);
    
    virtual bool mass_residual (bool request_jacobian,
                                DiffContext& context);
    
    
    
protected:
    
    unsigned int  Tvar;
    
    Real _k_thermal, _cp_thermal;
};


#endif // LIBMESH_USE_COMPLEX_NUMBERS


#endif /* defined(__FESystem__conduction_system__) */
