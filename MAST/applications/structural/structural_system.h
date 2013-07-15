//
//  structural_system.h
//  RealSolver
//
//  Created by Manav Bhatia on 6/24/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __FESystem__structural_system__
#define __FESystem__structural_system__

#ifndef LIBMESH_USE_COMPLEX_NUMBERS

// DiffSystem framework files
#include "libmesh/fem_system.h"


using namespace libMesh;


Real structural_solution_value(const Point& p,
                               const Parameters& parameters,
                               const std::string& sys_name,
                               const std::string& var_name);



void init_structural_variables(EquationSystems& es,
                               const std::string& system_name);


enum StructuralBoundaryConditionType {TEMPERATURE, PRESSURE, POINT_LOAD};


class StructuralSystem : public FEMSystem
{
public:
    // Constructor
    StructuralSystem(EquationSystems& es,
                     const std::string& name_in,
                     const unsigned int number_in)
    : FEMSystem(es, name_in, number_in)
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
    
    std::vector<unsigned int> _vars;
    
};


#endif // LIBMESH_USE_COMPLEX_NUMBERS


#endif /* defined(__FESystem__structural_system__) */
