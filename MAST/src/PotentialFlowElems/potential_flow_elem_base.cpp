//
//  potential_flow_elem_base.cpp
//  MAST
//
//  Created by Manav Bhatia on 10/7/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//


// MAST includes
#include "PotentialFlowElems/potential_flow_elem_base.h"

void
MAST::PotentialFlowElemBase::init_data ()
{
    // Check the input file for Reynolds number, application type,
    // approximation type
    
    // check if the simulation is viscous or inviscid
    // read the boundary conditions
    unsigned int n_bc, bc_id;
    // first the inviscid conditions
    // slip wall bc
    n_bc = _infile("n_slip_wall_bc", 0);
    if (n_bc > 0)
    {
        for (unsigned int i_bc=0; i_bc<n_bc; i_bc++)
        {
            bc_id = _infile("slip_wall_bc", 0, i_bc);
            _boundary_condition.insert
            (std::multimap<unsigned int, FluidBoundaryConditionType>::value_type
             (bc_id, SLIP_WALL));
        }
    }
    
    
    // symmetry wall bc
    n_bc = _infile("n_symmetry_wall_bc", 0);
    if (n_bc > 0)
    {
        for (unsigned int i_bc=0; i_bc<n_bc; i_bc++)
        {
            bc_id = _infile("symmetry_wall_bc", 0, i_bc);
            _boundary_condition.insert
            (std::multimap<unsigned int, FluidBoundaryConditionType>::value_type
             (bc_id, SYMMETRY_WALL));
        }
    }
    
    
    // far field bc
    n_bc = _infile("n_far_field_bc", 0);
    if (n_bc > 0)
    {
        for (unsigned int i_bc=0; i_bc<n_bc; i_bc++)
        {
            bc_id = _infile("far_field_bc", 0, i_bc);
            _boundary_condition.insert
            (std::multimap<unsigned int, FluidBoundaryConditionType>::value_type
             (bc_id, FAR_FIELD));
        }
    }
}


