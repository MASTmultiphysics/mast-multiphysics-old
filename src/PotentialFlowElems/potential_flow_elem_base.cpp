/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2014  Manav Bhatia
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */


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


