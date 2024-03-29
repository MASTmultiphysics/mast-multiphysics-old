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
#include "Solvers/residual_based_adaptive_time_solver.h"

#include "libmesh/diff_system.h"
#include "libmesh/euler2_solver.h"
#include "libmesh/numeric_vector.h"

ResidualBaseAdaptiveTimeSolver::ResidualBaseAdaptiveTimeSolver (sys_type& s)
:
AdaptiveTimeSolver(s),
n_iters_per_update(10),
_iter_counter(n_iters_per_update),
growth_exponent(1.2),
min_growth(0.25),
_t_old(s.time),
_x_dot_norm_old(0.),
_first_solve(true)
{
    // add two vectors to the system for storing old solution and
    // velocity estimation
    _x_old = &(s.add_vector("x_old"));
    _x_dot = &(s.add_vector("x_dot"));
    
    // We start with a reasonable time solver: implicit Euler
    core_time_solver.reset(new libMesh::Euler2Solver(s));
}



ResidualBaseAdaptiveTimeSolver::~ResidualBaseAdaptiveTimeSolver ()
{
    
}



void ResidualBaseAdaptiveTimeSolver::solve()
{
    libmesh_assert(core_time_solver.get());
    libmesh_assert(this->n_iters_per_update > 2); // need information from atleast three iterations to adapt

    // set the counter only for the first solve
    Real x_dot_norm = 1.0e10; // an arbitrary norm to begin with
    if (_first_solve) {
        _iter_counter = this->n_iters_per_update;
        _first_solve = false;
    }
    
    // now advance the time step
    core_time_solver->solve();
    core_time_solver->advance_timestep();
    _iter_counter--;
    
    *_x_dot = *_system.solution;
    _x_dot->add(-1., *_x_old);
    _x_dot->scale(1.0/(_system.time - _t_old));
    _x_dot->close();
    x_dot_norm = _x_dot->linfty_norm();
    
    
    if (!quiet)
        libMesh::out << "x_dot norm (old):  " << _x_dot_norm_old
        << "\nx_dot norm (curr.):  " << x_dot_norm << std::endl;
    
    // if the time step update is needed, do that now
    if (_iter_counter == 0)
    {
        if (!quiet)
            libMesh::out << "\n ===  Computing new time step ====" << std::endl;
        
        Real growth_factor = pow( _x_dot_norm_old/x_dot_norm , growth_exponent);
        
        if (growth_factor > this->max_growth)
        {
            growth_factor = this->max_growth;
            if (!quiet)
                libMesh::out << " Growth factor constrained by max allowable: "
                << growth_factor << std::endl;
        }
        
        if (growth_factor < this->min_growth)
        {
            growth_factor = this->min_growth;
            if (!quiet)
                libMesh::out << " Growth factor constrained by min allowable: "
                << growth_factor << std::endl;
        }
        
        this->_system.deltat *= growth_factor;
        
        // look for the maximum and minimum dt
        if (this->_system.deltat > this->max_deltat)
        {
            this->_system.deltat = this->max_deltat;
            if (!quiet)
                libMesh::out << " deltat constrained by max allowable: " << std::endl;
        }
        
        // look for the maximum and minimum dt
        if (this->_system.deltat < this->min_deltat)
        {
            this->_system.deltat = this->min_deltat;
            if (!quiet)
                libMesh::out << " deltat constrained by min allowable: " << std::endl;
        }
        
        libMesh::out << " new deltat :  " << this->_system.deltat << std::endl;
        
        _iter_counter = this->n_iters_per_update;
    }    

    // copy current values to "old"
    _t_old = _system.time;
    _x_dot_norm_old = x_dot_norm;
    *_x_old = *_system.solution;
    _x_old->close();
}

