//
//  residual_based_adaptive_time_solver.cpp
//  RealSolver
//
//  Created by Manav Bhatia on 4/29/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#include "solvers/time_solvers/residual_based_adaptive_time_solver.h"

#include "libmesh/diff_system.h"
#include "libmesh/euler_solver.h"
#include "libmesh/numeric_vector.h"

ResidualBaseAdaptiveTimeSolver::ResidualBaseAdaptiveTimeSolver (sys_type& s)
:
AdaptiveTimeSolver(s),
n_iters_per_update(10),
_iter_counter(n_iters_per_update),
growth_exponent(1.2),
min_growth(0.25),
_t1(0.),
_t2(0.),
_t3(0.),
_sol_xinf_t1(0.),
_sol_xinf_t2(0.),
_sol_xinf_t3(0.),
_first_solve(true)
{
    // We start with a reasonable time solver: implicit Euler
    core_time_solver.reset(new EulerSolver(s));
}



ResidualBaseAdaptiveTimeSolver::~ResidualBaseAdaptiveTimeSolver ()
{
    
}



void ResidualBaseAdaptiveTimeSolver::solve()
{
    libmesh_assert(core_time_solver.get());
    libmesh_assert(this->n_iters_per_update > 2); // need information from atleast three iterations to adapt

    // set the counter only for the first solve
    if (_first_solve)
    {
        _iter_counter = this->n_iters_per_update;
        _first_solve = false;
    }
    
    // get the norms for this time step
    _t1 = _t2;
    _t2 = _t3;
    _t3 = this->_system.time;
    _sol_xinf_t1 = _sol_xinf_t2;
    _sol_xinf_t2 = _sol_xinf_t3;
    _sol_xinf_t3 = this->_system.solution->linfty_norm();
    
    // if the time step update is needed, do that now
    if (_iter_counter == 0)
    {
        if (!quiet)
            libMesh::out << "\n ===  Computing new time step ====" << std::endl;

        const Real dx1 = fabs(_sol_xinf_t2 - _sol_xinf_t1)/(_t2 - _t1),
        dx2 = fabs(_sol_xinf_t3 - _sol_xinf_t2)/(_t3 - _t2);
        Real growth_factor = pow( dx1/dx2 , growth_exponent);
        
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

    // now advance the time step
    core_time_solver->solve();
    core_time_solver->advance_timestep();
    _iter_counter--;
}

