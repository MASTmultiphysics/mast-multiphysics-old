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

#ifndef __MAST_residual_based_adaptive_time_solver_h__
#define __MAST_residual_based_adaptive_time_solver_h__


// MAST includes
#include "Base/MAST_data_types.h"

// libMesh includes
#include "libmesh/adaptive_time_solver.h"

// Forward declerations
class FluidSystem;


/**
 * This class wraps another libMesh::UnsteadySolver derived class, and uses the
 * approximations of the velocity L-infty norm to adapt the time step.
 *
 * @author Manav Bhatia 2013
 */

// ------------------------------------------------------------
// Solver class definition
class ResidualBaseAdaptiveTimeSolver : public libMesh::AdaptiveTimeSolver
{
public:
    /**
     * The parent class
     */
    typedef libMesh::AdaptiveTimeSolver Parent;
    
    /**
     * Constructor. Requires a reference to the system
     * to be solved.
     */
    explicit
    ResidualBaseAdaptiveTimeSolver (sys_type& s);
    
    /**
     * Destructor.
     */
    ~ResidualBaseAdaptiveTimeSolver ();
    
    void solve();
    
    virtual void init_data() { core_time_solver->init_data(); }

    unsigned int n_iters_per_update, _iter_counter;

    Real growth_exponent, min_growth;

    Real _t_old, _x_dot_norm_old;

    bool _first_solve;
    
    /*!
     *    system solution from previous time-step
     */
    libMesh::NumericVector<Real>* _x_old;
    
    /*!
     *    system velocity estimate
     */
    libMesh::NumericVector<Real>* _x_dot;
};



#endif // __MAST_residual_based_adaptive_time_solver_h__
