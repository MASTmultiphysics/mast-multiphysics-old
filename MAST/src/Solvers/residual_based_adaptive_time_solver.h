//
//  residual_based_adaptive_time_solver.h
//  MAST
//
//  Created by Manav Bhatia on 4/29/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_residual_based_adaptive_time_solver_h__
#define __MAST_residual_based_adaptive_time_solver_h__


// Local includes
#include "libmesh/adaptive_time_solver.h"

// Forward declerations
class FluidSystem;


/**
 * This class wraps another UnsteadySolver derived class, and uses the
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

    libMesh::Real growth_exponent, min_growth;

    libMesh::Real _t_old, _x_dot_norm_old;

    bool _first_solve;
    
    /*!
     *    system solution from previous time-step
     */
    libMesh::NumericVector<libMesh::Number>* _x_old;
    
    /*!
     *    system velocity estimate
     */
    libMesh::NumericVector<libMesh::Number>* _x_dot;
};



#endif // __MAST_residual_based_adaptive_time_solver_h__
