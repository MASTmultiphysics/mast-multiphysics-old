//
//  residual_based_adaptive_time_solver.h
//  RealSolver
//
//  Created by Manav Bhatia on 4/29/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __RealSolver__residual_based_adaptive_time_solver__
#define __RealSolver__residual_based_adaptive_time_solver__


// Local includes
#include "libmesh/adaptive_time_solver.h"

// Forward declerations
class EulerSystem;


/**
 * This class wraps another UnsteadySolver derived class, and uses the
 * approximations of the velocity L-infty norm to adapt the time step.
 *
 * @author Manav Bhatia 2013
 */

// ------------------------------------------------------------
// Solver class definition
class ResidualBaseAdaptiveTimeSolver : public AdaptiveTimeSolver
{
public:
    /**
     * The parent class
     */
    typedef AdaptiveTimeSolver Parent;
    
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

    unsigned int n_iters_per_update, _iter_counter;

    Real growth_exponent, min_growth;

    Real _t1, _sol_xinf_t1;

    Real _t2, _sol_xinf_t2;

    Real _t3, _sol_xinf_t3;
    
    Real _xdot_linf_approx;

    bool _first_solve;
};



#endif /* defined(__RealSolver__residual_based_adaptive_time_solver__) */
