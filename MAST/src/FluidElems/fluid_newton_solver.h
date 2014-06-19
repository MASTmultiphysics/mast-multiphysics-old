//
//  fluid_newton_solver.h
//  MAST
//
//  Created by Manav Bhatia on 7/12/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_fluid_newton_solver_h__
#define __MAST_fluid_newton_solver_h__

// MAST includes
#include "Base/MAST_data_types.h"

// libMesh includes
#include "libmesh/newton_solver.h"





class FluidNewtonSolver: public libMesh::NewtonSolver
{
public:
    FluidNewtonSolver(sys_type& system);
    
    virtual ~FluidNewtonSolver();

    virtual unsigned int solve ();

protected:

    /*!
     *   This provides a specialized solution update strategy for each degree of freedom
     *   so that each variable continues to stay physically consistent.
     */
    virtual void line_search(Real& current_residual,
                             libMesh::NumericVector<Real> &newton_iterate,
                             const libMesh::NumericVector<Real> &linear_solution);
    
};



#endif // __MAST_fluid_newton_solver_h__
