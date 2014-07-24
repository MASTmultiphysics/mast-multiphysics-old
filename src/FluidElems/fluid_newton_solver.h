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
