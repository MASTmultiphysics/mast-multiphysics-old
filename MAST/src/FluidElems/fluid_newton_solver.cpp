//
//  fluid_newton_solver.cpp
//  MAST
//
//  Created by Manav Bhatia on 7/12/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

// MAST includes
#include "FluidElems/fluid_newton_solver.h"
#include "FluidElems/fluid_system.h"

// libMesh includes
#include "libmesh/numeric_vector.h"
#include "libmesh/diff_system.h"
#include "libmesh/mesh_base.h"




FluidNewtonSolver::FluidNewtonSolver(sys_type& system):
NewtonSolver::NewtonSolver(system)
{}


FluidNewtonSolver::~FluidNewtonSolver()
{}


void
FluidNewtonSolver::line_search(Real& current_residual,
                               libMesh::NumericVector<libMesh::Real> &newton_iterate,
                               const libMesh::NumericVector<libMesh::Real> &linear_solution)
{
  const libMesh::Real cp =  dynamic_cast<const FluidSystem&>
    (this->system()).flight_condition->gas_property.cp;
  const libMesh::Real min_density = 1.0e-3, min_cpT = 1.*cp, max_relative_change = 0.1;
    
    const unsigned int sys_id = this->system().number();
    dof_id_type dof_num;
    
    bool if_lagrange = true;
    for (unsigned int ivar=0; ivar<this->system().n_vars(); ivar++)
        if (this->system().variable_type(ivar).family != LAGRANGE)
            if_lagrange = false;
    
    
    // if the variables in the system use Lagrange elements,
    // then the nodal values are interpolated, and it is enough to
    // control their step size.
    if (if_lagrange)
    {
        libMesh::Real eta=1., var0, dvar, rho, kval;
        // iterate over each nodal dof and set the value
        MeshBase::node_iterator node_it  =
        this->system().get_mesh().pid_nodes_begin(libMesh::global_processor_id());
        MeshBase::node_iterator node_end =
        this->system().get_mesh().pid_nodes_end(libMesh::global_processor_id());
        
        for ( ; node_it != node_end; node_it++)
        {
            kval = 0.; // initialize kinetic energy
            for (unsigned int ivar=0; ivar<this->system().n_vars(); ivar++)
            {
                dof_num = (*node_it)->dof_number(sys_id, ivar, 0);
                var0 = newton_iterate.el(dof_num);
                dvar = linear_solution.el(dof_num);
                
                if (this->system().variable_name(ivar) == "rho")
                {
                    // make sure density remains positive
                    if (var0-dvar > 0)
                        eta = 1.;
                    else
                        eta = (var0-min_density)/dvar;
                    // check the max relative change
                    if (fabs(eta*dvar/var0) > max_relative_change)
                        eta = fabs(var0*max_relative_change/dvar);
                    var0 -= eta*dvar;
                    rho = var0;
                }
                else if (this->system().variable_name(ivar) == "rhoux" ||
                         this->system().variable_name(ivar) == "rhouy" ||
                         this->system().variable_name(ivar) == "rhouz")
                {
                    var0 -= eta*dvar;
                    kval += pow(var0/rho,2)/2.;
                }
                else if (this->system().variable_name(ivar) == "rhoe")
                {
                    // make sure that the kinetic energy is less than total energy
                    // for a positive temperature
                    if ((var0-dvar)/rho-kval > min_cpT)
                        eta = 1.;
                    else
                        eta = (var0-rho*(min_cpT+kval))/dvar;
                    var0 -= eta*dvar;
                }
                
                // update the solution
                newton_iterate.set(dof_num, var0);
            }
        }
    }
    newton_iterate.close();
}



unsigned int
FluidNewtonSolver::solve()
{
    START_LOG("solve()", "NewtonSolver");
    
    // Reset any prior solve result
    _solve_result = INVALID_SOLVE_RESULT;
    
    libMesh::NumericVector<libMesh::Real> &newton_iterate = *(_system.solution);
    
    AutoPtr<libMesh::NumericVector<libMesh::Real> > linear_solution_ptr = newton_iterate.zero_clone();
    libMesh::NumericVector<libMesh::Real> &linear_solution = *linear_solution_ptr;
    libMesh::NumericVector<libMesh::Real> &rhs = *(_system.rhs);
    
    newton_iterate.close();
    linear_solution.close();
    rhs.close();
    
#ifdef LIBMESH_ENABLE_CONSTRAINTS
    _system.get_dof_map().enforce_constraints_exactly(_system);
#endif
    
    SparseMatrix<libMesh::Real> &matrix = *(_system.matrix);
    
    // Prepare to take incomplete steps
    libMesh::Real last_residual=0.;
    
    // Set starting linear tolerance
    libMesh::Real current_linear_tolerance = initial_linear_tolerance;
    
    // Start counting our linear solver steps
    _inner_iterations = 0;
    
    // Now we begin the nonlinear loop
    for (_outer_iterations=0; _outer_iterations<max_nonlinear_iterations;
         ++_outer_iterations)
    {
        if (verbose)
            libMesh::out << "Assembling the System" << std::endl;
        
        _system.assembly(true, true);
        rhs.close();
        libMesh::Real current_residual = rhs.l2_norm();
        last_residual = current_residual;
        
        if (libmesh_isnan(current_residual))
        {
            libMesh::out << "  Nonlinear solver DIVERGED at step "
            << _outer_iterations
            << " with residual Not-a-Number"
            << std::endl;
            libmesh_convergence_failure();
            continue;
        }
        
        max_residual_norm = std::max (current_residual,
                                      max_residual_norm);
        
        // Compute the l2 norm of the whole solution
        libMesh::Real norm_total = newton_iterate.l2_norm();
        
        max_solution_norm = std::max(max_solution_norm, norm_total);
        
        if (verbose)
            libMesh::out << "Nonlinear Residual: "
            << current_residual << std::endl;
        
        // Make sure our linear tolerance is low enough
        current_linear_tolerance = std::min (current_linear_tolerance,
                                             current_residual * linear_tolerance_multiplier);
        
        // But don't let it be too small
        if (current_linear_tolerance < minimum_linear_tolerance)
        {
            current_linear_tolerance = minimum_linear_tolerance;
        }
        
        // At this point newton_iterate is the current guess, and
        // linear_solution is now about to become the NEGATIVE of the next
        // Newton step.
        
        // Our best initial guess for the linear_solution is zero!
        linear_solution.zero();
        
        if (verbose)
            libMesh::out << "Linear solve starting, tolerance "
            << current_linear_tolerance << std::endl;
        
        // Solve the linear system.
        const std::pair<unsigned int, Real> rval =
        linear_solver->solve (matrix, _system.request_matrix("Preconditioner"),
                              linear_solution, rhs, current_linear_tolerance,
                              max_linear_iterations);
        
        // We may need to localize a parallel solution
        _system.update ();
        // The linear solver may not have fit our constraints exactly
#ifdef LIBMESH_ENABLE_CONSTRAINTS
        _system.get_dof_map().enforce_constraints_exactly
        (_system, &linear_solution, /* homogeneous = */ true);
#endif
        
        const unsigned int linear_steps = rval.first;
        libmesh_assert_less_equal (linear_steps, max_linear_iterations);
        _inner_iterations += linear_steps;
        
        const bool linear_solve_finished =
        !(linear_steps == max_linear_iterations);
        
        if (verbose)
            libMesh::out << "Linear solve finished, step " << linear_steps
            << ", residual " << rval.second
            << std::endl;
        
        // Compute the l2 norm of the nonlinear update
        libMesh::Real norm_delta = linear_solution.l2_norm();
        
        if (verbose)
            libMesh::out << "Performing line search" << std::endl;

        // since we're not converged, backtrack if necessary
        this->line_search(current_residual,
                          newton_iterate, linear_solution);

        _system.update();
        
        // Check residual with updated Newton step
        _system.assembly (true, false);
        
        rhs.close();
        current_residual = rhs.l2_norm();
        // check for convergence
        if (verbose)
            libMesh::out << "  Current Residual: "
            << current_residual << std::endl;

        // don't fiddle around if we've already converged
        if (test_convergence(current_residual, norm_delta,
                             linear_solve_finished &&
                             current_residual <= last_residual))
        {
            if (!quiet)
                print_convergence(_outer_iterations, current_residual,
                                  norm_delta, linear_solve_finished &&
                                  current_residual <= last_residual);
            _outer_iterations++;
            break; // out of _outer_iterations for loop
        }
        
        
        // Check to see if backtracking failed,
        // and break out of the nonlinear loop if so...
        if (_solve_result == DiffSolver::DIVERGED_BACKTRACKING_FAILURE)
        {
            _outer_iterations++;
            break; // out of _outer_iterations for loop
        }
        
        if (_outer_iterations + 1 >= max_nonlinear_iterations)
        {
            libMesh::out << "  Nonlinear solver reached maximum step "
            << max_nonlinear_iterations << ", latest evaluated residual "
            << current_residual << std::endl;
            if (continue_after_max_iterations)
            {
                _solve_result = DiffSolver::DIVERGED_MAX_NONLINEAR_ITERATIONS;
                libMesh::out << "  Continuing..." << std::endl;
            }
            else
            {
                libmesh_convergence_failure();
            }
            continue;
        }
        
        // Compute the l2 norm of the whole solution
        norm_total = newton_iterate.l2_norm();
        
        max_solution_norm = std::max(max_solution_norm, norm_total);
        
        // Print out information for the
        // nonlinear iterations.
        if (verbose)
            libMesh::out << "  Nonlinear step: |du|/|u| = "
            << norm_delta / norm_total
            << ", |du| = " << norm_delta
            << std::endl;
    } // end nonlinear loop
    
    // The linear solver may not have fit our constraints exactly
#ifdef LIBMESH_ENABLE_CONSTRAINTS
    _system.get_dof_map().enforce_constraints_exactly(_system);
#endif
    
    // We may need to localize a parallel solution
    _system.update ();
    
    STOP_LOG("solve()", "NewtonSolver");
    
    // Make sure we are returning something sensible as the
    // _solve_result, except in the edge case where we weren't really asked to
    // solve.
    libmesh_assert (_solve_result != DiffSolver::INVALID_SOLVE_RESULT ||
                    !max_nonlinear_iterations);
    
    return _solve_result;
}





