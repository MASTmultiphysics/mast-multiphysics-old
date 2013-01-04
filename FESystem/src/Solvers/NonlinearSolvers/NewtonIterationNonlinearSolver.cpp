//
//  NewtonIterationNonlinearSolver.cpp
//  FESystem
//
//  Created by Manav Bhatia on 8/30/12.
//
//

// C++ includes
#include <iomanip>


// FESystem includes
#include "Solvers/NonlinearSolvers/NewtonIterationNonlinearSolver.h"
#include "Solvers/LinearSolvers/LinearSolverBase.h"
#include "Numerics/MatrixBase.h"
#include "Numerics/VectorBase.h"
#include "Numerics/LocalVector.h"
#include "Base/FESystemExceptions.h"
#include "Base/macros.h"


template <typename ValType>
FESystem::NonlinearSolvers::NewtonIterationNonlinearSolver<ValType>::NewtonIterationNonlinearSolver():
FESystem::NonlinearSolvers::NonlinearSolverBase<ValType>(),
base_res_l2(0.0),
base_res_l2_slope(0.0),
newton_step_res_l2(0.0),
if_backtracked(false),
jacobian(NULL),
temp_vec(NULL),
linear_solver(NULL)
{
    
}


template <typename ValType>
FESystem::NonlinearSolvers::NewtonIterationNonlinearSolver<ValType>::~NewtonIterationNonlinearSolver()
{
    if (this->temp_vec != NULL) delete this->temp_vec;
}
            


template <typename ValType>
void
FESystem::NonlinearSolvers::NewtonIterationNonlinearSolver<ValType>::clear()
{
    FESystem::NonlinearSolvers::NonlinearSolverBase<ValType>::clear();
    this->linear_solver = NULL;
    this->jacobian = NULL;
    this->base_res_l2 = 0.0;
    this->base_res_l2_slope = 0.0;
    this->newton_step_res_l2 = 0.0;
    this->if_backtracked = false;
}
            


template <typename ValType>
void
FESystem::NonlinearSolvers::NewtonIterationNonlinearSolver<ValType>::initialize(FESystem::Numerics::MatrixBase<ValType>& mat, FESystem::LinearSolvers::LinearSolverBase<ValType>& solver)
{
    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    FESystemAssert0(s.first == s.second, FESystem::Numerics::MatrixMustBeSquareForOperation);
    
    FESystem::NonlinearSolvers::NonlinearSolverBase<ValType>::initialize(s.first);
    if (this->temp_vec == NULL) this->temp_vec = new FESystem::Numerics::LocalVector<ValType>;
    this->temp_vec->resize(s.first);
    this->jacobian = &mat;
    this->linear_solver = &solver;
    this->linear_solver->setSystemMatrix(*(this->jacobian), this->if_reuse_linear_solver_data_structure);
}
            
            

template <typename ValType>
FESystem::Numerics::MatrixBase<ValType>&
FESystem::NonlinearSolvers::NewtonIterationNonlinearSolver<ValType>::getJacobianMatrix()
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    return *(this->jacobian);
}


template <typename ValType>
FESystem::NonlinearSolvers::NonlinearSolverCallBack
FESystem::NonlinearSolvers::NewtonIterationNonlinearSolver<ValType>::incrementSolution()
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    switch (this->latest_call_back)
    {
        case FESystem::NonlinearSolvers::WAITING_TO_START:
            this->latest_call_back = FESystem::NonlinearSolvers::SET_INITIAL_GUESS;
            break;
            
        case FESystem::NonlinearSolvers::SET_INITIAL_GUESS:
            this->latest_call_back = FESystem::NonlinearSolvers::EVALUATE_RESIDUAL_AND_JACOBIAN;
            break;

        case FESystem::NonlinearSolvers::EVALUATE_RESIDUAL:
        {
            if (this->residual->getL2Norm() <= this->convergence_tolerance)
                this->latest_call_back = FESystem::NonlinearSolvers::SOLUTION_CONVERGED;
            else if (this->current_iteration_number >= this->max_allowed_iterations)
                this->latest_call_back = FESystem::NonlinearSolvers::MAXIMUM_ITERATIONS_REACHED;
            else
            {
                this->newton_step_res_l2 = this->residual->getL2Norm();

                // if the residual L2 norm after full newton iterate is larger than the base norm, correct the step
                if ((this->newton_step_res_l2 > this->base_res_l2)  && !this->if_backtracked)  // Backtracking is allowed only once
                {
                    // calculate the correction in stepsize, and backtrack
                    FESystemDouble a0 = this->base_res_l2, a1 = this->base_res_l2_slope, a2 = this->newton_step_res_l2 - a0 - a1, // quadratic expression of the residual with respect to step size
                    lambda = -a1/2.0/a2;
                    
                    std::cout << "Newton Step Correction: "
                    << "  Base Res Norm: " << std::setw(15) << this->base_res_l2
                    << "  Full dx Res Norm: " << std::setw(15) << this->newton_step_res_l2
                    << "  Corrected step size:  " << std::setw(15) << lambda << std::endl;
                    
                    // now, subtract (1.0-lambda)dx from the current state (which is already -ve)
                    this->sol_vec->add(1.0-lambda, *(this->sol_increment_vec));
                    
                    this->if_backtracked = true;
                    this->latest_call_back = FESystem::NonlinearSolvers::EVALUATE_RESIDUAL;
                    return FESystem::NonlinearSolvers::EVALUATE_RESIDUAL;
                }

                // if it gets past this point, then no backtracking was requested
                this->if_backtracked = false;

                // increment iteration number
                this->current_iteration_number++;

                // otherwise, calculate the next step
                this->latest_call_back = FESystem::NonlinearSolvers::EVALUATE_JACOBIAN;
            }
        }
            break;

            
        case FESystem::NonlinearSolvers::EVALUATE_JACOBIAN:
        case FESystem::NonlinearSolvers::EVALUATE_RESIDUAL_AND_JACOBIAN:
        {
            FESystemDouble res = 0.0, dx = 0.0;
            
            // calculate the solution increment
            if ((this->current_iteration_number == 0) || !this->if_reuse_linear_solver_data_structure)
                this->linear_solver->clear();
            this->linear_solver->setSystemMatrix(*(this->jacobian), this->if_reuse_linear_solver_data_structure);
            this->linear_solver->solve(*(this->residual), *(this->sol_increment_vec));
                        
            // update the solution
            this->sol_vec->add(-1.0, *(this->sol_increment_vec));

            
            // calculate convergence criteria
            res = this->residual->getL2Norm();
            dx = this->sol_increment_vec->getL2Norm();
            
            this->jacobian->rightVectorMultiply(*(this->sol_increment_vec), *(this->temp_vec));
            
            // update ths residual values for backtracking in case the complete Newton step does not lower residual l2 norm
            this->base_res_l2 = this->residual->getL2Norm();
            this->base_res_l2_slope = -this->temp_vec->getL2Norm();

            
            std::cout << "Iter: " << std::setw(10) << this->current_iteration_number
            << "   Residual Norm: " << std::setw(15) << res
            << "   dX Norm: " << std::setw(15) << dx << std::endl;
            
            if ((res <= this->convergence_tolerance) || (dx <= this->convergence_tolerance))
                this->latest_call_back = FESystem::NonlinearSolvers::SOLUTION_CONVERGED;
            else if (this->current_iteration_number >= this->max_allowed_iterations)
                this->latest_call_back = FESystem::NonlinearSolvers::MAXIMUM_ITERATIONS_REACHED;
            else
                this->latest_call_back = FESystem::NonlinearSolvers::EVALUATE_RESIDUAL;
        }
            break;

        default:
            FESystemAssert1(false, FESystem::Exception::EnumerationNotHandled, this->latest_call_back);
            break;
    }
    
    return this->latest_call_back;
}
            


/***************************************************************************************/
// Template instantiations for some generic classes

INSTANTIATE_CLASS_FOR_ALL_DATA_TYPES(FESystem::NonlinearSolvers::NewtonIterationNonlinearSolver);


/***************************************************************************************/


