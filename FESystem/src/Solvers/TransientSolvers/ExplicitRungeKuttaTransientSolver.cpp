//
//  ExplicitRungeKuttaTransientSolver.cpp
//  FESystem
//
//  Created by Manav Bhatia on 5/29/12.
//  Copyright (c) 2012. All rights reserved.
//

// FESystem includes
#include "Solvers/TransientSolvers/ExplicitRungeKuttaTransientSolver.h"
#include "Solvers/LinearSolvers/LinearSolverBase.h"
#include "Numerics/LocalVector.h"
#include "Numerics/DenseMatrix.h"
#include "Base/macros.h"


template <typename ValType>
FESystem::TransientSolvers::ExplicitRungeKuttaTransientSolver<ValType>::ExplicitRungeKuttaTransientSolver():
FESystem::TransientSolvers::TransientSolverBase<ValType>(),
n_rk_steps_per_time_increment(0),
n_rk_steps_completed(0),
previous_time(0),
previous_state(NULL),
new_state_estimate(NULL)
{
    
}


template <typename ValType>
FESystem::TransientSolvers::ExplicitRungeKuttaTransientSolver<ValType>::~ExplicitRungeKuttaTransientSolver()
{
    
}


template <typename ValType>
void
FESystem::TransientSolvers::ExplicitRungeKuttaTransientSolver<ValType>::initialize(FESystemUInt o, FESystemUInt n_dofs, FESystemUInt n_rk_steps)
{
    // TODO: revisit for parallel and sparse
    FESystem::TransientSolvers::TransientSolverBase<ValType>::initialize(o,n_dofs);
    this->n_rk_steps_per_time_increment = n_rk_steps;
    this->n_rk_steps_completed = 0;
    this->previous_time = 0.0;
    
    switch (n_rk_steps)
    {
        case 4:
        {
            this->sub_step_coefficients_for_final_step.resize(4);
            this->sub_step_coefficients_for_final_step[0] = 1.0/6.0;
            this->sub_step_coefficients_for_final_step[1] = 1.0/3.0;
            this->sub_step_coefficients_for_final_step[2] = 1.0/3.0;
            this->sub_step_coefficients_for_final_step[3] = 1.0/6.0;

            this->sub_step_iterate_coefficients.resize(3);
            this->sub_step_iterate_coefficients[1] = 0.5;
            this->sub_step_iterate_coefficients[2] = 0.5;
            this->sub_step_iterate_coefficients[3] = 1.0;
            
        }
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
            break;
    }
    
    this->previous_state = new FESystem::Numerics::LocalVector<ValType>;
    this->previous_state->resize(this->order*this->n_dofs);
    this->new_state_estimate = new FESystem::Numerics::LocalVector<ValType>;
    this->new_state_estimate->resize(this->order*this->n_dofs);
    
}



template <typename ValType>
void
FESystem::TransientSolvers::ExplicitRungeKuttaTransientSolver<ValType>::clear()
{
    if (this->previous_state != NULL)
        delete this->previous_state;
    if (this->new_state_estimate != NULL)
        delete this->new_state_estimate;
    
    this->n_rk_steps_per_time_increment = 0;
    this->n_rk_steps_completed = 0;
    this->previous_time = 0.0;

    this->sub_step_coefficients_for_final_step.clear();
    this->sub_step_iterate_coefficients.clear();
    
    this->previous_state = NULL;
    this->new_state_estimate = NULL;
    FESystem::TransientSolvers::TransientSolverBase<ValType>::clear();
}



template <typename ValType>
FESystem::TransientSolvers::TransientSolverCallBack 
FESystem::TransientSolvers::ExplicitRungeKuttaTransientSolver<ValType>::incrementTimeStep()
{
    // make sure this is initialized
    FESystemAssert0(this->ifInitialized(), FESystem::Exception::InvalidState);
    
    // process based on the previous call back and status of solver
    switch (this->latest_call_back) 
    {
        case FESystem::TransientSolvers::TIME_STEP_CONVERGED:
        case FESystem::TransientSolvers::WAITING_TO_START:
            // first time step, setup the matrices
        {
            this->previous_state->copyVector(*this->current_state);
            this->new_state_estimate->copyVector(*this->current_state);
            
            this->current_velocity_function->zero();
            this->n_rk_steps_completed = 0;
            this->latest_call_back = FESystem::TransientSolvers::EVALUATE_X_DOT;
            return FESystem::TransientSolvers::EVALUATE_X_DOT;
        }
            break;
            
        case FESystem::TransientSolvers::EVALUATE_X_DOT:
        {
            // use the sub-iteration data to calculate the next step
            this->new_state_estimate->add(this->current_time_step*sub_step_coefficients_for_final_step[this->n_rk_steps_completed], *this->current_velocity_function);
                        
            if (this->n_rk_steps_completed < (this->n_rk_steps_per_time_increment-1)) // still in the current time step
            {
                this->current_state->copyVector(*this->previous_state);
                this->current_state->add(this->current_time_step*this->sub_step_iterate_coefficients[this->n_rk_steps_completed], *this->current_velocity_function);
                this->current_velocity_function->zero();

                // increment the number of time steps 
                this->current_time = this->previous_time + this->current_time_step * this->sub_step_iterate_coefficients[this->n_rk_steps_completed];
                this->n_rk_steps_completed++;

                this->latest_call_back = FESystem::TransientSolvers::EVALUATE_X_DOT;
                return FESystem::TransientSolvers::EVALUATE_X_DOT;                
            }
            else  // increment to the next time step
            {
                this->current_state->copyVector(*this->new_state_estimate);

                this->current_time = this->previous_time + this->current_time_step;
                this->previous_time = this->current_time;
                this->current_iteration_number++;

                this->latest_call_back = FESystem::TransientSolvers::TIME_STEP_CONVERGED;
                return FESystem::TransientSolvers::TIME_STEP_CONVERGED;
            }
        }
            break;
                        
        default:
            FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
            break;
    }
}



/***************************************************************************************/
// Template instantiations for some generic classes

INSTANTIATE_CLASS_FOR_ALL_DATA_TYPES(FESystem::TransientSolvers::ExplicitRungeKuttaTransientSolver);


/***************************************************************************************/

