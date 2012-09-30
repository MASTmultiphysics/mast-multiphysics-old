//
//  NewmarkTransientSolver.cpp
//  FESystem
//
//  Created by Manav Bhatia on 4/27/12.
//  Copyright (c) 2012 . All rights reserved.
//

// C++ includes
#include <iomanip>


// FESystem includes
#include "Solvers/TransientSolvers/NewmarkTransientSolver.h"
#include "Solvers/LinearSolvers/LinearSolverBase.h"
#include "Numerics/LocalVector.h"
#include "Numerics/DenseMatrix.h"
#include "Base/macros.h"
#include "Numerics/SparsityPattern.h"


template <typename ValType>
FESystem::TransientSolvers::NewmarkTransientSolver<ValType>::NewmarkTransientSolver():
FESystem::TransientSolvers::TransientSolverBase<ValType>(),
convergence_tolerance(1.0e-12),
max_nonlinear_iterations(20),
nonlinear_iteration_number(0),
if_explicit(false),
jacobian(NULL),
previous_state(NULL),
previous_velocity(NULL),
residual(NULL),
temp_vec(NULL),
temp_vec2(NULL)
{
    
}


template <typename ValType>
FESystem::TransientSolvers::NewmarkTransientSolver<ValType>::~NewmarkTransientSolver()
{
    
}





template <typename ValType>
void
FESystem::TransientSolvers::NewmarkTransientSolver<ValType>::initialize(FESystemUInt o, FESystemUInt n_dofs, const std::vector<FESystemDouble>& int_constants)
{
    // TODO: revisit for parallel and sparse
    FESystem::TransientSolvers::TransientSolverBase<ValType>::initialize(o,n_dofs);
    this->integration_constants = int_constants;
    
    // if all int_constants are zero, the integration is explicit
    FESystemDouble max_c=0.0;
    for (FESystemUInt i=0; i<int_constants.size(); i++)
        max_c = fmax(max_c, int_constants[i]);
    
    if (max_c == 0.0)
        this->if_explicit = true;
    else
    {
        this->if_explicit = false;
        this->previous_state = new FESystem::Numerics::LocalVector<ValType>; this->previous_state->resize(this->order*this->n_dofs);
        this->previous_velocity = new FESystem::Numerics::LocalVector<ValType>; this->previous_velocity->resize(this->order*this->n_dofs);
        this->residual = new FESystem::Numerics::LocalVector<ValType>; this->residual->resize(this->order*this->n_dofs);
        this->temp_vec = new FESystem::Numerics::LocalVector<ValType>; this->temp_vec->resize(this->order*this->n_dofs);
        this->temp_vec2 = new FESystem::Numerics::LocalVector<ValType>; this->temp_vec2->resize(this->order*this->n_dofs);
    }
    
    this->max_nonlinear_iterations = 20;
    this->convergence_tolerance = 1.0e-12;
    this->nonlinear_iteration_number = 0;
}




template <typename ValType>
void
FESystem::TransientSolvers::NewmarkTransientSolver<ValType>::clear()
{
    if (this->temp_vec != NULL) delete this->temp_vec;
    if (this->temp_vec2 != NULL) delete this->temp_vec2;
    if (this->previous_state != NULL) delete this->previous_state;
    if (this->previous_velocity != NULL) delete this->previous_velocity;
    if (this->residual != NULL) delete this->residual;

    
    this->integration_constants.clear();
    this->if_explicit = false;    
    
    this->jacobian = NULL;
    this->previous_state = NULL;
    this->previous_velocity = NULL;
    this->residual = NULL;
    this->temp_vec = NULL;
    this->temp_vec2 = NULL;

    this->max_nonlinear_iterations = 20;
    this->convergence_tolerance = 1.0e-12;
    this->nonlinear_iteration_number = 0;

    FESystem::TransientSolvers::TransientSolverBase<ValType>::clear();
}



template <typename ValType>
void
FESystem::TransientSolvers::NewmarkTransientSolver<ValType>::setJacobianMatrix(FESystem::Numerics::MatrixBase<ValType>& jac)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    FESystemAssert0(!this->if_explicit, FESystem::Exception::InvalidState);
    
    this->jacobian = &jac;
    if (this->if_explicit)
        this->initializeMatrixToJacobianTemplate(*this->jacobian);
}



template <typename ValType>
FESystem::Numerics::MatrixBase<ValType>&
FESystem::TransientSolvers::NewmarkTransientSolver<ValType>::getCurrentJacobianMatrix()
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    FESystemAssert0(!this->if_explicit, FESystem::Exception::InvalidState);
    
    return *(this->jacobian);
}



template <typename ValType>
void 
FESystem::TransientSolvers::NewmarkTransientSolver<ValType>::initializeMatrixSparsityPatterForSystem(const FESystem::Numerics::SparsityPattern& spatial_sparsity_pattern, 
                                                                                                  FESystem::Numerics::SparsityPattern& system_pattern) const
{
    // TODO: revisit for parallel 
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    system_pattern.clear();
    system_pattern.setNDofs(this->order*this->n_dofs);
    std::set<FESystemUInt> nonzeros;
    
    // now add 1 dof for each lower order time derivatives, and copy the sparsity pattern for each lower order coefficient of the function on rhs 
    // now set unit value on the diagonals of each submatrix
    if (this->if_explicit) // add the diagonal terms for the explicit methods. 
        for (FESystemUInt j=0; j<this->order-1; j++)
        {
            for (FESystemUInt i=0; i<this->n_dofs; i++)
            {
                nonzeros.clear();
                nonzeros.insert(j*this->n_dofs+i); // diagonal term
                nonzeros.insert((j+1)*this->n_dofs+i); // offdiagonal unit term
                system_pattern.addNonZeroColumnsForRow(j*this->n_dofs+i, nonzeros);
            }
        }
    
    // now copy the sparsity pattern for each coefficient term on RHS
    std::vector<std::pair<FESystemUInt, FESystemUInt> >::const_iterator it, end;
    FESystemUInt col_offset=0, row_offset=0;
    
    for (FESystemUInt j=0; j<this->order; j++)
        if (this->active_jacobian_terms[j])
        {
            col_offset = j*this->n_dofs;
            row_offset = (this->order-1)*this->n_dofs;
            for (FESystemUInt i=0; i<this->n_dofs; i++)
            {
                nonzeros.clear();
                
                it = spatial_sparsity_pattern.getAllNonzeroColumnsInRow(i).begin();
                end = spatial_sparsity_pattern.getAllNonzeroColumnsInRow(i).end(); 
                
                for ( ; it != end; it++)
                    nonzeros.insert(it->first+col_offset);
                
                nonzeros.insert((this->order-1)*this->n_dofs+i);// also add a diagonal term
                system_pattern.addNonZeroColumnsForRow(row_offset+i, nonzeros);
                if (!this->if_explicit)
                    system_pattern.addNonZeroColumnsForRow(row_offset+i-this->n_dofs, nonzeros);
            }
        }
    
    // finally, reinitialize the sparsity pattern
    system_pattern.reinit();
}



template <typename ValType>
FESystem::TransientSolvers::TransientSolverCallBack 
FESystem::TransientSolvers::NewmarkTransientSolver<ValType>::incrementTimeStep()
{
    // make sure this is initialized
    FESystemAssert0(this->ifInitialized(), FESystem::Exception::InvalidState);
    
    // process based on the previous call back and status of solver
    switch (this->latest_call_back) 
    {
        case FESystem::TransientSolvers::WAITING_TO_START:
        {
            // for beginning the set of iterations, copy the current state as the previous state
            this->previous_state->copyVector(*(this->current_state));
        }
            
        case FESystem::TransientSolvers::TIME_STEP_CONVERGED:
            // first time step, setup the matrices
        {
            if (this->if_explicit) 
            {
                // if max of A is zero, then this is an explicit intergration scheme, and get the x_dot only
                // x_dot_n+1 = -r = f_n+1
                // x_n+1 = x_n + dt x_dot_n+1
                this->latest_call_back = FESystem::TransientSolvers::EVALUATE_X_DOT;
                return FESystem::TransientSolvers::EVALUATE_X_DOT;
            }
            else
            {
                FESystemAssert0(this->linear_solver != NULL, FESystem::Exception::InvalidState); 
                // otherwise, an implicit scheme, and get both x_dot and x_dot_jacobian 
                this->latest_call_back = FESystem::TransientSolvers::EVALUATE_X_DOT_AND_X_DOT_JACOBIAN;
                return FESystem::TransientSolvers::EVALUATE_X_DOT_AND_X_DOT_JACOBIAN;
            }
        }
            break;
            
        case FESystem::TransientSolvers::EVALUATE_X_DOT:
        {
            if (this->if_explicit) 
                // if max of A is zero, then this is an explicit intergration scheme, and get the x_dot only
                // x_dot_n+1 = -r = f_n+1
                // x_n+1 = x_n + dt x_dot_n+1
                this->current_state->add(this->current_time_step, *this->current_velocity);
            
            // calculate residual to identify convergence
            this->evaluateResidual(*(this->previous_state), *(this->previous_velocity), *(this->current_state), *(this->current_velocity), *(this->residual));
            
            FESystemDouble l2 = this->residual->getL2Norm();
            
            std::cout << "Iter: " << std::setw(10) << this->current_iteration_number
            <<  "  Nonlin Iter: " << std::setw(5) << this->nonlinear_iteration_number << "  Residual: " << std::setw(15) << l2  << std::endl;
            
            
            // if converged, increment the time step
            if ((l2 < this->convergence_tolerance) || (this->nonlinear_iteration_number >= this->max_nonlinear_iterations))
            {
                std::cout << "Convergence Achieved" << std::endl;
                
                // copy the current state and velocity to the previous state and velocity
                this->previous_state->copyVector(*(this->current_state));
                this->previous_velocity->copyVector(*(this->current_velocity));

                this->current_time += this->current_time_step;
                this->current_iteration_number++;
                this->nonlinear_iteration_number = 0;
                
                this->latest_call_back = FESystem::TransientSolvers::TIME_STEP_CONVERGED;
                return FESystem::TransientSolvers::TIME_STEP_CONVERGED;
            }
            else
            {
                this->nonlinear_iteration_number++;
                this->latest_call_back = FESystem::TransientSolvers::EVALUATE_X_DOT_AND_X_DOT_JACOBIAN;
                return FESystem::TransientSolvers::EVALUATE_X_DOT_AND_X_DOT_JACOBIAN;
            }
        }
            break;
            
        case FESystem::TransientSolvers::EVALUATE_X_DOT_AND_X_DOT_JACOBIAN:
        {
            // only for implicit scheme, where both x_dot and x_dot_jacobian are used here
            
            if ((this->current_iteration_number == 0) || (!this->if_constant_system_matrices))
            {
                // copy the jacobian terms from 1st row to the 0th row
                switch (this->order) 
                {
                    case 1:
                    {
                        // scale the jacobian terms
                        this->jacobian->scale(-this->current_time_step*this->integration_constants[0]);
                        
                        // now add the mass term, if present, else shift the diagonal
                        if (this->if_identity_mass_matrix)
                            this->jacobian->shiftDiagonal(1.0);
                        else
                            this->jacobian->add(1.0, *(this->mass_matrix));
                    }
                        break;

                    case 2:
                    {
                        for (FESystemUInt i=0; i<this->order; i++)
                            if (this->active_jacobian_terms[i])
                                this->jacobian->setSubMatrixVals(0, this->n_dofs-1, i*this->n_dofs, (i+1)*this->n_dofs-1, 
                                                                 this->n_dofs, 2*this->n_dofs-1, i*this->n_dofs, (i+1)*this->n_dofs-1,
                                                                 *(this->jacobian));
                        
                        // scale the jacobian terms
                        for (FESystemUInt i=0; i<this->n_dofs; i++)
                        {
                            this->jacobian->scaleRow(i, -pow(this->current_time_step,2) * this->integration_constants[0]);
                            this->jacobian->scaleRow(this->n_dofs+i, -this->current_time_step * this->integration_constants[1]);
                        }

                        // now add the mass term, if present, else shift the diagonal
                        if (this->if_identity_mass_matrix)
                            this->jacobian->shiftDiagonal(1.0);
                        else
                        {
                            this->jacobian->addSubMatrixVals(0, this->n_dofs-1, 0, this->n_dofs-1,  // jacobian indices
                                                             0, this->n_dofs-1, 0, this->n_dofs-1,  // mass matrix indices
                                                             1.0, *(this->mass_matrix)); 
                            this->jacobian->addSubMatrixVals(this->n_dofs, 2*this->n_dofs-1, this->n_dofs, 2*this->n_dofs-1,  // jacobian indices
                                                             0, this->n_dofs-1, 0, this->n_dofs-1,  // mass matrix indices
                                                             1.0, *(this->mass_matrix));
                        }
                    }
                        break;
                        
                    default:
                        FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
                        break;
                }
                
                this->linear_solver->clear();
                this->linear_solver->setSystemMatrix(*this->jacobian);
            }
            
            this->evaluateResidual(*(this->previous_state), *(this->previous_velocity), *(this->current_state), *(this->current_velocity), *(this->residual));
            
            this->linear_solver->solve(*(this->residual), *(this->temp_vec));
            this->current_state->add(-1.0, *(this->temp_vec));
            
            // now, evaluate the exact time step
            this->latest_call_back = FESystem::TransientSolvers::EVALUATE_X_DOT;
            return FESystem::TransientSolvers::EVALUATE_X_DOT;
        }
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
            break;
    }
}






template <typename ValType>
void
FESystem::TransientSolvers::NewmarkTransientSolver<ValType>::evaluateResidual(const FESystem::Numerics::VectorBase<ValType>& prev_state, const FESystem::Numerics::VectorBase<ValType>& prev_velocity,
                                                                              const FESystem::Numerics::VectorBase<ValType>& curr_state, const FESystem::Numerics::VectorBase<ValType>& curr_velocity,
                                                                              FESystem::Numerics::VectorBase<ValType>& res)
{    
    switch (this->order)
    {
        case 1:
        {
            // velocity
            if (!this->if_identity_mass_matrix)
                this->mass_matrix->rightVectorMultiply(prev_velocity, res); // M x_dot_n
            
            res.scale(-this->current_time_step* (1.0 - this->integration_constants[0])); // -dt (1-gamma) M x_dot_n
            res.add(-this->current_time_step*this->integration_constants[0], curr_velocity); // -dt (1-gamma) M x_dot_n - dt gamma x_dot_n+1
            
            // now the states
            this->temp_vec->copyVector(curr_state); // x_n+1
            this->temp_vec->add(-1.0, prev_state); //  x_n+1 - x_n
            
            if (this->if_identity_mass_matrix)
                res.add(1.0, *(this->temp_vec)); // x_n+1 - x_n - dt (1-gamma) M x_dot_n - dt gamma x_dot_n+1;  (with M = 1)
            else
            {
                this->mass_matrix->rightVectorMultiply(*(this->temp_vec), *(this->temp_vec2)); // M_n+1 (x_n+1 - x_n)
                res.add(1.0, *(this->temp_vec2)); // M_n+1 (x_n+1 - x_n) - dt (1-gamma) M x_dot_n - dt gamma x_dot_n+1
            }
        }
            break;
            
        case 2:
        {
            res.copyVector(prev_velocity);

            // prepare the RHS forcing function
            if (this->if_identity_mass_matrix)
                res.scaleSubVector(0, this->n_dofs-1, -this->current_time_step); // [-dt  0; 0  1] {x_dot_n; x_ddot_n}
            else
            {
                this->mass_matrix->rightSubVectorMultiply(0, this->n_dofs-1, 0, this->n_dofs-1,
                                                          0, this->n_dofs-1, 0, this->n_dofs-1, res, *(this->temp_vec)); // [M  0; 0  1] {x_dot_n; x_ddot_n}
                this->temp_vec->scale(-this->current_time_step); // [-dt M  0; 0  1] {x_dot_n; x_ddot_n}
                
                res.setSubVectorVals(0, this->n_dofs-1, 0, this->n_dofs-1, *(this->temp_vec)); // copy to the res vector
            }
            
            res.addSubVectorVals(0, this->n_dofs-1, this->n_dofs, 2*this->n_dofs-1,
                                 -0.5*pow(this->current_time_step,2)*(1.0-2.0*this->integration_constants[0]), res); // [-dt M   -.5 dt^2 (1-2beta); 0  1] {x_dot_n; x_ddot_n}
            res.scaleSubVector(this->n_dofs, 2*this->n_dofs-1, -this->current_time_step*(1.0-this->integration_constants[1])); // [-dt M   -.5 dt^2 (1-2beta); 0  -dt(1-gamma)] {x_dot_n; x_ddot_n}
            
            // now the current velocity term
            this->temp_vec->copyVector(curr_velocity); // x_dot_n+1
            this->temp_vec->setSubVectorVals(0, this->n_dofs-1, this->n_dofs, 2*this->n_dofs-1, *(this->temp_vec)); // [0 1; 0 1] {x_dot_n+1; x_ddot_n+1}
            this->temp_vec->scaleSubVector(0, this->n_dofs-1, pow(this->current_time_step,2) * this->integration_constants[0]);  // [0 dt^2 beta; 0 1] {x_dot_n+1; x_ddot_n+1}
            this->temp_vec->scaleSubVector(this->n_dofs, 2*this->n_dofs-1, this->current_time_step * this->integration_constants[1]); // [0 dt^2 beta; 0 dt gamma] {x_dot_n+1; x_ddot_n+1}

            res.add(-1.0, *(this->temp_vec)); // -[dt M   .5 dt^2 (1-2beta); 0  dt(1-gamma)] {x_dot_n; x_ddot_n} - [0 dt^2 beta; 0 dt gamma] {x_dot_n+1; x_ddot_n+1}

            // now the states
            this->temp_vec->copyVector(curr_state); // {x_n+1; x_dot_n+1}
            this->temp_vec->add(-1.0, prev_state); // {x_n+1; x_dot_n+1} - {x_n; x_dot_n}
            
            if (this->if_identity_mass_matrix)
                res.add(1.0, *(this->temp_vec)); // {x_n+1; x_dot_n+1} - {x_n; x_dot_n} - [dt M  .5 dt^2 (1-2beta); 0  dt(1-gamma)] {x_dot_n; x_ddot_n} - [0 dt^2 beta; 0 dt gamma] {x_dot_n+1; x_ddot_n+1} ; with M = 1
            else
            {
                this->mass_matrix->rightSubVectorMultiply(0, this->n_dofs-1, 0, this->n_dofs-1,
                                                          0, this->n_dofs-1, 0, this->n_dofs-1, *(this->temp_vec), *(this->temp_vec2)); // [M 0; 0 1] ({x_n+1; x_dot_n+1} - {x_n; x_dot_n})
                this->mass_matrix->rightSubVectorMultiply(0, this->n_dofs-1, 0, this->n_dofs-1,
                                                          this->n_dofs, 2*this->n_dofs-1, this->n_dofs, 2*this->n_dofs-1, *(this->temp_vec), *(this->temp_vec2));// [M 0; 0 M] ({x_n+1; x_dot_n+1} - {x_n; x_dot_n})
                res.add(1.0, *(this->temp_vec2)); // [M 0; 0 M] ({x_n+1; x_dot_n+1} - {x_n; x_dot_n}) - [dt M   .5 dt^2 (1-2beta); 0  dt(1-gamma)] {x_dot_n; x_ddot_n} - [0 dt^2 beta; 0 dt gamma] {x_dot_n+1; x_ddot_n+1}
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

INSTANTIATE_CLASS_FOR_ALL_DATA_TYPES(FESystem::TransientSolvers::NewmarkTransientSolver);


/***************************************************************************************/

