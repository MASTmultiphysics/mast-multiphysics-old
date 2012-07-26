//
//  NewmarkTransientSolver.cpp
//  FESystem
//
//  Created by Manav Bhatia on 4/27/12.
//  Copyright (c) 2012 . All rights reserved.
//

// FESystem includes
#include "Solvers/TransientSolvers/NewmarkTransientSolver.h"
#include "Solvers/LinearSolvers/LinearSolverBase.h"
#include "Numerics/LocalVector.h"
#include "Numerics/DenseMatrix.h"
#include "Base/macros.h"
#include "Numerics/SparsityPattern.h"


template <typename ValType>
FESystem::Solvers::LinearNewmarkTransientSolver<ValType>::LinearNewmarkTransientSolver():
FESystem::Solvers::LinearTransientSolverBase<ValType>(),
if_explicit(false),
jacobian(NULL),
temp_vec(NULL)
{
    
}


template <typename ValType>
FESystem::Solvers::LinearNewmarkTransientSolver<ValType>::~LinearNewmarkTransientSolver()
{
    
}





template <typename ValType>
void
FESystem::Solvers::LinearNewmarkTransientSolver<ValType>::initialize(FESystemUInt o, FESystemUInt n_dofs, const std::vector<FESystemDouble>& int_constants)
{
    // TODO: revisit for parallel and sparse
    FESystem::Solvers::TransientSolverBase<ValType>::initialize(o,n_dofs);
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
        this->temp_vec = new FESystem::Numerics::LocalVector<ValType>;
        this->temp_vec2 = new FESystem::Numerics::LocalVector<ValType>;
        this->temp_vec->resize(this->order*this->n_dofs);
        this->temp_vec2->resize(this->order*this->n_dofs);
    }
    
}




template <typename ValType>
void
FESystem::Solvers::LinearNewmarkTransientSolver<ValType>::clear()
{
    if (this->temp_vec != NULL)
        delete this->temp_vec;

    this->integration_constants.clear();
    this->if_explicit = false;    
    
    this->jacobian = NULL;
    this->temp_vec = NULL;
    FESystem::Solvers::TransientSolverBase<ValType>::clear();
}



template <typename ValType>
void
FESystem::Solvers::LinearNewmarkTransientSolver<ValType>::setJacobianMatrix(FESystem::Numerics::MatrixBase<ValType>& jac)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    FESystemAssert0(!this->if_explicit, FESystem::Exception::InvalidState);
    
    this->jacobian = &jac;
    if (this->if_explicit)
        this->initializeMatrixToJacobianTemplate(*this->jacobian);
}



template <typename ValType>
FESystem::Numerics::MatrixBase<ValType>&
FESystem::Solvers::LinearNewmarkTransientSolver<ValType>::getCurrentJacobianMatrix()
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    FESystemAssert0(!this->if_explicit, FESystem::Exception::InvalidState);
    
    return *(this->jacobian);
}



template <typename ValType>
void 
FESystem::Solvers::LinearNewmarkTransientSolver<ValType>::initializeMatrixSparsityPatterForSystem(const FESystem::Numerics::SparsityPattern& spatial_sparsity_pattern, 
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
FESystem::Solvers::TransientSolverCallBack 
FESystem::Solvers::LinearNewmarkTransientSolver<ValType>::incrementTimeStep()
{
    // make sure this is initialized
    FESystemAssert0(this->ifInitialized(), FESystem::Exception::InvalidState);
    
    // process based on the previous call back and status of solver
    switch (this->latest_call_back) 
    {
        case FESystem::Solvers::TIME_STEP_CONVERGED:
        case FESystem::Solvers::WAITING_TO_START:
            // first time step, setup the matrices
        {
            if (this->if_explicit) 
            {
                // if max of A is zero, then this is an explicit intergration scheme, and get the x_dot only
                // x_dot_n+1 = -r = f_n+1
                // x_n+1 = x_n + dt x_dot_n+1
                this->latest_call_back = FESystem::Solvers::EVALUATE_X_DOT;
                return FESystem::Solvers::EVALUATE_X_DOT;
            }
            else
            {
                FESystemAssert0(this->linear_solver != NULL, FESystem::Exception::InvalidState); 
                // otherwise, an implicit scheme, and get both x_dot and x_dot_jacobian 
                this->latest_call_back = FESystem::Solvers::EVALUATE_X_DOT_AND_X_DOT_JACOBIAN;
                return FESystem::Solvers::EVALUATE_X_DOT_AND_X_DOT_JACOBIAN;
            }
        }
            break;
            
        case FESystem::Solvers::EVALUATE_X_DOT:
        {
            if (this->if_explicit) 
                // if max of A is zero, then this is an explicit intergration scheme, and get the x_dot only
                // x_dot_n+1 = -r = f_n+1
                // x_n+1 = x_n + dt x_dot_n+1
                this->state_vec->add(this->current_time_step, *this->state_velocity);
            else
                // this is the exact x_dot at the end of a time step, with the updated state vector for the time step
                // copy the current state velocity to the temporary vector for use in the solution
                this->temp_vec->copyVector(*this->state_velocity);
            
            // now increment the time step
            this->current_time += this->current_time_step;
            this->current_iteration_number++;
            
            this->latest_call_back = FESystem::Solvers::TIME_STEP_CONVERGED;
            return FESystem::Solvers::TIME_STEP_CONVERGED;

        }
            break;
            
        case FESystem::Solvers::EVALUATE_X_DOT_AND_X_DOT_JACOBIAN:
        {
            // only for implicit scheme, where both x_dot and x_dot_jacobian are used here
            
            if ((this->current_iteration_number == 0) || (!this->if_constant_system_matrices))
            {
                // copy the jacobian terms from 1st row to the 0th row
                switch (this->order) 
                {
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
                        
            switch (this->order) 
            {
                case 2:
                {
                    // prepare the RHS forcing function
                    if (this->if_identity_mass_matrix)
                        this->temp_vec->scaleSubVector(0, this->n_dofs-1, this->current_time_step);
                    else
                    {
                        this->temp_vec2->zero();
                        this->mass_matrix->rightSubVectorMultiply(0, this->n_dofs-1, 0, this->n_dofs-1, 
                                                                  0, this->n_dofs-1, 0, this->n_dofs-1, *(this->temp_vec), *(this->temp_vec2));
                        this->temp_vec2->scale(this->current_time_step);
                        
                        this->temp_vec->setSubVectorVals(0, this->n_dofs-1, 0, this->n_dofs-1, *(this->temp_vec2));
                    }

                    this->temp_vec->addSubVectorVals(0, this->n_dofs-1, this->n_dofs, 2*this->n_dofs-1, 0.5*pow(this->current_time_step,2)*(1.0-2.0*this->integration_constants[0]), *(this->temp_vec));
                    this->temp_vec->scaleSubVector(this->n_dofs, 2*this->n_dofs-1, this->current_time_step*(1.0-this->integration_constants[1]));

                    this->state_velocity->setSubVectorVals(0, this->n_dofs-1, this->n_dofs, 2*this->n_dofs-1, *(this->state_velocity));
                    this->state_velocity->scaleSubVector(0, this->n_dofs-1, pow(this->current_time_step,2) * this->integration_constants[0]);
                    this->state_velocity->scaleSubVector(this->n_dofs, 2*this->n_dofs-1, this->current_time_step * this->integration_constants[1]); 
                    
                    this->state_velocity->add(1.0, *(this->temp_vec));
                    this->temp_vec->zero();
                    
                    this->linear_solver->solve(*this->state_velocity, *this->temp_vec);
                    this->state_vec->add(1.0, *this->temp_vec);

                    // now, evaluate the exact time step
                    this->latest_call_back = FESystem::Solvers::EVALUATE_X_DOT;
                    return FESystem::Solvers::EVALUATE_X_DOT;
                }
                    break;
                    
                default:
                    FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
                    break;
            }
            
            // now increment the time step
            this->current_time += this->current_time_step;
            this->current_iteration_number++;
            
            this->latest_call_back = FESystem::Solvers::TIME_STEP_CONVERGED;
            return FESystem::Solvers::TIME_STEP_CONVERGED;
        }
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
            break;
    }
}



/***************************************************************************************/
// Template instantiations for some generic classes

INSTANTIATE_CLASS_FOR_ALL_DATA_TYPES(FESystem::Solvers::LinearNewmarkTransientSolver);


/***************************************************************************************/

