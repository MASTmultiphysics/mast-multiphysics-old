//
//  TransientSolverBase.cpp
//  FESystem
//
//  Created by Manav Bhatia on 4/27/12.
//  Copyright (c) 2012. All rights reserved.
//


// FESystem includes
#include "Solvers/TransientSolvers/TransientSolverBase.h"
#include "Base/FESystemExceptions.h"
#include "Base/DegreeOfFreedomMap.h"
#include "Numerics/LocalVector.h"
#include "Base/macros.h"
#include "Numerics/MatrixBase.h"
#include "Numerics/VectorBase.h"
#include "Numerics/SparsityPattern.h"


template <typename ValType>
FESystem::TransientSolvers::TransientSolverBase<ValType>::TransientSolverBase():
if_initialized(false),
if_initialized_initial_time_data(false),
if_identity_mass_matrix(false),
order(0),
n_dofs(0),
initial_time(0.0),
current_time(0.0),
current_time_step(0.0),
current_iteration_number(0),
latest_call_back(FESystem::TransientSolvers::WAITING_TO_START),
state_vec(NULL),
state_velocity(NULL)
{
    
}



template <typename ValType>
FESystem::TransientSolvers::TransientSolverBase<ValType>::~TransientSolverBase()
{
    this->clear();
}


template <typename ValType>
void 
FESystem::TransientSolvers::TransientSolverBase<ValType>::initialize(FESystemUInt o, FESystemUInt n_dofs)
{
    FESystemAssert0(!this->ifInitialized(), FESystem::Exception::InvalidState);
    
    this->order = o;
    this->n_dofs = n_dofs;
    
    this->state_vec = new FESystem::Numerics::LocalVector<ValType>;
    this->state_velocity = new FESystem::Numerics::LocalVector<ValType>;
    this->state_vec->resize(this->order*this->n_dofs);
    this->state_velocity->resize(this->order*this->n_dofs);
    
    this->active_jacobian_terms.resize(o);
    for (FESystemUInt i=0; i<o; i++) this->active_jacobian_terms[i] = false;
    
    
    this->if_initialized = true;
    this->if_initialized_initial_time_data = true;
}




template <typename ValType>
FESystemBoolean 
FESystem::TransientSolvers::TransientSolverBase<ValType>::ifInitialized() const
{
    return (this->if_initialized && this->if_initialized_initial_time_data);
}



template <typename ValType>
void
FESystem::TransientSolvers::TransientSolverBase<ValType>::clear()
{
    this->if_initialized = false;
    this->order = 0;
    this->initial_time = 0.0;
    this->current_time = 0.0;
    this->current_time_step = 0.0;
    this->current_iteration_number = 0;
    this->latest_call_back = FESystem::TransientSolvers::WAITING_TO_START;
    
    // delete the vectors if they have been initialized
    if (this->state_vec != NULL)
        delete this->state_vec;
    if (this->state_velocity != NULL)
        delete this->state_velocity;
    
    this->state_vec = NULL;
    this->state_velocity = NULL;
}


template <typename ValType>
void 
FESystem::TransientSolvers::TransientSolverBase<ValType>::setInitialTimeData(typename RealOperationType(ValType) t0, typename RealOperationType(ValType) dt, FESystem::Numerics::VectorBase<ValType>& vec)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    FESystemAssert2(vec.getSize() == this->order*this->n_dofs, FESystem::Exception::DimensionsDoNotMatch, vec.getSize(), this->order*this->n_dofs);
    
    this->initial_time = t0;
    this->current_time = t0;
    this->current_time_step = dt;
    this->current_iteration_number = 0;
    this->latest_call_back = FESystem::TransientSolvers::WAITING_TO_START;
    
    this->state_vec->copyVector(vec);
    
    this->if_initialized_initial_time_data = true;
}



template <typename ValType>
FESystem::Numerics::VectorBase<ValType>&
FESystem::TransientSolvers::TransientSolverBase<ValType>::getCurrentStateVector()
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    return *(this->state_vec);
}



template <typename ValType>
FESystem::Numerics::VectorBase<ValType>&
FESystem::TransientSolvers::TransientSolverBase<ValType>::getCurrentStateVelocityVector()
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    return *(this->state_velocity);
}


template <typename ValType>
typename RealOperationType(ValType)
FESystem::TransientSolvers::TransientSolverBase<ValType>::getCurrentTime()
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    return this->current_time;
}


template <typename ValType>
typename RealOperationType(ValType)
FESystem::TransientSolvers::TransientSolverBase<ValType>::getCurrentStepSize()
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    return this->current_time_step;
}



template <typename ValType>
FESystemUInt
FESystem::TransientSolvers::TransientSolverBase<ValType>::getCurrentIterationNumber()
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    return this->current_iteration_number;
}



template <typename ValType>
void 
FESystem::TransientSolvers::TransientSolverBase<ValType>::initializeStateVector(FESystem::Numerics::VectorBase<ValType>& vec)
{
    // TODO: revisit for parallel and sparsity
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    vec.resize(this->state_vec->getSize());
}




template <typename ValType>
void
FESystem::TransientSolvers::TransientSolverBase<ValType>::setMassMatrix(FESystemBoolean if_identity, FESystem::Numerics::MatrixBase<ValType>* mass_mat_ptr)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    this->if_identity_mass_matrix = if_identity;
    if (!if_identity)
    { FESystemAssert0(mass_mat_ptr == NULL, FESystem::Exception::InvalidValue); }
    else
        this->mass_matrix = mass_mat_ptr;
}



template <typename ValType>
void
FESystem::TransientSolvers::TransientSolverBase<ValType>::setActiveJacobianTerm(std::vector<FESystemBoolean>& active_terms)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    FESystemAssert2(active_terms.size() == this->order, FESystem::Exception::DimensionsDoNotMatch, active_terms.size(), this->order);
    for (FESystemUInt i=0; i<this->order; i++) this->active_jacobian_terms[i] = active_terms[i];
}



template <typename ValType>
void 
FESystem::TransientSolvers::TransientSolverBase<ValType>::resizeMatrixToJacobianTemplate(FESystem::Numerics::MatrixBase<ValType>& state_jac)
{
    // TODO: revisit for parallel and sparsity
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    state_jac.resize(this->order*this->n_dofs, this->order*this->n_dofs);
}



template <typename ValType>
void 
FESystem::TransientSolvers::TransientSolverBase<ValType>::initializeMatrixToJacobianTemplate(FESystem::Numerics::MatrixBase<ValType>& state_jac)
{
    // TODO: revisit for parallel and sparsity
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    state_jac.zero();
    
    // now set unit value on the diagonals of each submatrix
    for (FESystemUInt j=0; j<this->order-1; j++)
        for (FESystemUInt i=0; i<this->n_dofs; i++)
            state_jac.setVal(j*this->n_dofs+i, (j+1)*this->n_dofs+i, 1.0);
}




template <typename ValType>
void 
FESystem::TransientSolvers::TransientSolverBase<ValType>::initializeMatrixSparsityPatterForSystem(const FESystem::Numerics::SparsityPattern& spatial_sparsity_pattern,  FESystem::Numerics::SparsityPattern& system_pattern) const
{
    // TODO: revisit for parallel and sparsity
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    system_pattern.clear();
    system_pattern.setNDofs(this->order*this->n_dofs);
    std::set<FESystemUInt> nonzeros;
    
    // now add 1 dof for each lower order time derivatives, and copy the sparsity pattern for each lower order coefficient of the function on rhs 
    // now set unit value on the diagonals of each submatrix
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
    
    for (FESystemUInt j=0; j<this->order-1; j++)
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
            }
        }
    
    // finally, reinitialize the sparsity pattern
    system_pattern.reinit();
}




template <typename ValType>
void 
FESystem::TransientSolvers::TransientSolverBase<ValType>::updateVectorValuesForDerivativeOrder(FESystemUInt o, const FESystem::Numerics::VectorBase<ValType>& dofs,
                                                                                      FESystem::Numerics::VectorBase<ValType>& state)
{
    // TODO: revisit for parallel and sparsity
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    FESystemAssert2(o<this->order, FESystem::Exception::IndexOutOfBound, o, this->order);
    
    state.setSubVectorVals(o*this->n_dofs,(o+1)*this->n_dofs-1, 0, this->n_dofs-1, dofs);
}



template <typename ValType>
void 
FESystem::TransientSolvers::TransientSolverBase<ValType>::updateJacobianValuesForDerivativeOrder(FESystemUInt p, FESystemUInt q,
                                                                                        const FESystem::Numerics::MatrixBase<ValType>& dof_jac, 
                                                                                        FESystem::Numerics::MatrixBase<ValType>& state_jac)
{
    // TODO: revisit for parallel and sparsity
    FESystemAssert0(p>0, FESystem::Exception::InvalidValue);
    FESystemAssert2(p<=this->order, FESystem::Exception::IndexOutOfBound, p, this->order+1);
    FESystemAssert2(q<this->order, FESystem::Exception::IndexOutOfBound, q, this->order);
    
    state_jac.setSubMatrixVals((p-1)*this->n_dofs,p*this->n_dofs-1, 
                               q*this->n_dofs,(q+1)*this->n_dofs-1, 
                               0, this->n_dofs-1, 0, this->n_dofs-1, 
                               dof_jac);
}


template <typename ValType>
void 
FESystem::TransientSolvers::TransientSolverBase<ValType>::extractVectorValuesForDerivativeOrder(FESystemUInt o, const FESystem::Numerics::VectorBase<ValType>& state, 
                                                                                       FESystem::Numerics::VectorBase<ValType>& dofs)
{
    // TODO: revisit for parallel and sparsity
    FESystemAssert2(o<this->order, FESystem::Exception::IndexOutOfBound, o, this->order);
    state.getSubVectorVals(o*this->n_dofs,(o+1)*this->n_dofs-1, 0, this->n_dofs-1, dofs);
}



template <typename ValType>
void
FESystem::TransientSolvers::TransientSolverBase<ValType>::copyDerivativeValuesFromStateToVelocityVector(const FESystem::Numerics::VectorBase<ValType>& state, 
                                                                                               FESystem::Numerics::VectorBase<ValType>& velocity) const
{
    // TODO: revisit for parallel and sparsity
    FESystemAssert0(this->ifInitialized(), FESystem::Exception::InvalidState);
    FESystemAssert2(state.getSize() == velocity.getSize(), FESystem::Exception::DimensionsDoNotMatch, state.getSize(), state.getSize());
    FESystemAssert2(state.getSize() == this->order*this->n_dofs, FESystem::Exception::DimensionsDoNotMatch, state.getSize(), this->order*this->n_dofs);
    
    if (this->order>1) // needed only if the order is greater than 1
        velocity.setSubVectorVals(0, (this->order-1)*this->n_dofs-1, this->n_dofs, this->order*this->n_dofs-1, state);
}



template <typename ValType> 
FESystem::TransientSolvers::LinearTransientSolverBase<ValType>::LinearTransientSolverBase():
FESystem::TransientSolvers::TransientSolverBase<ValType>(),
linear_solver(NULL),
if_constant_system_matrices(false)
{
    
}



template <typename ValType> 
FESystem::TransientSolvers::LinearTransientSolverBase<ValType>::~LinearTransientSolverBase()
{
    
}


template <typename ValType> 
void
FESystem::TransientSolvers::LinearTransientSolverBase<ValType>::clear()
{
    this->linear_solver = NULL;
    this->if_constant_system_matrices = false;
    FESystem::TransientSolvers::TransientSolverBase<ValType>::clear();
}


template <typename ValType> 
void 
FESystem::TransientSolvers::LinearTransientSolverBase<ValType>::setLinearSolver(FESystem::LinearSolvers::LinearSolverBase<ValType>& solver, FESystemBoolean if_constant_matrices)
{
    this->linear_solver = &solver;
    this->if_constant_system_matrices = if_constant_matrices;
}



/***************************************************************************************/
// Template instantiations for some generic classes

INSTANTIATE_CLASS_FOR_ALL_DATA_TYPES(FESystem::TransientSolvers::TransientSolverBase);
INSTANTIATE_CLASS_FOR_ALL_DATA_TYPES(FESystem::TransientSolvers::LinearTransientSolverBase);


/***************************************************************************************/


