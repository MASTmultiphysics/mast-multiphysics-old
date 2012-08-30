//
//  NonlinearSolverBase.cpp
//  FESystem
//
//  Created by Manav Bhatia on 8/30/12.
//
//

// FESystem includes
#include "Solvers/NonlinearSolvers/NonlinearSolverBase.h"
#include "Base/macros.h"
#include "Numerics/LocalVector.h"


template <typename ValType>
FESystem::NonlinearSolvers::NonlinearSolverBase<ValType>::NonlinearSolverBase():
if_initialized(false),
max_allowed_iterations(20),
convergence_tolerance(FESystem::Base::getMachineEpsilon<typename RealOperationType(ValType)>()),
n_dofs(0),
current_iteration_number(0),
latest_call_back(WAITING_TO_START),
sol_increment_vec(NULL),
sol_vec(NULL),
residual(NULL)
{
    
}


template <typename ValType>
FESystem::NonlinearSolvers::NonlinearSolverBase<ValType>::~NonlinearSolverBase()
{
    if (this->sol_increment_vec != NULL) delete this->sol_increment_vec;
    if (this->sol_vec != NULL) delete this->sol_vec;
    if (this->residual != NULL) delete this->residual;
}



template <typename ValType>
FESystemBoolean
FESystem::NonlinearSolvers::NonlinearSolverBase<ValType>::ifInitialized() const
{
    return this->if_initialized;
}




template <typename ValType>
void
FESystem::NonlinearSolvers::NonlinearSolverBase<ValType>::setConvergenceLimits(FESystemUInt max_iters, FESystemDouble tol)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    this->max_allowed_iterations = max_iters;
    this->convergence_tolerance = tol;
}



template <typename ValType>
void
FESystem::NonlinearSolvers::NonlinearSolverBase<ValType>::clear()
{
    if (this->sol_increment_vec != NULL) delete this->sol_increment_vec;
    if (this->sol_vec != NULL) delete this->sol_vec;
    if (this->residual != NULL) delete this->residual;
    
    this->sol_increment_vec = NULL;
    this->sol_vec = NULL;
    this->residual = NULL;
    
    this->if_initialized = false;
    this->max_allowed_iterations = 20;
    this->convergence_tolerance = FESystem::Base::getMachineEpsilon<typename RealOperationType(ValType)>();
    this->n_dofs = 0;
    this->current_iteration_number = 0;
    this->latest_call_back = FESystem::NonlinearSolvers::WAITING_TO_START;
}



template <typename ValType>
FESystem::NonlinearSolvers::NonlinearSolverCallBack
FESystem::NonlinearSolvers::NonlinearSolverBase<ValType>::getCurrentCallBack() const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    return this->latest_call_back;
}



template <typename ValType>
FESystem::Numerics::VectorBase<ValType>&
FESystem::NonlinearSolvers::NonlinearSolverBase<ValType>::getCurrentSolution()
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    return *(this->sol_vec);
}



template <typename ValType>
FESystem::Numerics::VectorBase<ValType>&
FESystem::NonlinearSolvers::NonlinearSolverBase<ValType>::getResidualVector()
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    return *(this->residual);
}



template <typename ValType>
unsigned int
FESystem::NonlinearSolvers::NonlinearSolverBase<ValType>::getCurrentIterationNumber()
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    return this->current_iteration_number;
}



template <typename ValType>
void
FESystem::NonlinearSolvers::NonlinearSolverBase<ValType>::initialize(FESystemUInt n)
{
    FESystemAssert0(!this->if_initialized, FESystem::Exception::InvalidState);
    
    this->sol_increment_vec = new FESystem::Numerics::LocalVector<ValType>;
    this->sol_vec = new FESystem::Numerics::LocalVector<ValType>;
    this->residual = new FESystem::Numerics::LocalVector<ValType>;
    
    this->n_dofs = n;
    
    this->sol_increment_vec->resize(n);
    this->sol_vec->resize(n);
    this->residual->resize(n);
    
    this->if_initialized = true;
}


/***************************************************************************************/
// Template instantiations for some generic classes

INSTANTIATE_CLASS_FOR_ALL_DATA_TYPES(FESystem::NonlinearSolvers::NonlinearSolverBase);


/***************************************************************************************/


