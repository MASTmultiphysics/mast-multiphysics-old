/*
 *  LinearSolverBase.cpp
 *  FESystem
 *
 *  Created by Manav Bhatia on 1/4/11.
 *  Copyright 2011 . All rights reserved.
 *
 */

// FEsystem includes
#include "Solvers/LinearSolvers/LinearSolverBase.h"
#include "Base/macros.h"

template <typename ValType>
FESystem::Solvers::LinearSolverBase<ValType>::LinearSolverBase():
if_initialized(false),
system_matrix(NULL)
{
    
}


template <typename ValType>
FESystem::Solvers::LinearSolverBase<ValType>::~LinearSolverBase()
{
    this->system_matrix = NULL;
}


template <typename ValType>
void
FESystem::Solvers::LinearSolverBase<ValType>::clear()
{
    this->system_matrix = NULL;
    this->if_initialized = false;
}


template <typename ValType>
void 
FESystem::Solvers::LinearSolverBase<ValType>::setSystemMatrix(const FESystem::Numerics::MatrixBase<ValType>& mat)
{
    FESystemAssert0(!this->if_initialized, FESystem::Exception::InvalidState);
    
    this->system_matrix = &mat;
    
    this->if_initialized = true;
}



template <typename ValType>
const FESystem::Numerics::MatrixBase<ValType>& 
FESystem::Solvers::LinearSolverBase<ValType>::getSystemMatrix() const
{
    FESystemAssert0( this->system_matrix != NULL ,
                    FESystem::Solvers::SystemMatrixNotInitialized);
    
    return *(this->system_matrix);
}


/***************************************************************************************/
// Template instantiations for some generic classes

INSTANTIATE_CLASS_FOR_ALL_DATA_TYPES(FESystem::Solvers::LinearSolverBase);


/***************************************************************************************/
