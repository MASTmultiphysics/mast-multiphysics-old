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
FESystem::LinearSolvers::LinearSolverBase<ValType>::LinearSolverBase():
if_initialized(false),
if_reuse_data(false),
system_matrix(NULL)
{
    
}


template <typename ValType>
FESystem::LinearSolvers::LinearSolverBase<ValType>::~LinearSolverBase()
{
    this->system_matrix = NULL;
}


template <typename ValType>
void
FESystem::LinearSolvers::LinearSolverBase<ValType>::clear()
{
    this->system_matrix = NULL;
    this->if_initialized = false;
    this->if_reuse_data = false;
}


template <typename ValType>
void 
FESystem::LinearSolvers::LinearSolverBase<ValType>::setSystemMatrix(const FESystem::Numerics::MatrixBase<ValType>& mat, FESystemBoolean if_reuse_data_structure)
{
    FESystemAssert0(!this->if_initialized || if_reuse_data_structure, FESystem::Exception::InvalidState);
    
    this->system_matrix = &mat;
    this->if_reuse_data = if_reuse_data_structure;
    
    this->initializeDataStructures();
    
    this->if_initialized = true;
}



template <typename ValType>
const FESystem::Numerics::MatrixBase<ValType>& 
FESystem::LinearSolvers::LinearSolverBase<ValType>::getSystemMatrix() const
{
    FESystemAssert0( this->system_matrix != NULL ,
                    FESystem::LinearSolvers::SystemMatrixNotInitialized);
    
    return *(this->system_matrix);
}


/***************************************************************************************/
// Template instantiations for some generic classes

INSTANTIATE_CLASS_FOR_ALL_DATA_TYPES(FESystem::LinearSolvers::LinearSolverBase);


/***************************************************************************************/
