/*
 *  LinearLeastSquareSolver.cpp
 *  FESystem
 *
 *  Created by Manav Bhatia on 1/8/11.
 *  Copyright 2011 . All rights reserved.
 *
 */

// FEsystem includes
#include "Solvers/LinearSolvers/LinearLeastSquareSolver.h"
#include "Solvers/LinearSolvers/LinearSolverBase.h"
#include "Numerics/MatrixBase.h"
#include "Numerics/VectorBase.h"
#include "Base/macros.h"


template <typename ValType>
FESystem::LinearSolvers::LinearLeastSquareSolver<ValType>::LinearLeastSquareSolver():
LinearSolverBase<ValType>(),
linear_solver(NULL),
least_square_system_matrix(NULL)
{
    
}


template <typename ValType>
FESystem::LinearSolvers::LinearLeastSquareSolver<ValType>::~LinearLeastSquareSolver()
{
    
}




template <typename ValType>
void 
FESystem::LinearSolvers::LinearLeastSquareSolver<ValType>::clear()
{
    this->linear_solver = NULL;
    this->least_square_system_matrix.reset();
    this->vec.reset();
    this->mat.reset();

    // call the parent's clear routine
    FESystem::LinearSolvers::LinearSolverBase<ValType>::clear();
}


template <typename ValType>
void 
FESystem::LinearSolvers::LinearLeastSquareSolver<ValType>::setLinearSolver(FESystem::LinearSolvers::LinearSolverBase<ValType>& solver)
{
    FESystemAssert0(this->linear_solver == NULL, FESystem::Exception::InvalidState);
    
    this->linear_solver = &solver;
}


template <typename ValType>
FESystem::LinearSolvers::LinearSolverBase<ValType>& 
FESystem::LinearSolvers::LinearLeastSquareSolver<ValType>::getLinearSolver()
{
    FESystemAssert0( this->linear_solver != NULL ,
                    FESystem::LinearSolvers::LinearSolverNotInitialized);
    return *(this->linear_solver);
}



template <typename ValType>
void 
FESystem::LinearSolvers::LinearLeastSquareSolver<ValType>::initializeDataStructures()
{
    FESystemAssert0( this->linear_solver != NULL, FESystem::LinearSolvers::LinearSolverNotInitialized);
        
    std::pair<FESystemUInt, FESystemUInt> s = this->getSystemMatrix().getSize();
    const FESystem::Numerics::MatrixBase<ValType>& m_val = this->getSystemMatrix();
    
    this->least_square_system_matrix.reset(FESystem::Numerics::MatrixCreate<ValType>(m_val.getType()).release());
    this->least_square_system_matrix->resize(s.second, s.second);
    m_val.matrixTransposeRightMultiply(1.0, m_val, *(this->least_square_system_matrix));
    
    this->getLinearSolver().setSystemMatrix(*(this->least_square_system_matrix), false);
}





template <typename ValType>
void
FESystem::LinearSolvers::LinearLeastSquareSolver<ValType>::solve(const FESystem::Numerics::VectorBase<ValType>& rhs,
                                                            FESystem::Numerics::VectorBase<ValType>& sol)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    std::pair<FESystemUInt, FESystemUInt> s = this->getSystemMatrix().getSize();
    
    // the rhs is to be multiplied by the transpose of the system matrix
    FESystemAssert3(rhs.getSize() == s.first, 
                    FESystem::Numerics::MatrixVectorSizeMismatch,
                    s.first, s.second, rhs.getSize()); 
    
    // the sol is to be multiplied by the matrix itelf
    FESystemAssert3(sol.getSize() == s.second, 
                    FESystem::Numerics::MatrixVectorSizeMismatch,
                    s.first, s.second, sol.getSize()); 
    
    // initialize the scratch vector if it is not done
    if (this->vec.get() == NULL)
        this->vec.reset(FESystem::Numerics::VectorCreate<ValType>(rhs.getType()).release());
    vec->resize(s.second);

    this->getSystemMatrix().leftVectorMultiply(rhs, *(this->vec));
    
    this->getLinearSolver().solve(*(this->vec), sol);
}



template <typename ValType>
void
FESystem::LinearSolvers::LinearLeastSquareSolver<ValType>::solve(const FESystem::Numerics::MatrixBase<ValType>& rhs,
                                                            FESystem::Numerics::MatrixBase<ValType>& sol)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    std::pair<FESystemUInt, FESystemUInt> s_rhs = rhs.getSize(), s_sol = sol.getSize(), s_sys = this->system_matrix->getSize();
    
    FESystemAssert4((s_rhs.first == s_sol.first) && (s_rhs.second == s_sol.second),
                    FESystem::Numerics::MatrixSizeMismatch,
                    s_rhs.first,  s_rhs.second,
                    s_sol.first,  s_sol.second);
    
    FESystemAssert4(s_rhs.first == s_sys.first,
                    FESystem::Numerics::MatrixMultiplicationSizeMismatch,
                    s_sys.first,  s_sys.second,
                    s_rhs.first,  s_rhs.second);
    
    
    // initialize the scratch matrix if it is not done
    if (this->mat.get() == NULL)
        this->mat.reset(FESystem::Numerics::MatrixCreate<ValType>(rhs.getType()).release());
    this->mat->resize(s_rhs.first, s_rhs.second);
    
    // A^T b
    this->getSystemMatrix().matrixTransposeRightMultiply(1.0,rhs,*(this->mat));
    // backsubstitute
    this->getLinearSolver().solve(*(this->mat), sol);
}



/***************************************************************************************/
// Template instantiations for some generic classes

INSTANTIATE_CLASS_FOR_ALL_DATA_TYPES(FESystem::LinearSolvers::LinearLeastSquareSolver);


/***************************************************************************************/




