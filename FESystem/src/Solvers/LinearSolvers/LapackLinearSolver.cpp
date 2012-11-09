/*
 *  LapackLinearSolver.cpp
 *  FESystem
 *
 *  Created by Manav Bhatia on 1/4/11.
 *  Copyright 2011 . All rights reserved.
 *
 */

// C++ includes
//#include "clapack.h"

// declaration of lapack routines
extern "C" {
extern int dgetrf_(int*, int*, double*, int*, int*, int*);
extern int dgetrs_(char*, int*, int*, double*, int*, int*, double*, int*, int*);
extern int sgetrf_(int*, int*, float*, int*, int*, int*);
extern int sgetrs_(char*, int*, int*, float*, int*, int*, float*, int*, int*);
}

// FEsystem includes
#include "Solvers/LinearSolvers/LapackLinearSolver.h"
#include "Numerics/MatrixBase.h"
#include "Numerics/VectorBase.h"
#include "Base/macros.h"


template <typename ValType>
FESystem::LinearSolvers::LapackLinearSolver<ValType>::LapackLinearSolver():
LinearSolverBase<ValType>(),
system_matrix_work_copy(NULL)
{
    
}


template <typename ValType>
FESystem::LinearSolvers::LapackLinearSolver<ValType>::~LapackLinearSolver()
{
    
}


template <typename ValType>
void 
FESystem::LinearSolvers::LapackLinearSolver<ValType>::clear()
{
    this->system_matrix_work_copy.clear();
    this->ipiv.clear();
    // call the parent's method too
    FESystem::LinearSolvers::LinearSolverBase<ValType>::clear();
}




template <>
void 
FESystem::LinearSolvers::LapackLinearSolver<FESystemDouble>::initializeDataStructures()
{
    // and initialize the rest of the data structures
    std::pair<FESystemUInt, FESystemUInt> s = this->getSystemMatrix().getSize();
    const FESystem::Numerics::MatrixBase<FESystemDouble>& m_val = this->getSystemMatrix();
    
    this->system_matrix_work_copy.resize(s.first*s.second);
    this->ipiv.resize(s.first);
    
    for (FESystemUInt j=0; j<s.second; j++)
        for (FESystemUInt i=0; i<s.first; i++)
            this->system_matrix_work_copy[s.first * j + i] = m_val.getVal(i,j);
    
    // Factorize matrix
    FESystemInt m=s.first, n=s.second, info = 0;
    dgetrf_(&m, &n, &(this->system_matrix_work_copy[0]), &m, &(this->ipiv[0]), &info);
    
}



template <>
void 
FESystem::LinearSolvers::LapackLinearSolver<FESystemFloat>::initializeDataStructures()
{    
    // and initialize the rest of the data structures
    std::pair<FESystemUInt, FESystemUInt> s = this->getSystemMatrix().getSize();
    const FESystem::Numerics::MatrixBase<FESystemFloat>& m_val = this->getSystemMatrix();
    
    this->system_matrix_work_copy.resize(s.first*s.second);
    this->ipiv.resize(s.first);
    
    for (FESystemUInt j=0; j<s.second; j++)
        for (FESystemUInt i=0; i<s.first; i++)
            this->system_matrix_work_copy[s.first * j + i] = m_val.getVal(i,j);
    
    // Factorize matrix
    FESystemInt m=s.first, n=s.second, info = 0;
    sgetrf_(&m, &n, &(this->system_matrix_work_copy[0]), &m, &(this->ipiv[0]), &info);
    
}





template <>
void
FESystem::LinearSolvers::LapackLinearSolver<FESystemDouble>::solve(const FESystem::Numerics::VectorBase<FESystemDouble>& rhs,
                                                             FESystem::Numerics::VectorBase<FESystemDouble>& sol)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    FESystemAssert3(rhs.getSize() == this->getSystemMatrix().getSize().first, 
                    FESystem::Numerics::MatrixVectorSizeMismatch,
                    this->getSystemMatrix().getSize().first, 
                    this->getSystemMatrix().getSize().second,
                    rhs.getSize()); 
    
    FESystemAssert3(sol.getSize() == this->getSystemMatrix().getSize().first, 
                    FESystem::Numerics::MatrixVectorSizeMismatch,
                    this->getSystemMatrix().getSize().first, 
                    this->getSystemMatrix().getSize().second,
                    sol.getSize()); 
    
    sol.copyVector(rhs);
    FESystemDouble* b = sol.getVectorValues();
    
    // Solve using the factorized matrix and return the solution
    std::pair<FESystemUInt, FESystemUInt> s = this->getSystemMatrix().getSize();
    FESystemInt m=s.first, n=s.second, info = 0, nrhs=1;
    char trans = 'N';
    dgetrs_(&trans, &n, &nrhs, &(this->system_matrix_work_copy[0]), &m, &(this->ipiv[0]), b, &m, &info);
    
}



template <>
void
FESystem::LinearSolvers::LapackLinearSolver<FESystemFloat>::solve(const FESystem::Numerics::VectorBase<FESystemFloat>& rhs,
                                                             FESystem::Numerics::VectorBase<FESystemFloat>& sol)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    FESystemAssert3(rhs.getSize() == this->getSystemMatrix().getSize().first, 
                    FESystem::Numerics::MatrixVectorSizeMismatch,
                    this->getSystemMatrix().getSize().first, 
                    this->getSystemMatrix().getSize().second,
                    rhs.getSize()); 
    
    FESystemAssert3(sol.getSize() == this->getSystemMatrix().getSize().first, 
                    FESystem::Numerics::MatrixVectorSizeMismatch,
                    this->getSystemMatrix().getSize().first, 
                    this->getSystemMatrix().getSize().second,
                    sol.getSize()); 
    
    sol.copyVector(rhs);
    FESystemFloat* b = sol.getVectorValues();
    
    // Solve using the factorized matrix and return the solution
    std::pair<FESystemUInt, FESystemUInt> s = this->getSystemMatrix().getSize();
    FESystemInt m=s.first, n=s.second, info = 0, nrhs=1;
    char trans = 'N';
    sgetrs_(&trans, &n, &nrhs, &(this->system_matrix_work_copy[0]), &m, &(this->ipiv[0]), b, &m, &info);
    
}



template <>
void
FESystem::LinearSolvers::LapackLinearSolver<FESystemDouble>::solve(const FESystem::Numerics::MatrixBase<FESystemDouble>& rhs,
                                                       FESystem::Numerics::MatrixBase<FESystemDouble>& sol)
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
    
    char trans = 'N';
    FESystemInt m=s_sys.first, n=s_sys.second, info = 0, nrhs=1;

    // create a temporary vector
    std::auto_ptr<FESystem::Numerics::VectorBase<FESystemDouble> > vec(FESystem::Numerics::VectorCreate<FESystemDouble>(FESystem::Numerics::LOCAL_VECTOR).release());
    vec->resize(s_rhs.first);
    
    // iterate over each column and solve
    for (FESystemUInt i=0; i<s_rhs.second; i++)
    {
        rhs.getColumnVals(i, 0, s_rhs.first-1, *vec);
        FESystemDouble* b = vec->getVectorValues();
        
        // Solve using the factorized matrix and return the solution
        info = 0; nrhs=1;
        dgetrs_(&trans, &n, &nrhs, &(this->system_matrix_work_copy[0]), &m, &(this->ipiv[0]), b, &m, &info);
        sol.setColumnVals(i,0,s_rhs.first-1,*vec);
    }
}




template <>
void
FESystem::LinearSolvers::LapackLinearSolver<FESystemFloat>::solve(const FESystem::Numerics::MatrixBase<FESystemFloat>& rhs,
                                                             FESystem::Numerics::MatrixBase<FESystemFloat>& sol)
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
    
    char trans = 'N';
    FESystemInt m=s_sys.first, n=s_sys.second, info = 0, nrhs=1;
    
    // create a temporary vector
    std::auto_ptr<FESystem::Numerics::VectorBase<FESystemFloat> > vec(FESystem::Numerics::VectorCreate<FESystemFloat>(FESystem::Numerics::LOCAL_VECTOR).release());
    vec->resize(s_rhs.first);
    
    // iterate over each column and solve
    for (FESystemUInt i=0; i<s_rhs.second; i++)
    {
        rhs.getColumnVals(i, 0, s_rhs.first-1, *vec);
        FESystemFloat* b = vec->getVectorValues();
        
        // Solve using the factorized matrix and return the solution
        info = 0; nrhs=1;
        sgetrs_(&trans, &n, &nrhs, &(this->system_matrix_work_copy[0]), &m, &(this->ipiv[0]), b, &m, &info);
        sol.setColumnVals(i,0,s_rhs.first-1,*vec);
    }
}


/***************************************************************************************/
// Template instantiations for some generic classes

INSTANTIATE_CLASS_FOR_ONLY_REAL_DATA_TYPES(FESystem::LinearSolvers::LapackLinearSolver);


/***************************************************************************************/


