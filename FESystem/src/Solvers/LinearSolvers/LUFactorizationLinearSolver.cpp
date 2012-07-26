//
//  LUFactorizationLinearSolver.cpp
//  FESystem
//
//  Created by Manav Bhatia on 6/17/12.
//  Copyright (c) 2012. All rights reserved.
//


// FEsystem includes
#include "Solvers/LinearSolvers/LUFactorizationLinearSolver.h"
#include "Solvers/Factorizations/LUFactorization.h"
#include "Solvers/Factorizations/TriangularBacksubstitution.h"
#include "Numerics/MatrixBase.h"
#include "Numerics/VectorBase.h"
#include "Base/macros.h"


template <typename ValType>
FESystem::Solvers::LUFactorizationLinearSolver<ValType>::LUFactorizationLinearSolver():
FESystem::Solvers::LinearSolverBase<ValType>(),
lu_factorization(new FESystem::Solvers::LUFactorization<ValType>),
l_triangular_backsubstitute(new FESystem::Solvers::TriangularBacksubstitution<ValType>),
u_triangular_backsubstitute(new FESystem::Solvers::TriangularBacksubstitution<ValType>)
{
    
}



template <typename ValType>
FESystem::Solvers::LUFactorizationLinearSolver<ValType>::~LUFactorizationLinearSolver()
{
    
}


template <typename ValType>
void 
FESystem::Solvers::LUFactorizationLinearSolver<ValType>::clear()
{
    this->lu_factorization->clear();
    this->l_triangular_backsubstitute->clear();
    this->u_triangular_backsubstitute->clear();
    this->vec.reset();
    this->mat.reset();
    // call the parent's method too
    FESystem::Solvers::LinearSolverBase<ValType>::clear();
}




template <typename ValType>
void 
FESystem::Solvers::LUFactorizationLinearSolver<ValType>::setSystemMatrix(const FESystem::Numerics::MatrixBase<ValType>& mat)
{
    FESystemAssert0(!this->if_initialized, FESystem::Exception::InvalidState);
    
    // and initialize the rest of the data structures
    this->lu_factorization->setMatrix(&mat);
    this->l_triangular_backsubstitute->setMatrix(this->lu_factorization->getLMatrix());
    this->u_triangular_backsubstitute->setMatrix(this->lu_factorization->getUMatrix());

    this->l_triangular_backsubstitute->setTriangularMatrixType(FESystem::Solvers::LOWER_TRIANGULAR);
    this->u_triangular_backsubstitute->setTriangularMatrixType(FESystem::Solvers::UPPER_TRIANGULAR);
    
    FESystem::Solvers::LinearSolverBase<ValType>::setSystemMatrix(mat); 
}





template <typename ValType>
void
FESystem::Solvers::LUFactorizationLinearSolver<ValType>::solve(const FESystem::Numerics::VectorBase<ValType>& rhs,
                                                               FESystem::Numerics::VectorBase<ValType>& sol)
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
    
    
    // initialize the scratch vector if it is not done
    if (this->vec.get() == NULL)
        this->vec.reset(FESystem::Numerics::VectorCreate<ValType>(rhs.getType()).release());
    this->vec->resize(rhs.getSize());

    // L^{-1} b
    this->l_triangular_backsubstitute->backSubstitute(rhs, *(this->vec));
    // U^{-1} b
    this->u_triangular_backsubstitute->backSubstitute(*(this->vec), sol);
    
}




template <typename ValType>
void
FESystem::Solvers::LUFactorizationLinearSolver<ValType>::solve(const FESystem::Numerics::MatrixBase<ValType>& rhs,
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
    
    
    // to be implemented
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}



/***************************************************************************************/
// Template instantiations for some generic classes

INSTANTIATE_CLASS_FOR_ALL_DATA_TYPES(FESystem::Solvers::LUFactorizationLinearSolver);


/***************************************************************************************/


