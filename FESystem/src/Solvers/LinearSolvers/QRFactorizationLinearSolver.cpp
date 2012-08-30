//
//  QRFactorizationLinearSolver.cpp
//  FESystem
//
//  Created by Manav Bhatia on 4/12/12.
//  Copyright (c) 2012. All rights reserved.
//


// FEsystem includes
#include "Solvers/LinearSolvers/QRFactorizationLinearSolver.h"
#include "Solvers/Factorizations/HouseHolderTriangulation.h"
#include "Solvers/Factorizations/TriangularBacksubstitution.h"
#include "Numerics/MatrixBase.h"
#include "Numerics/VectorBase.h"
#include "Base/macros.h"


template <typename ValType>
FESystem::LinearSolvers::QRFactorizationLinearSolver<ValType>::QRFactorizationLinearSolver():
FESystem::LinearSolvers::LinearSolverBase<ValType>(),
qr_factorization(new FESystem::FactorizationSolvers::HouseholderTriangulation<ValType>),
triangular_backsubstitute(new FESystem::FactorizationSolvers::TriangularBacksubstitution<ValType>)
{
    
}



template <typename ValType>
FESystem::LinearSolvers::QRFactorizationLinearSolver<ValType>::~QRFactorizationLinearSolver()
{
    
}


template <typename ValType>
void 
FESystem::LinearSolvers::QRFactorizationLinearSolver<ValType>::clear()
{
    this->qr_factorization->clear();
    this->triangular_backsubstitute->clear();
    this->vec.reset();
    this->mat.reset();
    // call the parent's method too
    FESystem::LinearSolvers::LinearSolverBase<ValType>::clear();
}




template <typename ValType>
void 
FESystem::LinearSolvers::QRFactorizationLinearSolver<ValType>::setSystemMatrix(const FESystem::Numerics::MatrixBase<ValType>& mat)
{
    FESystemAssert0(!this->if_initialized, FESystem::Exception::InvalidState);
    
    // and initialize the rest of the data structures
    this->qr_factorization->setMatrix(&mat);
    this->qr_factorization->factorize();
    this->triangular_backsubstitute->setMatrix(this->qr_factorization->getRMatrix());
    this->triangular_backsubstitute->setTriangularMatrixType(FESystem::FactorizationSolvers::UPPER_TRIANGULAR);
    
    FESystem::LinearSolvers::LinearSolverBase<ValType>::setSystemMatrix(mat); 
}





template <typename ValType>
void
FESystem::LinearSolvers::QRFactorizationLinearSolver<ValType>::solve(const FESystem::Numerics::VectorBase<ValType>& rhs,
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

    // Q^T b
    this->qr_factorization->getQMatrix().leftVectorMultiply(rhs,*(this->vec));
    // backsubstitute
    this->triangular_backsubstitute->backSubstitute(*(this->vec), sol);
}




template <typename ValType>
void
FESystem::LinearSolvers::QRFactorizationLinearSolver<ValType>::solve(const FESystem::Numerics::MatrixBase<ValType>& rhs,
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
    
    // Q^T b
    this->qr_factorization->getQMatrix().matrixTransposeRightMultiply(1.0,rhs,*(this->mat));
    // backsubstitute
    this->triangular_backsubstitute->backSubstitute(*(this->mat), sol);
}



/***************************************************************************************/
// Template instantiations for some generic classes

INSTANTIATE_CLASS_FOR_ALL_DATA_TYPES(FESystem::LinearSolvers::QRFactorizationLinearSolver);


/***************************************************************************************/


