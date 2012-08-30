//
//  PseudoTimeSteppingLinearSolver.cpp
//  FESystem
//
//  Created by Manav Bhatia on 4/24/12.
//  Copyright (c) 2012. All rights reserved.
//


// FESystem includes
#include "Solvers/LinearSolvers/PseudoTimeSteppingLinearSolver.h"
#include "Numerics/LocalVector.h"
#include "Numerics/MatrixBase.h"
#include "Base/macros.h"


template <typename ValType> 
FESystem::LinearSolvers::PseudoTimeSteppingLinearSolver<ValType>::PseudoTimeSteppingLinearSolver():
FESystem::LinearSolvers::LinearSolverBase<ValType>(),
tolerance(1.0e-6),
max_iters(100)
{
    
}

template <typename ValType> 
FESystem::LinearSolvers::PseudoTimeSteppingLinearSolver<ValType>::~PseudoTimeSteppingLinearSolver()
{
    
}
            

template <typename ValType> 
void 
FESystem::LinearSolvers::PseudoTimeSteppingLinearSolver<ValType>::clear()
{
    this->residual_vec.reset();
    // call the parent's method too
    FESystem::LinearSolvers::LinearSolverBase<ValType>::clear();
}
            

template <typename ValType> 
void 
FESystem::LinearSolvers::PseudoTimeSteppingLinearSolver<ValType>::setSystemMatrix(const FESystem::Numerics::MatrixBase<ValType>& mat)
{
    FESystemAssert0(!this->if_initialized, FESystem::Exception::InvalidState);

    this->residual_vec.reset(new FESystem::Numerics::LocalVector<ValType>);
    
    this->residual_vec->resize(mat.getSize().first);
    
    FESystem::LinearSolvers::LinearSolverBase<ValType>::setSystemMatrix(mat); 
}

template <typename ValType> 
void
FESystem::LinearSolvers::PseudoTimeSteppingLinearSolver<ValType>::solve(const FESystem::Numerics::VectorBase<ValType>& rhs,
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

    FESystemDouble conv=0.0, pseudo_time_step=0.001;
    FESystemUInt n_iters=0;
    FESystemBoolean if_converged=false;
    
    // set the initial guess
    sol.setAllVals(1.0);
    
    while (!if_converged)
    {
        this->residual_vec->zero();
        this->getSystemMatrix().rightVectorMultiply(sol, *(this->residual_vec)); // A x 
        this->residual_vec->add(-1.0, rhs); // Ax-b
        conv = this->residual_vec->getL2Norm();
        sol.add(-pseudo_time_step, *(this->residual_vec));
        n_iters++;
        if ((conv < this->tolerance) || (n_iters == this->max_iters))
            if_converged = true;
        std::cout << n_iters << " " << conv << std::endl;
    }
}
            

template <typename ValType> 
void 
FESystem::LinearSolvers::PseudoTimeSteppingLinearSolver<ValType>::solve(const FESystem::Numerics::MatrixBase<ValType>& rhs,
                                                                   FESystem::Numerics::MatrixBase<ValType>& sol)
{
    
}


/***************************************************************************************/
// Template instantiations for some generic classes

INSTANTIATE_CLASS_FOR_ALL_DATA_TYPES(FESystem::LinearSolvers::PseudoTimeSteppingLinearSolver);


/***************************************************************************************/
