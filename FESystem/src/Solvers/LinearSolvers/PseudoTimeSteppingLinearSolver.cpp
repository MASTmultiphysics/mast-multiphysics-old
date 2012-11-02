//
//  PseudoTimeSteppingLinearSolver.cpp
//  FESystem
//
//  Created by Manav Bhatia on 4/24/12.
//  Copyright (c) 2012. All rights reserved.
//


// FESystem includes
#include "Solvers/LinearSolvers/PseudoTimeSteppingLinearSolver.h"
#include "Solvers/TransientSolvers/ExplicitRungeKuttaTransientSolver.h"
#include "Numerics/LocalVector.h"
#include "Numerics/MatrixBase.h"
#include "Base/macros.h"


template <typename ValType> 
FESystem::LinearSolvers::PseudoTimeSteppingLinearSolver<ValType>::PseudoTimeSteppingLinearSolver():
FESystem::LinearSolvers::LinearSolverBase<ValType>(),
tolerance(1.0e-6),
max_iters(10e10)
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

    FESystemDouble conv=0.0, pseudo_time_step=1.0e-4;
    FESystemUInt n_iters=0;
    FESystemBoolean if_converged=false;

    sol.setAllVals(1.0);

    FESystem::TransientSolvers::ExplicitRungeKuttaTransientSolver<ValType> ode_solver;
    ode_solver.initialize(1, this->getSystemMatrix().getSize().first, 4);
    ode_solver.setInitialTimeData(0.0, pseudo_time_step, sol);
    
    // set the initial guess
    
    FESystem::TransientSolvers::TransientSolverCallBack call_back;

    while (!if_converged)
    {
        call_back = ode_solver.incrementTimeStep();
        switch (call_back)
        {
            case FESystem::TransientSolvers::TIME_STEP_CONVERGED:
                break;
                
            case FESystem::TransientSolvers::EVALUATE_X_DOT:
                // the Jacobian is not updated since it is constant with respect to time
            {
                this->getSystemMatrix().rightVectorMultiply(ode_solver.getCurrentStateVector(), ode_solver.getCurrentVelocityFunctionVector()); // Ax
                ode_solver.getCurrentVelocityFunctionVector().add(-1.0, rhs); // Ax-b
                ode_solver.getCurrentVelocityFunctionVector().scale(-1.0); // -(Ax-b)
                conv = ode_solver.getCurrentVelocityFunctionVector().getL2Norm();
            }
                break;
                
            default:
                FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
                break;
        }
        
        n_iters++;
        if ((conv < this->tolerance))// || (n_iters == this->max_iters))
            if_converged = true;
        std::cout << n_iters << " " << conv << std::endl;
    }
    
    sol.copyVector(ode_solver.getCurrentStateVector());
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
