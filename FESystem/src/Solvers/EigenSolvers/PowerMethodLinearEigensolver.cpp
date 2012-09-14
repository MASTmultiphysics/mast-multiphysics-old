//
//  PowerMethodLinearEigensolver.cpp
//  FESystem
//
//  Created by Manav Bhatia on 3/26/11.
//  Copyright 2011 . All rights reserved.
//

#include "Solvers/EigenSolvers/PowerMethodLinearEigensolver.h"
#include "Numerics/VectorBase.h"
#include "Numerics/MatrixBase.h"
#include "Base/macros.h"

template <typename ValType> 
FESystem::EigenSolvers::PowerMethodLinearEigenSolver<ValType>::PowerMethodLinearEigenSolver():
FESystem::EigenSolvers::LinearEigenSolverBase<ValType>()
{
    
}


template <typename ValType> 
FESystem::EigenSolvers::PowerMethodLinearEigenSolver<ValType>::~PowerMethodLinearEigenSolver()
{
    
}



template <typename ValType> 
void
FESystem::EigenSolvers::PowerMethodLinearEigenSolver<ValType>::solve()
{
    FESystemAssert0(this->matrices_are_set, FESystem::EigenSolvers::MatrixNotSet);
    
    switch (this->getEigenProblemType()) {
        case FESystem::EigenSolvers::HERMITIAN:
        case FESystem::EigenSolvers::NONHERMITIAN:    
            this->completePowerIterations(this->getAMatrix(), *(this->eig_val_vec), *(this->eig_vec_mat_right));
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
            break;
    }
    
}




/***************************************************************************************/
// Template instantiations for some generic classes

INSTANTIATE_CLASS_FOR_ONLY_REAL_DATA_TYPES(FESystem::EigenSolvers::PowerMethodLinearEigenSolver);


/***************************************************************************************/
