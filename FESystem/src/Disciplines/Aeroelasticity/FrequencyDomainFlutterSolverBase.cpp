//
//  FrequencyDomainFlutterSolverBase.cpp
//  FESystem
//
//  Created by Manav Bhatia on 7/24/12.
//  Copyright (c) 2012. All rights reserved.
//

// FESystem includes
#include "Disciplines/Aeroelasticity/FrequencyDomainFlutterSolverBase.h"
#include "Base/FESystemExceptions.h"
#include "Solvers/EigenSolvers/LinearEigenSolverBase.h"
#include "Numerics/DenseMatrix.h"
#include "Base/macros.h"


template <typename ValType>
FESystem::Aeroelasticity::FrequencyDomainFlutterSolverBase<ValType>::FrequencyDomainFlutterSolverBase():
FESystem::Aeroelasticity::FlutterSolverBase<ValType>(),
A_mat(NULL),
B_mat(NULL),
eigen_solver(NULL)
{
    this->A_mat = new FESystem::Numerics::DenseMatrix<ValType>;
    this->B_mat = new FESystem::Numerics::DenseMatrix<ValType>;
}



template <typename ValType>
FESystem::Aeroelasticity::FrequencyDomainFlutterSolverBase<ValType>::~FrequencyDomainFlutterSolverBase()
{
    if (this->A_mat != NULL)
        delete this->A_mat;

    if (this->B_mat != NULL)
        delete this->B_mat;    
}

    

/***************************************************************************************/
// Template instantiations for some generic classes

INSTANTIATE_CLASS_FOR_ONLY_COMPLEX_DATA_TYPES(FESystem::Aeroelasticity::FrequencyDomainFlutterSolverBase);


/***************************************************************************************/
