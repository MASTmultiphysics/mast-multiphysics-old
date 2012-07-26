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
if_initialized(false),
n_basis(0),
fluid_rho(0.0),
aero_b_ref(0.0),
current_k_ref(0.0),
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
//    if (this->eigen_solver != NULL)
//        this->eigen_solver->clear();
    
    if (this->A_mat != NULL)
        delete this->A_mat;

    if (this->B_mat != NULL)
        delete this->B_mat;    
}



template <typename ValType>
void
FESystem::Aeroelasticity::FrequencyDomainFlutterSolverBase<ValType>::initialize(FESystemUInt n, 
                                                                                typename RealOperationType(ValType) b_ref, 
                                                                                typename RealOperationType(ValType) rho, 
                                                                                FESystem::Solvers::LinearEigenSolverBase<ValType>& es)
{
    FESystemAssert0(!this->if_initialized, FESystem::Exception::InvalidState);
    this->n_basis = n;
    this->fluid_rho = rho;
    this->aero_b_ref = b_ref;
    this->eigen_solver = &es;

//    this->eigen_solver->clear();
    
    // initializes the two matrices for use as the two matrices of the generalized non-Hermitian eigensolution
    this->initializeMatrices(n, *(this->A_mat), *(this->B_mat));
}




template <typename ValType>
typename RealOperationType(ValType) 
FESystem::Aeroelasticity::FrequencyDomainFlutterSolverBase<ValType>::getCurrentReducedFrequency()
{
    return this->current_k_ref;
}
    

/***************************************************************************************/
// Template instantiations for some generic classes

INSTANTIATE_CLASS_FOR_ONLY_COMPLEX_DATA_TYPES(FESystem::Aeroelasticity::FrequencyDomainFlutterSolverBase);


/***************************************************************************************/
