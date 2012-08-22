//
//  FlutterSolverBase.cpp
//  FESystem
//
//  Created by Manav Bhatia on 7/24/12.
//  Copyright (c) 2012. All rights reserved.
//

// FESystem includes
#include "Disciplines/Aeroelasticity/FlutterSolverBase.h"
#include "Base/FESystemExceptions.h"
#include "Base/macros.h"


template <typename ValType>
FESystem::Aeroelasticity::FlutterSolverBase<ValType>::FlutterSolverBase():
if_initialized(false),
b_ref(0.0),
fluid_rho(0.0),
if_structural_matrices_available(false),
n_flutter_solutions_found(0),
current_call_back(FESystem::Aeroelasticity::INVALID_CALL_BACK),
structural_mass_matrix(NULL),
structural_damping_matrix(NULL),
structural_stiffness_matrix(NULL)
{
    
}



template <typename ValType>
FESystem::Aeroelasticity::FlutterSolverBase<ValType>::~FlutterSolverBase()
{
    
}



template <typename ValType>
void
FESystem::Aeroelasticity::FlutterSolverBase<ValType>::initialize()
{
    
}
            

template <typename ValType>
void
FESystem::Aeroelasticity::FlutterSolverBase<ValType>::setStructuralMatrices(const FESystem::Numerics::MatrixBase<typename RealOperationType(ValType)>*  stiff,
                                                                            const FESystem::Numerics::MatrixBase<typename RealOperationType(ValType)>*  mass,
                                                                            const FESystem::Numerics::MatrixBase<typename RealOperationType(ValType)>*  damp)
{
    FESystemAssert0(!this->if_initialized, FESystem::Exception::InvalidState);
    FESystemAssert0(!this->if_structural_matrices_available, FESystem::Exception::InvalidState);
    
    this->structural_mass_matrix = stiff;
    this->structural_damping_matrix = damp;
    this->structural_stiffness_matrix = stiff;
}
            


/***************************************************************************************/
// Template instantiations for some generic classes

INSTANTIATE_CLASS_FOR_ONLY_COMPLEX_DATA_TYPES(FESystem::Aeroelasticity::FlutterSolverBase);


/***************************************************************************************/

