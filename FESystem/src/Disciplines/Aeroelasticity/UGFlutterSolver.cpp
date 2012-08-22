//
//  UGFlutterSolver.cpp
//  FESystem
//
//  Created by Manav Bhatia on 7/24/12.
//  Copyright (c) 2012. All rights reserved.
//


// FESystem includes
#include "Disciplines/Aeroelasticity/UGFlutterSolver.h"
#include "Numerics/MatrixBase.h"
#include "Numerics/VectorBase.h"
#include "Solvers/EigenSolvers/LinearEigenSolverBase.h"
#include "Base/FESystemExceptions.h"
#include "Base/FESystemTypes.h"
#include "Base/macros.h"


template <typename ValType>
FESystem::Aeroelasticity::UGFlutterSolver<ValType>::UGFlutterSolver():
FESystem::Aeroelasticity::FrequencyDomainFlutterSolverBase<ValType>()
{
    
}


template <typename ValType>
FESystem::Aeroelasticity::UGFlutterSolver<ValType>::~UGFlutterSolver()
{
    
}


template <typename ValType>
void
FESystem::Aeroelasticity::UGFlutterSolver<ValType>::setReducedFrequency(std::vector<typename RealOperationType(ValType)>& k_vals)
{
    
}



template <typename ValType>
void
FESystem::Aeroelasticity::UGFlutterSolver<ValType>::setReducedFrequency(typename RealOperationType(ValType) k_low, typename RealOperationType(ValType) k_up, FESystemUInt n_divs)
{
    
}



template <typename ValType>
void
FESystem::Aeroelasticity::UGFlutterSolver<ValType>::setAerodynamicMatrix(const FESystem::Numerics::MatrixBase<ValType>& q_mat)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    FESystemAssert0(this->current_call_back == FESystem::Aeroelasticity::UPDATE_AERODYNAMIC_MATRICES, FESystem::Exception::InvalidState);
    
    this->aerodynamic_matrix = &q_mat;
}


template <typename ValType>
FESystem::Aeroelasticity::FlutterSolutionCallBack
FESystem::Aeroelasticity::UGFlutterSolver<ValType>::solve()
{
    switch (this->current_call_back)
    {
        case FESystem::Aeroelasticity::WAITING_TO_BEGIN:
        {
            this->current_call_back = FESystem::Aeroelasticity::UPDATE_STRUCTURAL_MATRICES;
            return FESystem::Aeroelasticity::UPDATE_STRUCTURAL_MATRICES;
        }
            break;

        case FESystem::Aeroelasticity::UPDATE_STRUCTURAL_MATRICES:
        {
            // set the reduced frequency to the lowest value
            if (this->if_given_k_ref_vals)
                this->current_k_ref = *(this->k_ref_vals.rbegin());
            else
                this->current_k_ref = this->k_ref_range.second;
            
            this->current_call_back = FESystem::Aeroelasticity::UPDATE_AERODYNAMIC_MATRICES;
            return FESystem::Aeroelasticity::UPDATE_AERODYNAMIC_MATRICES;
        }
            break;

        case UPDATE_AERODYNAMIC_MATRICES: // matrices have been set, solve the eigenproblem
        {
            this->eigenSolution();  // prepare matrices for eigensolution
            
            // get the eigenvalues and calculate the g and V
            typename RealOperationType(ValType) g_eig, v_eig;
            const FESystem::Numerics::VectorBase<ValType>& eig_vals = this->eigen_solver->getEigenValues();
            for (FESystemUInt eig_id=0; eig_id<eig_vals.getSize(); eig_id++)
            {
                g_eig = FESystem::Base::imag<ValType, typename RealOperationType(ValType)>(eig_vals.getVal(eig_id));
                v_eig = FESystem::Base::real<ValType, typename RealOperationType(ValType)>(eig_vals.getVal(eig_id));
                if (v_eig > 0.0)
                    v_eig = 1.0/sqrt(v_eig);
                else // VG sometimes gets into this issue of negative values for some combination of system parameter
                    v_eig = FESystem::Base::getMachineMax<typename RealOperationType(ValType)>();
            }
            
            // set the next set of reduced frequency
            FESystemBoolean if_terminate = true;
            // check if the solution needs to terminate
            if (if_terminate)
            {
                // perform any necessary operations
                this->current_call_back = FESystem::Aeroelasticity::COMPLETED_FLUTTER_SOLUTION;
            }
            else
            {
                // move to the next k_ref value
                this->current_call_back = FESystem::Aeroelasticity::UPDATE_AERODYNAMIC_MATRICES;
            }
            
        }
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
            break;
    }
}



template <typename ValType>
void
FESystem::Aeroelasticity::UGFlutterSolver<ValType>::eigenSolution()
{
    // B = K
    if (this->current_iter_num == 0)
        this->B_mat->copyRealMatrix(*(this->structural_stiffness_matrix));
    
    // A = (k/b)^2 M + rho/2 * Q
    this->A_mat->copyRealMatrix(*(this->structural_mass_matrix));
    this->A_mat->scale(pow(this->current_k_ref/this->b_ref,2));
    this->A_mat->add(0.5*this->fluid_rho, *(this->aerodynamic_matrix));
  
    // set the matrix and solve
    this->eigen_solver->setMatrix(this->A_mat, this->B_mat);
    this->eigen_solver->solve();
    
}



/***************************************************************************************/
// Template instantiations for some generic classes

INSTANTIATE_CLASS_FOR_ONLY_COMPLEX_DATA_TYPES(FESystem::Aeroelasticity::UGFlutterSolver);


/***************************************************************************************/

