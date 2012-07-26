//
//  UGFlutterSolver.cpp
//  FESystem
//
//  Created by Manav Bhatia on 7/24/12.
//  Copyright (c) 2012. All rights reserved.
//


// FESystem includes
#include "Disciplines/Aeroelasticity/UGFlutterSolver.h"


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
FESystem::Aeroelasticity::UGFlutterSolver<ValType>::solve()
{
    switch (this->current_call_back)
    {
        case CALCULATE_MATRICES: // matrices have been set, solve the eigenproblem
        {
            this->eigen_solver->solve();
            // get the eigenvalues and calculate the g and V
            FESystem::Numerics::VectorBase<ValType>& eig_vals = this->eigen_solver->getEigenValues();
            for (FESystemUInt eig_id=0; eig_id<eig_vals.getSize(); eig_id++)
            {
                g_eig = FESystem::Base::imag(eig_vals.getVal(eig_id));
                v_eig = FESystem::Base::real(eig_vals.getVal(eig_id));
                if (v_eig > 0.0)
                    v_eig = 1.0/sqrt(v_eig);
                else // VG sometimes gets into this issue of negative values for some iffy parameters
                    v_eig = FESystem::Base::getMachineMax<typename RealOperationType(ValType)>();
            }
        }
            break;
            
        default:
            break;
    }
}
