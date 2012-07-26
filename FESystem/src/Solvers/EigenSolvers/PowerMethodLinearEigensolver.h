//
//  PowerMethodLinearEigensolver.h
//  FESystem
//
//  Created by Manav Bhatia on 3/26/11.
//  Copyright 2011 . All rights reserved.
//

#ifndef __fesystem_power_method_eigensolver_h__
#define __fesystem_power_method_eigensolver_h__

// FESystem includes
#include "Solvers/EigenSolvers/LinearEigenSolverBase.h"

namespace FESystem
{
    namespace Solvers
    {
        // Forward declerations
        namespace Numerics {template <typename ValType> class MatrixBase;}
        
        template <typename ValType> 
        class PowerMethodLinearEigenSolver: 
        public FESystem::Solvers::LinearEigenSolverBase<ValType>
        {
        public:
            
            PowerMethodLinearEigenSolver();
            
            virtual ~PowerMethodLinearEigenSolver();
                        
            virtual void solve();
            
        protected:
                        
        };

    }
}


#endif // __fesystem_power_method_eigensolver_h__
