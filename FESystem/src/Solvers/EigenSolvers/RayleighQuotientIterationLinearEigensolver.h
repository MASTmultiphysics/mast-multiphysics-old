//
//  RayleighQuotientIterationLinearEigensolver.h
//  FESystem
//
//  Created by Manav Bhatia on 5/9/11.
//  Copyright 2011 . All rights reserved.
//

#ifndef __fesystem_rayleigh_quotient_iteration_linear_eigensolver_h__
#define __fesystem_rayleigh_quotient_iteration_linear_eigensolver_h__

// C++ includes
#include <memory>

// FESystem includes
#include "Solvers/EigenSolvers/LinearEigenSolverBase.h"

namespace FESystem
{
    namespace EigenSolvers
    {
        // Forward declerations
        namespace Numerics {template <typename ValType> class MatrixBase;}
        
        template <typename ValType> 
        class RayleighQuotientIterationLinearEigenSolver: public FESystem::EigenSolvers::LinearEigenSolverBase<ValType>
        {
        public:
            
            RayleighQuotientIterationLinearEigenSolver();
            
            virtual ~RayleighQuotientIterationLinearEigenSolver();
            
            void setShift(ValType v);
            
            virtual void solve();
            
        protected:
            
            void shiftAndInvertPowerIterations();
            
            ValType solver_shift;
            
            std::auto_ptr<FESystem::Numerics::MatrixBase<ValType> > shifted_mat; 
            
        };
        
    }
}



#endif // __fesystem_rayleigh_quotient_iteration_linear_eigensolver_h__


