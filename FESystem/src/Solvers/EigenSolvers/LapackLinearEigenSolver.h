//
//  LapackLinearEigenSolver.h
//  FESystem
//
//  Created by Manav Bhatia on 5/3/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem_lapack_linear_eigensolver_h__
#define __fesystem_lapack_linear_eigensolver_h__

// FESystem includes
#include "Solvers/EigenSolvers/LinearEigenSolverBase.h"

namespace FESystem
{
    // Forward declerations
    namespace Numerics {template <typename ValType> class MatrixBase;}
    
    namespace Solvers
    {
        /*!
         *    This class uses the QR algorithm for solution of given eigenproblems. The algorithm 
         *    calculates all eigenvectors and eigenvalues for the system. 
         *
         *     \todo Improve the efficiency by application of shifts, and locking of converged eigenvalues.
         *     \todo The eigenvector and eigenvalues are currently only real numbers, and consideration of nonhermitian problems will need more careful implementation. 
         */
        
        template <typename ValType> 
        class LapackLinearEigenSolver: public FESystem::Solvers::LinearEigenSolverBase<ValType>
        {
        public:
            
            LapackLinearEigenSolver();
            
            virtual ~LapackLinearEigenSolver();
            
            //void setShift(ValType v);
            
            virtual void solve();
            
        protected:
            
        };
        
        template <> void FESystem::Solvers::LapackLinearEigenSolver<FESystemDouble>::solve();
        template <> void FESystem::Solvers::LapackLinearEigenSolver<FESystemFloat>::solve();
        template <> void FESystem::Solvers::LapackLinearEigenSolver<FESystemComplexDouble>::solve();
        template <> void FESystem::Solvers::LapackLinearEigenSolver<FESystemComplexFloat>::solve();
    }
}

#endif // __fesystem_lapack_linear_eigensolver_h__
