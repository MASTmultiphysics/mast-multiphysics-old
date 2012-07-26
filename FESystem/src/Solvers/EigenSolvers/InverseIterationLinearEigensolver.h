//
//  InverseIterationLinearEigensolver.h
//  FESystem
//
//  Created by Manav Bhatia on 3/27/11.
//  Copyright 2011 . All rights reserved.
//

#ifndef __fesystem_inverse_iteration_linear_eigensolver_h__
#define __fesystem_inverse_iteration_linear_eigensolver_h__

// C++ includes
#include <memory>

// FESystem includes
#include "Solvers/EigenSolvers/LinearEigenSolverBase.h"

namespace FESystem
{
    // Forward declerations
    namespace Numerics {template <typename ValType> class MatrixBase;}

    namespace Solvers
    {
        template <typename ValType> 
        class InverseIterationLinearEigenSolver: 
        public FESystem::Solvers::LinearEigenSolverBase<ValType>
        {
        public:
            
            InverseIterationLinearEigenSolver();
            
            virtual ~InverseIterationLinearEigenSolver();
            
            void setShift(ValType v);
            
            virtual void solve();
            
        protected:
            
            void shiftAndInvertPowerIterations(FESystem::Numerics::MatrixBase<ValType>& mat);
            
            ValType solver_shift;
            
            /*!
             *  \brief storage for the shifted matrix
             */
            std::auto_ptr<FESystem::Numerics::MatrixBase<ValType> > shifted_mat; 
            
        };
        
    }
}



#endif // __fesystem_inverse_iteration_linear_eigensolver_h__

