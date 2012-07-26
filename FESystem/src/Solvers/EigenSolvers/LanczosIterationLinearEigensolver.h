//
//  LanczosIterationLinearEigensolver.h
//  FESystem
//
//  Created by Manav Bhatia on 7/14/11.
//  Copyright 2011 . All rights reserved.
//

#ifndef __fesystem_lanczos_iteration_linear_eigensolver_h__
#define __fesystem_lanczos_iteration_linear_eigensolver_h__

// C++ includes
#include <memory>

// FESystem includes
#include "Solvers/EigenSolvers/LinearEigenSolverBase.h"

namespace FESystem
{
    namespace Solvers
    {
        // Forward declerations
        namespace Numerics {template <typename ValType> class MatrixBase;}
        
        /*!
         *    This class uses the Lanczos method for solution of given eigenproblems. 
         *
         */
        
        template <typename ValType> 
        class LanczosIterationLinearEigenSolver: 
        public FESystem::Solvers::LinearEigenSolverBase<ValType>
        {
        public:
            
            LanczosIterationLinearEigenSolver();
            
            virtual ~LanczosIterationLinearEigenSolver();
            
            //void setShift(ValType v);
                        
            virtual void solve();
            
        protected:
            
            /*!
             *   This is to create the matrix \p{eig_mat} and vector \p{eig_vec} for use during
             *   iterations
             */
            void initializeMatrices();

            
            /*!
             * \brief number of Krylov basis vector at the end of Lanczos iterations
             */
            FESystemUInt n_krylov_basis;

            /*!
             *  \brief storage for the Krylov basis as they are computed in the iterations
             */
            std::auto_ptr<FESystem::Numerics::MatrixBase<ValType> > krylov_basis_mat; 

            /*!
             *  \brief storage for the tridiagonal / Hessenberg matrix as it is updated in the iterations
             */
            std::auto_ptr<FESystem::Numerics::MatrixBase<ValType> > hessenberg_mat; 

        };
        
    }
}



#endif // __fesystem_lanczos_iteration_linear_eigensolver_h__

