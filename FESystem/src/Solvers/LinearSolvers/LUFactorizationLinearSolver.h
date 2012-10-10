//
//  LUFactorizationLinearSolver.h
//  FESystem
//
//  Created by Manav Bhatia on 6/17/12.
//  Copyright (c) 2012. All rights reserved.
//


#ifndef __fesystem_lu_factorization_linear_solver_h__
#define __fesystem_lu_factorization_linear_solver_h__


// C++ includes
#include <memory>

// FESystem includes
#include "Solvers/LinearSolvers/LinearSolverBase.h"


namespace FESystem
{
    // Forward declerations
    namespace FactorizationSolvers {template <typename ValType> class LUFactorization;}
    namespace FactorizationSolvers {template <typename ValType> class TriangularBacksubstitution;}
    

    namespace LinearSolvers
    {
        /*!
         *    Provides LU factorization based linear solver for system of equations \f$ A x= b \f$. The method first 
         *    creates a factorization \f$ A = LU \f$ where \f$ L \f$ is a lower triangular matrix, and then solves for each 
         *    RHS using two steps: (1) \f$ b1 = L^{-1} \f$ and (2) \f$ x = U^{-1} b1 \f$, where both steps use triangular 
         *    back substitution. 
         */
        template <typename ValType> 
        class LUFactorizationLinearSolver: public FESystem::LinearSolvers::LinearSolverBase<ValType>
        {
        public:
            /*!
             *   Constructor
             */
            LUFactorizationLinearSolver();
            
            virtual ~LUFactorizationLinearSolver();
            
            /*!
             *   Clears the data structures for used with another system matrix.
             */
            virtual void clear();
            
            /*!
             *   Sets the system matrix \f$ A \f$ in the system of equation through parameter \p mat and
             *   initializes the necessary data structures for a solution. If the solver is already associated with a 
             *   different matrix, the clear() method must be called before reassigning the system matrix. 
             */
            virtual void setSystemMatrix(const FESystem::Numerics::MatrixBase<ValType>& mat);
            
            /*!
             *   The setSystemMatrix must be called before this method. The method solves \f$ A x = b \f$ where 
             *   \f$ b \f$ is given as the parameter \p rhs and \f$ x \f$ is given as the parameter \p sol. 
             */
            virtual void solve(const FESystem::Numerics::VectorBase<ValType>& rhs,
                               FESystem::Numerics::VectorBase<ValType>& sol);
            
            /*!
             *   The setSystemMatrix must be called before this method. The method solves \f$ A X = B \f$ where 
             *   \f$ B \f$ is a matrix given as the parameter \p rhs and \f$ x \f$ is a matrix given as the parameter \p sol. 
             */
            virtual void solve(const FESystem::Numerics::MatrixBase<ValType>& rhs,
                               FESystem::Numerics::MatrixBase<ValType>& sol);
            
        protected:
            
            /*!
             *   Temporary vector storage
             */
            std::auto_ptr<FESystem::Numerics::VectorBase<ValType> > vec;
            
            
            /*!
             *   Temporary matrix storage
             */
            std::auto_ptr<FESystem::Numerics::MatrixBase<ValType> > mat;
            
            
            /*!
             *    A smart pointer to the QR factorization object
             */
            std::auto_ptr<FESystem::FactorizationSolvers::LUFactorization<ValType> >  lu_factorization;
            
            /*!
             *    A smart pointer to the triangular backsubsitution object
             */
            std::auto_ptr<FESystem::FactorizationSolvers::TriangularBacksubstitution<ValType> >  l_triangular_backsubstitute;
            std::auto_ptr<FESystem::FactorizationSolvers::TriangularBacksubstitution<ValType> >  u_triangular_backsubstitute;
            
        };
    }
}





#endif  //__fesystem_lu_factorization_linear_solver_h__

