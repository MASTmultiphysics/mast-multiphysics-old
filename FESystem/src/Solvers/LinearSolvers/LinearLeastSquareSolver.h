/*
 *  LinearLeastSquareSolver.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 1/8/11.
 *  Copyright 2011 . All rights reserved.
 *
 */

#ifndef __fesystem_linear_least_square_solver_h__
#define __fesystem_linear_least_square_solver_h__

// C++ includes
#include <memory>

// FESystem includes
#include "Base/FESystemExceptions.h"
#include "Solvers/LinearSolvers/LinearSolverBase.h"

namespace FESystem 
{
    // Forward declerations
    namespace Numerics {template <typename ValType> class VectorBase;}
    namespace Numerics {template <typename ValType> class MatrixBase;}

    namespace Solvers
    {
        
        /*!
         *   This is a template class to solve a linear least square system of equations
         *   \f$ A x = b \f$. The matrix \f$ A \f$ is not rectangular, and the following operations are used to solve
         *   the system of equations: (1) \f$ A1 = A^T A \f$, (2) \f$ b1 = A^T b \f$, and (3) solve \f$ A1 x = b1 \f$. 
         *
         *   The least square solver uses a linear solver that must be provided using the setLinearSolver method. 
         * 
         *   The matrix \f$ A \f$ is given to the solver using the setSystemMatrix method,
         *   while the solve method is used to solve the equation. The solve method can be called a number of 
         *   times with multiple RHS. This class inherits from the base class LinearSolverBase<ValType>
         */
        template <typename ValType>
        class LinearLeastSquareSolver: public LinearSolverBase<ValType>
        {
        public:
            /*!
             *   Constructor
             */
            LinearLeastSquareSolver();
            
            ~LinearLeastSquareSolver();
            
            /*!
             *   clears the data structures for use with a second matrix
             */
            void clear();
            
            /*!
             *   Sets the linear solver to be used by the equation. 
             */
            virtual void setLinearSolver(FESystem::Solvers::LinearSolverBase<ValType>& solver);
            
            /*!
             *   Returns a reference to the linear solver object
             */
            FESystem::Solvers::LinearSolverBase<ValType>& getLinearSolver();

            /*!
             *   Sets the system matrix \f$ A \f$ in the system of equation through parameter \p mat and
             *   initializes the necessary data structures for a solution. If the solver is already associated with a 
             *   different matrix, the clear() method must be called before reassigning the system matrix. 
             */
            void setSystemMatrix(const FESystem::Numerics::MatrixBase<ValType>& mat);
            
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
             *    A pointer to the linear solver object, which is set using the setLinearSolver method. 
             */
            FESystem::Solvers::LinearSolverBase<ValType>* linear_solver;
            
            /*!
             *    A smart pointer to the least square system matrix \f$ A1 = A^T A \f$.
             */
            std::auto_ptr<FESystem::Numerics::MatrixBase<ValType> > least_square_system_matrix;
        };
        
        DeclareException0(LinearSolverNotInitialized, 
                          << "Linear Solver Not Initialized Before Usage\n");
    }
}




#endif // __fesystem_linear_least_square_solver_h__
