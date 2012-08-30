/*
 *  DirectLinearSolver.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 1/4/11.
 *  Copyright 2011 . All rights reserved.
 *
 */

#ifndef __fesystem_Lapack_direct_linear_solver_h__
#define __fesystem_Lapack_direct_linear_solver_h__

// C++ includes
#include <vector>

// FESystem includes
#include "Solvers/LinearSolvers/LinearSolverBase.h"


namespace FESystem
{
    namespace LinearSolvers
    {        
        /*!
         *   This is a template class that uses the Lapack solver to solve a linear system of equations
         *   \f$ A x = b \f$. The matrix \f$ A \f$ is given to the solver using the setSystemMatrix method,
         *   while the solve method is used to solve the equation. The solve method can be called a number of 
         *   times with multiple RHS. This class inherits from the base class LinearSolverBase<ValType>
         */
        template <typename ValType>
        class LapackLinearSolver: public FESystem::LinearSolvers::LinearSolverBase<ValType>
        {
        public:
            /*!
             *    Contructor
             */
            LapackLinearSolver();
            
            ~LapackLinearSolver();
            
            /*!
             *   Clears the data structures for used with another system matrix.
             */
            void clear();
            
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
             *   The matrix data is copied and stored in this vector to be passed on to Lapack. 
             */
            std::vector<ValType> system_matrix_work_copy;
            
            /*!
             *   Work array necessary for Lapack
             */
            std::vector<FESystemInt> ipiv;
            
        };
        
        DeclareException0(SystemMatrixWorkCopyNotInitialized, 
                          << "System Matrix Work Copy Not Initialized Before Usage\n");

        // template specialization
        template <>
        void FESystem::LinearSolvers::LapackLinearSolver<FESystemDouble>::setSystemMatrix(const FESystem::Numerics::MatrixBase<FESystemDouble> &mat);
        template <>
        void FESystem::LinearSolvers::LapackLinearSolver<FESystemDouble>::solve(const FESystem::Numerics::VectorBase<FESystemDouble> &rhs, FESystem::Numerics::VectorBase<FESystemDouble> &sol);
        template <>
        void FESystem::LinearSolvers::LapackLinearSolver<FESystemDouble>::solve(const FESystem::Numerics::MatrixBase<FESystemDouble> &rhs, FESystem::Numerics::MatrixBase<FESystemDouble> &sol);

        template <>
        void FESystem::LinearSolvers::LapackLinearSolver<FESystemFloat>::setSystemMatrix(const FESystem::Numerics::MatrixBase<FESystemFloat> &mat);
        template <>
        void FESystem::LinearSolvers::LapackLinearSolver<FESystemFloat>::solve(const FESystem::Numerics::VectorBase<FESystemFloat> &rhs, FESystem::Numerics::VectorBase<FESystemFloat> &sol);
        template <>
        void FESystem::LinearSolvers::LapackLinearSolver<FESystemFloat>::solve(const FESystem::Numerics::MatrixBase<FESystemFloat> &rhs, FESystem::Numerics::MatrixBase<FESystemFloat> &sol);
    }
}


#endif 
