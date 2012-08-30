//
//  PseudoTimeSteppingLinearSolver.h
//  FESystem
//
//  Created by Manav Bhatia on 4/24/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem_pseudo_time_stepping_linear_solver_h__
#define __fesystem_pseudo_time_stepping_linear_solver_h__

// C++ includes
#include <memory>

// FESystem includes
#include "Solvers/LinearSolvers/LinearSolverBase.h"


namespace FESystem
{
    namespace LinearSolvers
    {
        /*!
         *    Solves the system \f$ A x = b \f$ using pseudo-time stepping where the equation is solved as 
         *    \f$ d x/ dt - A x + b = 0 \f$. It is important to note that this works only when the pseudo-time stepping
         *    system represents the transient analysis system. This is true, in general, for fluid-dynamic analysis 
         *    expressed in the conservative form. However, the structural system of equations, represented in the first-order 
         *    state-space form would need to be premultiplied with the mass matrix (with appropriate BCs applied) so
         *    that the system matrix is appropriately conditioned. 
         */
        template <typename ValType> 
        class PseudoTimeSteppingLinearSolver: public FESystem::LinearSolvers::LinearSolverBase<ValType>
        {
        public:
            /*!
             *   Constructor
             */
            PseudoTimeSteppingLinearSolver();
            
            virtual ~PseudoTimeSteppingLinearSolver();
            
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
             *    tolerance for convergence of the solver
             */
            FESystemDouble tolerance;
            
            
            /*!
             *    maximum iterations
             */
            FESystemUInt max_iters;
            
            
            /*!
             *   Temporary vector storage
             */
            std::auto_ptr<FESystem::Numerics::VectorBase<ValType> > residual_vec;
        };
    }
}





#endif  // __fesystem_pseudo_time_stepping_linear_solver_h__

