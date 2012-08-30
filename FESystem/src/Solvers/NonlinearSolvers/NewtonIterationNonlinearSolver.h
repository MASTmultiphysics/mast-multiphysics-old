//
//  NewtonIterationNonlinearSolver.h
//  FESystem
//
//  Created by Manav Bhatia on 8/30/12.
//
//

#ifndef __FESystem__NewtonIterationNonlinearSolver__
#define __FESystem__NewtonIterationNonlinearSolver__

// FESystem includes
#include "Solvers/NonlinearSolvers/NonlinearSolverBase.h"


namespace FESystem
{
    
    // Forward declerations
    namespace LinearSolvers {template <typename  ValType> class LinearSolverBase;}
    

    namespace NonlinearSolvers
    {
        template <typename ValType>
        class NewtonIterationNonlinearSolver: public FESystem::NonlinearSolvers::NonlinearSolverBase<ValType>
        {
        public:
            /*!
             *   constructor
             */
            NewtonIterationNonlinearSolver();
            
            virtual ~NewtonIterationNonlinearSolver();
            
            /*!
             *    clears the data structures
             */
            virtual void clear();
            
            /*!
             *    initializes the solver and associated data structures with the Jacobian matrix specified
             */
            virtual void initialize(FESystem::Numerics::MatrixBase<ValType>& mat);
                        
            
            /*!
             *   Returns a reference to the Jacobian matrix
             */
            FESystem::Numerics::MatrixBase<ValType>& getJacobianMatrix();
            
            /*!
             *    provides the linear solver object for this transient solver. The boolean flag should be set to true if the system matrices
             *    are constant with respect to time.
             */
            void setLinearSolver(FESystem::LinearSolvers::LinearSolverBase<ValType>& solver);
            
            /*
             *   increments to the next time step
             */
            virtual FESystem::NonlinearSolvers::NonlinearSolverCallBack incrementSolution();
            
        protected:
            
            /*!
             *    Jacobian matrix used in this solution
             */
            FESystem::Numerics::MatrixBase<ValType>* jacobian;
            
            /*!
             *    Pointer to linear sovler
             */
            FESystem::LinearSolvers::LinearSolverBase<ValType>* linear_solver;
            
        };
    }
}


#endif /* defined(__FESystem__NewtonIterationNonlinearSolver__) */
