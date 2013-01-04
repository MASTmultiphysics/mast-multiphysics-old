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
            virtual void initialize(FESystem::Numerics::MatrixBase<ValType>& mat, FESystem::LinearSolvers::LinearSolverBase<ValType>& solver);
                        
            
            /*!
             *   Returns a reference to the Jacobian matrix
             */
            FESystem::Numerics::MatrixBase<ValType>& getJacobianMatrix();
                        
            /*
             *   increments to the next time step
             */
            virtual FESystem::NonlinearSolvers::NonlinearSolverCallBack incrementSolution();
            
        protected:
            
            /*!
             *   convergence tolerance
             */
            typename RealOperationType(ValType) base_res_l2, base_res_l2_slope, newton_step_res_l2;
            
            /*!
             *    Boolean to check if the backtracking has been done once for the current Newton step
             */
            FESystemBoolean if_backtracked;

            /*!
             *    Jacobian matrix used in this solution
             */
            FESystem::Numerics::MatrixBase<ValType>* jacobian;

            /*!
             *    temporary scratch vector
             */
            FESystem::Numerics::VectorBase<ValType>* temp_vec;

            /*!
             *    Pointer to linear sovler
             */
            FESystem::LinearSolvers::LinearSolverBase<ValType>* linear_solver;
            
        };
    }
}


#endif /* defined(__FESystem__NewtonIterationNonlinearSolver__) */
