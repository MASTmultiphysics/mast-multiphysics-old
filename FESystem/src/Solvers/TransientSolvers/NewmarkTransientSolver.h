//
//  NewmarkTransientSolver.h
//  FESystem
//
//  Created by Manav Bhatia on 4/27/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem_newmark_transient_solver_h__
#define __fesystem_newmark_transient_solver_h__

// C++ includes
#include <vector>

// FESystem includes
#include "Solvers/TransientSolvers/TransientSolverBase.h"


namespace FESystem
{
    namespace TransientSolvers
    {
        /*!
         *   provides the implementation of Newmark transient solver
         */
        template <typename ValType>
        class NewmarkTransientSolver: public FESystem::TransientSolvers::TransientSolverBase<ValType>
        {
        public:
            /*!
             *   Constructor
             */
            NewmarkTransientSolver();
            
            virtual ~NewmarkTransientSolver();
            
            /*!
             *   initializes the solver to a system with n_dofs unknowns and highest order of derivative \p o. This 
             *   will create local dense matrices and vectors. The constants of integration for the Newmark solver are defined 
             *   in the vector int_constants. The first element is for the 0th order state, followed by the 1st order time differential term,
             *   and so on.
             */
            void initialize(FESystemUInt o, FESystemUInt n_dofs, const std::vector<typename RealOperationType(ValType)>& int_constants);
                        
            /*!
             *   this method clears the data structures of this object. This should be called 
             *   each time the used finishes using this object.
             */
            virtual void clear(); 

            /*
             *   increments to the next time step
             */
            virtual FESystem::TransientSolvers::TransientSolverCallBack incrementTimeStep();
            
//            /*!
//             *   Rewinds the time step to the beginning of the previous time iteration. This may be necessary for cases where the same time 
//             *   step needs to be evaluated with updated loads, for example in fluid structure interaction. 
//             */
//            virtual void rewindTimeStep();
            
            /*
             *   sets the Jacobian matrix
             */
            virtual void setJacobianMatrix(FESystem::Numerics::MatrixBase<ValType>& jac);

            /*
             *   returns the Jacobian matrix
             */
            virtual FESystem::Numerics::MatrixBase<ValType>& getCurrentJacobianMatrix();

            
            /*!
             *
             */
            virtual void initializeMatrixSparsityPatterForSystem(const FESystem::Numerics::SparsityPattern& spatial_sparsity_pattern,  FESystem::Numerics::SparsityPattern& system_pattern) const;
            
            void setConvergenceTolerance(FESystemDouble tol, FESystemUInt max_itrs);
            
        protected:
 
            void evaluateResidual(const FESystem::Numerics::VectorBase<ValType>& prev_state, const FESystem::Numerics::VectorBase<ValType>&  prev_velocity,
                                  const FESystem::Numerics::VectorBase<ValType>&  curr_state, const FESystem::Numerics::VectorBase<ValType>&  curr_velocity,
                                  FESystem::Numerics::VectorBase<ValType>&  res);

            /*!
             *   convergence tolerance
             */
            typename RealOperationType(ValType) convergence_tolerance, base_res_l2, base_res_l2_slope, newton_step_res_l2;
            
            /*!
             *    Boolean to check if the backtracking has been done once for the current Newton step
             */
            FESystemBoolean if_backtracked;
            
            
            /*!
             *   Maximum allowable iterations
             */
            FESystemUInt  nonlinear_iteration_number, max_nonlinear_iterations;
            
            /*!
             *   constants for integration rule definition, starting with the lowest order and going towards the (o-1) order
             */
            std::vector<typename RealOperationType(ValType)> integration_constants;
            
            /*!
             *   Jacobian matrix for the solver
             */
            FESystem::Numerics::MatrixBase<ValType>* jacobian;
            
            /*!
             *   temporary storage vector
             */
            FESystem::Numerics::VectorBase<ValType> *residual, *newton_step, *temp_vec, *temp_vec2;

        };
    }
}


#endif // __fesystem_newmark_transient_solver_h__
