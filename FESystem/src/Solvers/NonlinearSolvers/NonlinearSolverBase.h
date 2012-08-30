//
//  NonlinearSolverBase.h
//  FESystem
//
//  Created by Manav Bhatia on 8/30/12.
//
//

#ifndef __FESystem__NonlinearSolverBase__
#define __FESystem__NonlinearSolverBase__

// C++ includes
#include <vector>

// FESystem includes
#include "Base/FESystemTypes.h"


namespace FESystem
{
    // Forward declerations
    namespace Numerics {template <typename ValType> class MatrixBase;}
    namespace Numerics {template <typename ValType> class VectorBase;}
    namespace Numerics {class SparsityPattern;}
    namespace Base {class DegreeOfFreedomMap;}
    
    namespace NonlinearSolvers
    {
        
        /*!
         *   enumeration defines the action required by the solver
         */
        enum NonlinearSolverCallBack
        {
            MAXIMUM_ITERATIONS_REACHED,
            SET_INITIAL_GUESS,
            EVALUATE_RESIDUAL,
            EVALUATE_JACOBIAN,
            EVALUATE_RESIDUAL_AND_JACOBIAN,
            SOLUTION_CONVERGED,
            WAITING_TO_START
        };
        
        
        /*!   this class provides an interface to a solver for the solution of
         *    transient system.
         */
        template <typename ValType>
        class NonlinearSolverBase
        {
        public:
            
            /*!
             *     constructor, with the order of the highest derivative defined as \p o.
             */
            NonlinearSolverBase();
            
            
            virtual ~NonlinearSolverBase();
            
            /*!
             *   Returns the state of initialization of this object
             */
            FESystemBoolean ifInitialized() const;
            
            /*!
             *   this method clears the data structures of this object. This should be called
             *   each time the used finishes using this object.
             */
            virtual void clear();
            
            /*!
             *   Sets the tolerance limits for the solver
             */
            void setConvergenceLimits(FESystemUInt max_iters, FESystemDouble tol);

            /*!
             *    Initialize solution with the given initial guess in \p vec
             */
            void setInitialGuess(FESystem::Numerics::VectorBase<ValType>& vec);
            
            /*!
             *   returns the current call back
             */
            FESystem::NonlinearSolvers::NonlinearSolverCallBack getCurrentCallBack() const;
            
            /*
             *   increments to the next time step
             */
            virtual FESystem::NonlinearSolvers::NonlinearSolverCallBack incrementSolution() = 0;
                        
            /*!
             *   Returns a reference to the current state vector
             */
            FESystem::Numerics::VectorBase<ValType>& getCurrentSolution();
            
            /*!
             *   Returns a reference to the current state velocity (or x dot) vector
             */
            FESystem::Numerics::VectorBase<ValType>& getResidualVector();
                        
            /*!
             *  Returns the current time iteration of the solver integration
             */
            virtual unsigned int getCurrentIterationNumber();
                        
        protected:
            
            /*!
             *   initializes the solver to a system with n_dofs unknowns and highest order of derivative \p o. This
             *   will create local dense matrices and vectors.
             */
            void initialize(FESystemUInt n);
            
            /*!
             *    initial time data
             */
            FESystemBoolean if_initialized;

            /*!
             *   convergence tolerance
             */
            typename RealOperationType(ValType) convergence_tolerance;

            /*!
             *   Maximum allowable iterations
             */
            FESystemUInt max_allowed_iterations;

                        
            /*!
             *    number of degrees of freedom for the spatial discretization
             */
            FESystemUInt n_dofs;
                                    
            /*!
             *    current iteration_number
             */
            FESystemUInt current_iteration_number;
            
            /*!
             *    latest action requested by the solver
             */
            FESystem::NonlinearSolvers::NonlinearSolverCallBack latest_call_back;
            
            /*!
             *    Stores the delta solution per iteration
             */
            FESystem::Numerics::VectorBase<ValType>* sol_increment_vec;

            /*!
             *    Stores the latest solution
             */
            FESystem::Numerics::VectorBase<ValType>* sol_vec;
            
            /*!
             *    Stores the latest residual vector
             */
            FESystem::Numerics::VectorBase<ValType>* residual;
        };
    }
}


#endif /* defined(__FESystem__NonlinearSolverBase__) */
