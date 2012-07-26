//
//  ExplicitRungeKuttaTransientSolver.h
//  FESystem
//
//  Created by Manav Bhatia on 5/29/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem__explicit_runge_kutta_transient_solver_h__
#define __fesystem__explicit_runge_kutta_transient_solver_h__


// C++ includes
#include <vector>

// FESystem includes
#include "Solvers/TransientSolvers/TransientSolverBase.h"


namespace FESystem
{
    namespace Solvers
    {
        /*!
         *   provides the implementation of Newmark transient solver
         */
        template <typename ValType>
        class ExplicitRungeKuttaTransientSolver: public FESystem::Solvers::LinearTransientSolverBase<ValType>
        {
        public:
            /*!
             *   Constructor
             */
            ExplicitRungeKuttaTransientSolver();
            
            ~ExplicitRungeKuttaTransientSolver();
            
            /*!
             *   initializes the solver to a system with n_dofs unknowns and highest order of derivative \p o. This 
             *   will create local dense matrices and vectors. The constants of integration for the Newmark solver are defined 
             *   in the vector int_constants. The first element is for the 0th order state, followed by the 1st order time differential term,
             *   and so on.
             */
            void initialize(FESystemUInt o, FESystemUInt n_dofs, FESystemUInt n_rk_steps);
                        
            /*!
             *   this method clears the data structures of this object. This should be called 
             *   each time the used finishes using this object.
             */
            virtual void clear(); 
            
            /*
             *   increments to the next time step
             */
            virtual FESystem::Solvers::TransientSolverCallBack incrementTimeStep();
                        
//            /*!
//             *   Rewinds the time step to the beginning of the previous time iteration. This may be necessary for cases where the same time 
//             *   step needs to be evaluated with updated loads, for example in fluid structure interaction. 
//             */
//            virtual void rewindTimeStep();

        protected:
            
            /*!
             *    number of RK steps per time increment
             */
            FESystemUInt n_rk_steps_per_time_increment;
            
            /*!
             *    number of subiterations completed
             */
            FESystemUInt n_rk_steps_completed;
            
            /*!
             *   stores the previous time value during the sub-iterations
             */
            FESystemDouble previous_time;
            
            /*!
             *   Storage for state from previous time step
             */
            FESystem::Numerics::VectorBase<ValType>* previous_state;
            
            /*!
             *   Storage for velocity from previous time step
             */
            FESystem::Numerics::VectorBase<ValType>* new_state_estimate;

            /*!
             *    coefficients of each sub-iterate towards the final increment
             */
            std::vector<FESystemDouble> sub_step_coefficients_for_final_step;

            /*!
             *    sub-step coefficients for evaluating next sub-iterate
             */
            std::vector<FESystemDouble> sub_step_iterate_coefficients;

        };
    }
}




#endif  // __fesystem__explicit_runge_kutta_transient_solver_h__
