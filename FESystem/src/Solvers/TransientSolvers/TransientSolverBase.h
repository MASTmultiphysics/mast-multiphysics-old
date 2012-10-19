//
//  TransientSolverBase.h
//  FESystem
//
//  Created by Manav Bhatia on 4/25/12.
//  Copyright (c) 2012. All rights reserved.
//


#ifndef __fesystem_transient_solver_base_h__
#define __fesystem_transient_solver_base_h__

// C++ includes
#include <vector>

// FESystem includes
#include "Base/FESystemTypes.h"


namespace FESystem 
{
    // Forward declerations
    namespace LinearSolvers {template <typename  ValType> class LinearSolverBase;}
    namespace Numerics {template <typename ValType> class MatrixBase;}
    namespace Numerics {template <typename ValType> class VectorBase;}
    namespace Numerics {class SparsityPattern;}
    namespace Base {class DegreeOfFreedomMap;}
    
    namespace TransientSolvers
    {
        
        /*!
         *   enumeration defines the action required by the solver
         */
        enum TransientSolverCallBack
        {
            EVALUATE_X_DOT,
            EVALUATE_X_DOT_JACOBIAN,
            EVALUATE_X_DOT_AND_X_DOT_JACOBIAN,
            TIME_STEP_CONVERGED,
            WAITING_TO_START
        };
        
        
        /*!   this class provides an interface to a solver for the solution of 
         *    transient system. 
         */
        template <typename ValType>
        class TransientSolverBase 
        {
        public:
            
            /*!
             *     constructor, with the order of the highest derivative defined as \p o. 
             */ 
            TransientSolverBase();
            
            
            virtual ~TransientSolverBase();

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
             *    Initialize solution with the given initial time and step size and the initial conditions in \p vec
             */
            void setInitialTimeData(typename RealOperationType(ValType) t0, typename RealOperationType(ValType) dt, FESystem::Numerics::VectorBase<ValType>& vec);

            /*
             *   increments to the next time step
             */
            virtual FESystem::TransientSolvers::TransientSolverCallBack incrementTimeStep() = 0;
            
//            /*!
//             *   Rewinds the time step to the beginning of the previous time iteration. This may be necessary for cases where the same time 
//             *   step needs to be evaluated with updated loads, for example in fluid structure interaction. 
//             */
//            virtual void rewindTimeStep() = 0;

            /*!
             *   Returns a reference to the current state vector
             */
            FESystem::Numerics::VectorBase<ValType>& getCurrentStateVector();

            /*!
             *   Returns a reference to the current state velocity (or x dot) vector 
             */
            FESystem::Numerics::VectorBase<ValType>& getCurrentStateVelocityVector();

            /*!
             *   Returns a reference to the previous state vector
             */
            FESystem::Numerics::VectorBase<ValType>& getPreviousStateVector();
            
            /*!
             *   Returns a reference to the previous state velocity (or x dot) vector
             */
            FESystem::Numerics::VectorBase<ValType>& getPreviousStateVelocityVector();

            /*
             *   Returns the current time of the solver integration
             */
            virtual typename RealOperationType(ValType) getCurrentTime();
            
            /*!
             *   Returns the current time step of the solver integration
             */
            virtual typename RealOperationType(ValType) getCurrentStepSize();
            
            /*!
             *  Returns the current time iteration of the solver integration
             */
            virtual FESystemUInt getCurrentIterationNumber();
            
            /*!
             *   Initializes the state vector to the appropriate dimension for this transient solver. The
             *   size will be order times the number of degrees of freedom in the system
             */
            virtual void initializeStateVector(FESystem::Numerics::VectorBase<ValType>& vec);
            
            /*!
             *   Sets the options of whether or not the Mass matrix (coefficient matrix of highest derivative term on LHS) is an Identity matrix. This is true by default.
             *   If that is not the case, then the user needs to provide a pointer to the mass matrix.
             */
            void setMassMatrix(FESystemBoolean if_identity, FESystem::Numerics::MatrixBase<ValType>* mass_mat_ptr=NULL);
            
            /*!
             *   Sets the active terms for the Jacobian matrix. These correspond to the RHS of equation for the higherst time derivative
             */
            void setActiveJacobianTerm(std::vector<FESystemBoolean>& active_terms);

            /*!
             *   Resizes the matrix to the appropriate dimensions for this system. This is meant to be a utility function for dense matrices, and should not be used for
             *   SparseMatrix objects. For SparseMatrix objects, use the sparsity pattern update utility method. 
             */
            virtual void resizeMatrixToJacobianTemplate(FESystem::Numerics::MatrixBase<ValType>& state_jac);

            
            /*!
             *   Initializes the matrix to the Jacobian template for the given order of the system. 
             */
            virtual void initializeMatrixToJacobianTemplate(FESystem::Numerics::MatrixBase<ValType>& state_jac);

            /*!
             *   Initializes the sparsity pattern for this transient solver. The sparsity_pattern for this transient system is returned in \p system_pattern. 
             */
            virtual void initializeMatrixSparsityPatterForSystem(const FESystem::Numerics::SparsityPattern& spatial_sparsity_pattern,  FESystem::Numerics::SparsityPattern& system_pattern) const;

            /*!
             *   Update vector values from the dof vector for derivative order \p o into the state vector for the solver
             */
            virtual void updateVectorValuesForDerivativeOrder(FESystemUInt o, const FESystem::Numerics::VectorBase<ValType>& dof_vec,
                                                              FESystem::Numerics::VectorBase<ValType>& state_vec);

            /*!
             *   Update Jacobian values from the matrix that defines the quantity d x^p / d x^q . Where x^p is d^p x / dt^p and x^q is d^q x/ dt^q. 
             */
            virtual void updateJacobianValuesForDerivativeOrder(FESystemUInt p, FESystemUInt q,
                                                                const FESystem::Numerics::MatrixBase<ValType>& dof_jac, 
                                                                FESystem::Numerics::MatrixBase<ValType>& state_jac);
            
            /*!
             *   Extracts vector values from the state vector for derivative order \p o into the dof vector
             */
            virtual void extractVectorValuesForDerivativeOrder(FESystemUInt o, const FESystem::Numerics::VectorBase<ValType>& state_vec, 
                                                               FESystem::Numerics::VectorBase<ValType>& dof_vec);
            
            /*!
             *   This is a utility method to copy the derivative values of order 1 to (o-1), where o is the highest derivative of this 
             *   solver, from the given state vector \p state to the given state velocity vector \p velocity. This is essentially the same as
             *   \$[ velocity(0:(o-1)*n_dofs-1) = state(n_dofs:o*n_dofs-1) \$]
             */
            virtual void copyDerivativeValuesFromStateToVelocityVector(const FESystem::Numerics::VectorBase<ValType>& state, FESystem::Numerics::VectorBase<ValType>& velocity) const;
            
            /*!
             *    provides the linear solver object for this transient solver. The boolean flag should be set to true if the system matrices
             *    are constant with respect to time.
             */
            void setLinearSolver(FESystem::LinearSolvers::LinearSolverBase<ValType>& solver, FESystemBoolean if_constant_matrices);

        protected:

            /*!
             *   initializes the solver to a system with n_dofs unknowns and highest order of derivative \p o. This 
             *   will create local dense matrices and vectors. 
             */
            void initialize(FESystemUInt o, FESystemUInt n_dofs);            

            /*!
             *    initial time data
             */
            FESystemBoolean if_initialized;
                        
            /*!
             *    if the initial time data has been initialized
             */
            FESystemBoolean if_initialized_initial_time_data;
            
            /*!
             *    boolean about whether the mass matrix is assumed to be an Identity matrix. The mass matrix is the coefficient of the highest time derivative term on LHS
             */
            FESystemBoolean if_identity_mass_matrix;
            
            /*!
             *   Mass matrix pointer if it is not identity
             */
            FESystem::Numerics::MatrixBase<ValType>* mass_matrix;
            
            /*!
             *    order of highest time derivative
             */
            FESystemUInt order;

            /*!
             *    Boolean about whether or not the respective derivative term in the Jacobian exists. This is set to true by default, and the user can 
             *    change it if needed.
             */
            std::vector<FESystemBoolean> active_jacobian_terms;
            
            /*!
             *    number of degrees of freedom for the spatial discretization
             */
            FESystemUInt n_dofs;

            /*!
             *    Initial time value
             */
            typename RealOperationType(ValType) initial_time;

            /*!
             *    Current time value
             */
            typename RealOperationType(ValType) current_time;

            /*!
             *    current time step
             */
            typename RealOperationType(ValType) current_time_step;

            /*!
             *    current iteration_number
             */
            FESystemUInt current_iteration_number;
            
            /*!
             *    latest action requested by the solver
             */
            FESystem::TransientSolvers::TransientSolverCallBack latest_call_back;

            /*!
             *    Stores the solution at the last solved time step, used to store the initial conditions at the beginning of solution. 
             */
            FESystem::Numerics::VectorBase<ValType> *current_state, *current_velocity, *previous_state, *previous_velocity;
            
            /*!
             *    Pointer to linear sovler
             */
            FESystem::LinearSolvers::LinearSolverBase<ValType>* linear_solver;
            
            
            /*!
             *    In case the system matrices are constant with respect to time, the linear solver will be initialized only once for the first time step.
             */
            FESystemBoolean if_constant_system_matrices;

        };
    }
}


#endif // __fesystem_transient_solver_base_h__


