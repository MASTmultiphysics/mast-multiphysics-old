/*
 *  LinearSolverBase.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 1/4/11.
 *  Copyright 2011 . All rights reserved.
 *
 */

#ifndef __fesystem_linear_solver_base_h__
#define __fesystem_linear_solver_base_h__

// FESystem includes
#include "Base/FESystemTypes.h"
#include "Base/FESystemExceptions.h"


namespace FESystem
{
    // Forward declerations
    namespace Numerics {template <typename ValType> class VectorBase;}
    namespace Numerics {template <typename ValType> class MatrixBase;}
    
    namespace LinearSolvers
    {
        /*!
         *    Provides a base class to define the interface for linear solver of system of equations \f$Ax= b\f$.
         */
        template <typename ValType> 
        class LinearSolverBase
        {
        public:
            /*!
             *   Constructor
             */
            LinearSolverBase();
            
            virtual ~LinearSolverBase();

            /*!
             *   clears the data structures for use with a second matrix
             */
            virtual void clear();
            
            /*!
             *   Sets the system matrix \f$ A \f$ in the system of equation through parameter \p mat and
             *   initializes the necessary data structures for a solution. If the solver is already associated with a 
             *   different matrix, the clear() method must be called before reassigning the system matrix. 
             */
            virtual void setSystemMatrix(const FESystem::Numerics::MatrixBase<ValType>& mat);

            /*!
             *    Returns a constant reference to the system matrix. The matrix should be set using setSystemMatrix 
             *    before this method is called. 
             */
            const FESystem::Numerics::MatrixBase<ValType>& getSystemMatrix() const;

            /*!
             *   The setSystemMatrix must be called before this method. The method solves \f$ A x = b \f$ where 
             *   \f$ b \f$ is given as the parameter \p rhs and \f$ x \f$ is given as the parameter \p sol. 
             */
            virtual void solve(const FESystem::Numerics::VectorBase<ValType>& rhs,
                               FESystem::Numerics::VectorBase<ValType>& sol) = 0;

            /*!
             *   The setSystemMatrix must be called before this method. The method solves \f$ A X = B \f$ where 
             *   \f$ B \f$ is a matrix given as the parameter \p rhs and \f$ X \f$ is a matrix given as the parameter \p sol. 
             */
            virtual void solve(const FESystem::Numerics::MatrixBase<ValType>& rhs,
                               FESystem::Numerics::MatrixBase<ValType>& sol) = 0;

        protected:
            
            /*!
             *   stores the state of the solver, whether or not it is initialized
             */
            FESystemBoolean if_initialized;
                    
            /*!
             *    A pointer to the system matrix. 
             */
            const FESystem::Numerics::MatrixBase<ValType>* system_matrix;
        };
        
        DeclareException0(SystemMatrixNotInitialized, 
                          << "System Matrix Not Initialized Before Usage\n");

        DeclareException0(SystemMatrixNotClearedBeforeAssignment, 
                          << "System Matrix Not Cleared Before Assignment\n");

    }
}


#endif 

