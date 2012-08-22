//
//  FlutterSolverBase.h
//  FESystem
//
//  Created by Manav Bhatia on 7/24/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem_flutter_solver_base_h__
#define __fesystem_flutter_solver_base_h__

// FESystem includes
#include "Base/FESystemTypes.h"


namespace FESystem
{
    namespace Numerics {template <typename ValType> class MatrixBase;}
    
    namespace Aeroelasticity
    {
        enum FlutterSolutionCallBack
        {
            WAITING_TO_BEGIN,
            UPDATE_AERODYNAMIC_MATRICES,
            UPDATE_STRUCTURAL_MATRICES,
            COMPLETED_FLUTTER_SOLUTION,
            FOUND_FLUTTER_ROOT,
            INVALID_CALL_BACK
        };
        
        /*!
         *   Base class for flutter solver
         */
        template <typename ValType>
        class FlutterSolverBase
        {
        public:
            /*!
             *   Constructor
             */
            FlutterSolverBase();
            
            virtual ~FlutterSolverBase();
            
            /*!
             *   Initializes the flutter solver for the given set of parameters
             */
            virtual void initialize();

            /*!
             *   Sets the structural matrices. The matrices are given as pointers. The stiffness matrix is always required, but the
             *   damping pointer could be NULL for a flutter solution
             */
            void setStructuralMatrices(const FESystem::Numerics::MatrixBase<typename RealOperationType(ValType)>*  stiff,
                                       const FESystem::Numerics::MatrixBase<typename RealOperationType(ValType)>*  mass,
                                       const FESystem::Numerics::MatrixBase<typename RealOperationType(ValType)>*  damp);
            
            /*!
             *    Solves for the flutter point
             */
            virtual FESystem::Aeroelasticity::FlutterSolutionCallBack solve()=0;
            
        protected:
            
            /*!
             *   Boolean to check if the solver is initialized
             */
            FESystemBoolean if_initialized;
            
            /*!
             *   reference chord length
             */
            typename RealOperationType(ValType) b_ref;

            /*!
             *   reference chord length
             */
            typename RealOperationType(ValType) fluid_rho;

            /*!
             *    Boolean to check if the structural matrices have been set
             */
            FESystemBoolean if_structural_matrices_available;
            
            /*!
             *    number of flutter roots found
             */
            FESystemUInt n_flutter_solutions_found;
            
            /*!
             *    current call back during flutter solution
             */
            FESystem::Aeroelasticity::FlutterSolutionCallBack current_call_back;
            
            /*!
             *   Pointer to the structural mass matrix
             */
            const FESystem::Numerics::MatrixBase<typename RealOperationType(ValType)>*  structural_mass_matrix;
            
            /*!
             *   Pointer to the structural damping matrix
             */
            const FESystem::Numerics::MatrixBase<typename RealOperationType(ValType)>*  structural_damping_matrix;
            
            /*!
             *   Pointer to the structural stiffness matrix
             */
            const FESystem::Numerics::MatrixBase<typename RealOperationType(ValType)>*  structural_stiffness_matrix;
            
        };
    }
}



#endif // __fesystem_flutter_solver_base_h__
