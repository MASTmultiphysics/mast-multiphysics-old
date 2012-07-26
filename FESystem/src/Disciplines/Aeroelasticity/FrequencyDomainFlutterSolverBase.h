//
//  FrequencyDomainFlutterSolverBase.h
//  FESystem
//
//  Created by Manav Bhatia on 7/24/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem_frequency_domain_flutter_solver_base_h__
#define __fesystem_frequency_domain_flutter_solver_base_h__


// FESystem includes
#include "Disciplines/Aeroelasticity/FlutterSolverBase.h"
#include "Base/FESystemTypes.h"


namespace FESystem
{    
    namespace Solvers {template <typename ValType> class LinearEigenSolverBase;}
    namespace Numerics {template <typename ValType> class MatrixBase;}
    namespace Numerics {template <typename ValType> class VectorBase;}
    
    namespace Aeroelasticity
    {
        /*!
         *   Frequency domani flutter solver base class. 
         */
        template <typename ValType>
        class FrequencyDomainFlutterSolverBase: public FESystem::Aeroelasticity::FlutterSolverBase<ValType>
        {
        public:
            /*!
             *    Constructor
             */
            FrequencyDomainFlutterSolverBase();
            
            virtual ~FrequencyDomainFlutterSolverBase();
            
            /*!
             *    Initializes the flutter solver for the specified parameters: \p n is the number of generalized coordinates
             *    \p b_ref reference chord, \p rho aerodynamic fluid density, \p es Complex generalized non-Hermitian eigensolver
             */
            void initialize(FESystemUInt n, typename RealOperationType(ValType) b_ref, typename RealOperationType(ValType) rho, 
                            FESystem::Solvers::LinearEigenSolverBase<ValType>& es);
            
            /*!
             *   Gets the current value of reduce frequency for generation of the unsteady aerodynamic matrices: k = omega b / U
             */
            typename RealOperationType(ValType) getCurrentReducedFrequency();
            
            /*!
             *   sets the matrices for preparation of the eigensolution
             */
            virtual void setMatrices(FESystem::Numerics::MatrixBase<typename RealOperationType(ValType)>* mass,
                                     FESystem::Numerics::MatrixBase<typename RealOperationType(ValType)>* damp,
                                     FESystem::Numerics::MatrixBase<typename RealOperationType(ValType)>* stiff,
                                     FESystem::Numerics::MatrixBase<ValType>* aero) = 0;
            
        protected:
            
            /*!
             *    method to initialize matrices
             */
            virtual void initializeMatrices(FESystemUInt n_basis, FESystem::Numerics::MatrixBase<ValType>& A_mat, 
                                            FESystem::Numerics::MatrixBase<ValType>& B_mat)=0;
            
            /*!
             *   boolean to store if this solver is initialized
             */
            FESystemBoolean if_initialized;
            
            /*!
             *    the number of generalized coordinates
             */
            FESystemUInt n_basis;

            /*!
             *    aerodynamic fluid density
             */
            typename RealOperationType(ValType) fluid_rho;
            
            /*!
             *    aerodynamic reference chord
             */
            typename RealOperationType(ValType) aero_b_ref;

            /*!
             *    aerodynamic reduced frequency (omega*b/U)
             */
            typename RealOperationType(ValType) current_k_ref;

            /*!
             *    complex matrices for eigensolver: A matrix
             */
            FESystem::Numerics::MatrixBase<ValType>* A_mat;
            
            /*!
             *    complex matrices for eigensolver: B matrix
             */
            FESystem::Numerics::MatrixBase<ValType>* B_mat;
            
            /*!
             *    complex eigensolver
             */
            FESystem::Solvers::LinearEigenSolverBase<ValType>* eigen_solver; 
        };
    }
}


#endif // __fesystem_frequency_domain_flutter_solver_base_h__
