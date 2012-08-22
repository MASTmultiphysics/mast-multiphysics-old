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
            
        protected:
            
            /*!
             *    method to initialize matrices
             */
            virtual void initializeMatrices(FESystemUInt n_basis, FESystem::Numerics::MatrixBase<ValType>& A_mat, 
                                            FESystem::Numerics::MatrixBase<ValType>& B_mat)=0;
            
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
