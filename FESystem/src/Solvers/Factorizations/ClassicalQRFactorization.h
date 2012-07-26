//
//  QRFactorization.h
//  FESystem
//
//  Created by Manav Bhatia on 3/27/11.
//  Copyright 2011 . All rights reserved.
//

#ifndef __fesystem_classical_qr_factorization_h__
#define __fesystem_classical_qr_factorization_h__


// FESystem include 
#include "Solvers/Factorizations/MatrixQRFactorizationBase.h"


namespace FESystem
{
    namespace Solvers
    {
        // Forward declerations
        template <typename ValType> class VectorBase;
        template <typename ValType> class MatrixBase;
        
        template <typename ValType> 
        class ClassicalQRFactorization: public FESystem::Solvers::MatrixQRFactorizationBase<ValType>
        {
        public:
            
            ClassicalQRFactorization();
            
            virtual ~ClassicalQRFactorization();
            
            virtual void factorize();

        protected:
            
        };
        
    }
}


#endif // __fesystem_classical_qr_factorization_h__
