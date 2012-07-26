//
//  ModifiedQRFactorization.h
//  FESystem
//
//  Created by Manav Bhatia on 4/8/11.
//  Copyright 2011 . All rights reserved.
//

#ifndef __fesystem_modified_qr_factorization_h__
#define __fesystem_modified_qr_factorization_h__


// FESystem include 
#include "Solvers/Factorizations/MatrixQRFactorizationBase.h"


namespace FESystem
{
    // Forward declerations
    namespace Numerics {template <typename ValType> class VectorBase;}
    namespace Numerics {template <typename ValType> class MatrixBase;}

    namespace Solvers
    {
        
        template <typename ValType> 
        class ModifiedQRFactorization: public FESystem::Solvers::MatrixQRFactorizationBase<ValType>
        {
        public:
            
            ModifiedQRFactorization();
            
            virtual ~ModifiedQRFactorization();
            
            virtual void factorize();
            
        protected:
            
        };
        
    }
}


#endif // __fesystem_classical_qr_factorization_h__
