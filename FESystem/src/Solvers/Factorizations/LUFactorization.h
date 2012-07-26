//
//  LUFactorization.h
//  FESystem
//
//  Created by Manav Bhatia on 6/17/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem_lu_factorization_h__
#define __fesystem_lu_factorization_h__

// C++ includes
#include <memory>

// FESystem include 
#include "Base/FESystemExceptions.h"

namespace FESystem
{
    // Forward declerations
    namespace Numerics {template <typename ValType> class VectorBase;}
    namespace Numerics {template <typename ValType> class MatrixBase;}
    namespace Numerics {class SparsityPattern;}
    
    namespace Solvers
    {
        
        template <typename ValType> 
        class LUFactorization
        {
        public:
            
            LUFactorization();
            
            virtual ~LUFactorization();
            
            void setMatrix(const FESystem::Numerics::MatrixBase<ValType>* m);
            
            const FESystem::Numerics::MatrixBase<ValType>& getMatrix() const;

            const FESystem::Numerics::MatrixBase<ValType>& getLMatrix() const;
            
            const FESystem::Numerics::MatrixBase<ValType>& getUMatrix() const;
            
            virtual void clear();
                        
        protected:
            
            void initializeMatrices();
            
            FESystemBoolean factorization_complete;
            
            const FESystem::Numerics::MatrixBase<ValType>* mat;
            
            std::auto_ptr<FESystem::Numerics::MatrixBase<ValType> > l_mat;
            std::auto_ptr<FESystem::Numerics::MatrixBase<ValType> > u_mat;
            std::auto_ptr<FESystem::Numerics::SparsityPattern> l_sparsity_pattern;
            std::auto_ptr<FESystem::Numerics::SparsityPattern> u_sparsity_pattern;
        };
    }
}


#endif // __fesystem_lu_factorization_h__

