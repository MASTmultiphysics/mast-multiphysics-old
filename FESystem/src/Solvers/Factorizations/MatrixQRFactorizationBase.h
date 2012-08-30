//
//  MatrixQRFactorizationBase.h
//  FESystem
//
//  Created by Manav Bhatia on 4/8/11.
//  Copyright 2011 . All rights reserved.
//

#ifndef __fesystem_matrix_qr_factorization_base_h__
#define __fesystem_matrix_qr_factorization_base_h__

// C++ includes
#include <memory>

// FESystem include 
#include "Base/FESystemExceptions.h"

namespace FESystem
{
    // Forward declerations
    namespace Numerics {template <typename ValType> class VectorBase;}
    namespace Numerics {template <typename ValType> class MatrixBase;}

    namespace FactorizationSolvers
    {
        
        template <typename ValType> 
        class MatrixQRFactorizationBase
        {
        public:
            
            MatrixQRFactorizationBase();
            
            virtual ~MatrixQRFactorizationBase();
            
            void setMatrix(const FESystem::Numerics::MatrixBase<ValType>* m);
            
            const FESystem::Numerics::MatrixBase<ValType>& getMatrix() const;
            
            FESystem::Numerics::MatrixBase<ValType>& getQMatrix();
            
            FESystem::Numerics::MatrixBase<ValType>& getRMatrix();
            
            virtual void clear();
            
            virtual void factorize() = 0;
            
        protected:
            
            void initializeMatrices();
            
            FESystemBoolean factorization_complete;
            
            const FESystem::Numerics::MatrixBase<ValType>* mat;
            
            std::auto_ptr<FESystem::Numerics::MatrixBase<ValType> > Q_mat;
            std::auto_ptr<FESystem::Numerics::MatrixBase<ValType> > R_mat;
        };
        
        DeclareException0(MatrixNotSetBeforeFactorization, 
                          << "Matrix must be set before factorization\n");
        
        DeclareException0(FactorizationNotComplete, 
                          << "Factorization must be complete before access to Q or R matrices.\n");
        
    }
}


#endif // __fesystem_matrix_factorization_base_h__
