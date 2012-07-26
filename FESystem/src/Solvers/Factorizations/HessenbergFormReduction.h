//
//  HessenbergFormReduction.h
//  FESystem
//
//  Created by Manav Bhatia on 5/8/11.
//  Copyright 2011 . All rights reserved.
//

#ifndef __fesystem_hessenberg_form_reduction_h__
#define __fesystem_hessenberg_form_reduction_h__

// C++ includes
#include <memory>

// FESystem includes
#include "Base/FESystemTypes.h"


namespace FESystem
{    
    // Forward declerations
    namespace Numerics {template <typename ValType> class MatrixBase;}

    namespace Solvers
    {
        
        /*! 
         *   This class calculates the Hessenberg form of a given matrix \f$ A \f$
         *    \f[ A = Q H Q^* \f]
         *   where \f$ Q \f$ is the orthonormal matrix
         *   and \f$ H \f$ is the Hessenberg matrix
         */
        template <typename ValType>
        class HessenbergFormReduction
        {
            public:
            
            HessenbergFormReduction();
            
            ~HessenbergFormReduction();
            
            void initializeMatrices();
            
            void setMatrix(FESystem::Numerics::MatrixBase<ValType>* m);
            
            FESystem::Numerics::MatrixBase<ValType>& getMatrix();
            
            FESystem::Numerics::MatrixBase<ValType>& getQMatrix();
            
            FESystem::Numerics::MatrixBase<ValType>& getHMatrix();
            
            void factorize();
            
            protected:
            
            FESystem::Numerics::MatrixBase<ValType>* mat;
            
            FESystemBoolean factorization_complete;
            
            std::auto_ptr<FESystem::Numerics::MatrixBase<ValType> > Q_mat;
            std::auto_ptr<FESystem::Numerics::MatrixBase<ValType> > H_mat;
            
        };
    }
}


#endif // __fesystem_hessenberg_form_reduction_h__
