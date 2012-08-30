//
//  HouseHolderTriangulation.h
//  FESystem
//
//  Created by Manav Bhatia on 3/27/11.
//  Copyright 2011 . All rights reserved.
//
#ifndef __fesystem_householder_triangulation_h__
#define __fesystem_householder_triangulation_h__

// C++ includes
#include <memory>

// FESystem include 
#include "Solvers/Factorizations/MatrixQRFactorizationBase.h"

namespace FESystem
{
    // Forward declerations
    namespace Numerics {template <typename ValType> class VectorBase;}
    namespace Numerics {template <typename ValType> class MatrixBase;}

    namespace FactorizationSolvers
    {
        
        template <typename ValType> 
        class HouseholderTriangulation: public FESystem::FactorizationSolvers::MatrixQRFactorizationBase<ValType>
        {
        public:
            
            HouseholderTriangulation();
            
            virtual ~HouseholderTriangulation();
            
            virtual void factorize();
            
        protected:
            
        };
    }
}


#endif // __fesystem_householder_triangulation_h__
