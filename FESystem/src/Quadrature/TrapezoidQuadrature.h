//
//  TrapezoidQuadrature.h
//  FESystem
//
//  Created by Manav Bhatia on 4/9/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem_trapezoid_quadrature_h__
#define __fesystem_trapezoid_quadrature_h__

// FESystem includes
#include "Quadrature/QuadratureBase.h"


namespace FESystem
{
    namespace Quadrature
    {
        /*!
         *   The TrapezoidQuadrature derives from the QuadratureBase and implements the quadrature rules 
         *   for a trapezoidal integration. 
         */
        class TrapezoidQuadrature: public FESystem::Quadrature::QuadratureBase
        {
        public:
            /*!
             *   Constructor
             */
            TrapezoidQuadrature();
            
            
            virtual ~TrapezoidQuadrature();
            
            /*!  
             *   This implements the procedure to initialize the data structure for the specified dimension 
             *   and integration order. 
             */ 
            virtual void init(FESystemUInt dim, FESystemUInt order);
            
        protected:
            
            
        };
    }
}



#endif // __fesystem_trapezoid_quadrature_h__
