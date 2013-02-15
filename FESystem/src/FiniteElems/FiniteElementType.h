//
//  FiniteElementType.h
//  FESystem
//
//  Created by Manav Bhatia on 2/14/13.
//
//

#ifndef FESystem_FiniteElementType_h
#define FESystem_FiniteElementType_h

namespace FESystem
{
    namespace FiniteElement
    {
        /*!
         *    Type of finite element
         */
        enum FiniteElementType
        {
            FE_LAGRANGE,
            FE_LEGENDRE
        };
        
        enum DofDistributionType
        {
            POINT_INDEPENDENT,
            UNIFORM,
            UNIFORM_WITH_BOUNDARY,
            SPECIFIED_POINTS
        };
        
        enum ContinuityType
        {
            C0_CONTINUOUS,
            DISCONTINUOUS
        };
    }
}

#endif
