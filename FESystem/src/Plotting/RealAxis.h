//
//  RealAxis.h
//  FESystem
//
//  Created by Manav Bhatia on 4/14/11.
//  Copyright 2011 . All rights reserved.
//


#ifndef __fesystem_real_axis_h__
#define __fesystem_real_axis_h__

// C++ includes
#include <string>

// FESystem includes
#include "Plotting/AxisBase.h"


namespace FESystem
{
    namespace Plotting
    {
        template <typename ValType> 
        class RealAxis: public FESystem::Plotting::AxisBase<ValType>
        {
        public: 
            RealAxis();
            
            virtual ~RealAxis();
            
            FESystem::Plotting::AxisKind getAxisKind();
            
        protected:
            
        };
    }
}

#endif //__fesystem_real_axis_h__

