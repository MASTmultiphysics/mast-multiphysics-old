//
//  LogAxis.h
//  FESystem
//
//  Created by Manav Bhatia on 4/18/11.
//  Copyright 2011 . All rights reserved.
//

#ifndef __fesystem_log_axis_h__
#define __fesystem_log_axis_h__


// C++ includes
#include <string>

// FESystem includes
#include "Plotting/AxisBase.h"


namespace FESystem
{
    namespace Plotting
    {
        template <typename ValType> 
        class LogAxis: public FESystem::Plotting::AxisBase<ValType>
        {
        public: 
            LogAxis();
            
            virtual ~LogAxis();
            
            FESystem::Plotting::AxisKind getAxisKind();
            
        protected:
            
        };
    }
}



#endif // __fesystem_log_axis_h__
