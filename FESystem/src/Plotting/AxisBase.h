//
//  AxisBase.h
//  FESystem
//
//  Created by Manav Bhatia on 4/14/11.
//  Copyright 2011 . All rights reserved.
//


#ifndef __fesystem_axis_base_h__
#define __fesystem_axis_base_h__

// C++ includes
#include <string>
#include <memory>

// FESystem includes
#include "Base/FESystemExceptions.h"


namespace FESystem
{
    namespace Plotting
    {
        enum AxisKind 
        {
            REAL_AXIS,
            LOG_AXIS
        };

        
        template <typename ValType> 
        class AxisBase
        {
        public: 
            AxisBase();
            
            virtual ~AxisBase();
            
            void setLabel(const std::string& s);
            
            void setRange(ValType v_low, ValType v_up);
            
            std::pair<ValType, ValType> getRange() const;

            FESystemBoolean ifAutomaticRange() const;
            
            const std::string& getLabel() const;

            virtual FESystem::Plotting::AxisKind getAxisKind() = 0;
            
        protected:
            
            std::string label;
            
            FESystemBoolean if_automatic_range;
            
            ValType low;
            
            ValType up;
        };
        
        
        template <typename ValType>
        std::auto_ptr<FESystem::Plotting::AxisBase<ValType> > 
        AxisCreate(FESystem::Plotting::AxisKind a);
        
        
        DeclareException2(InvalidAxisRange, 
                          FESystemDouble, FESystemDouble, 
                          << "Upper limit must be greater than lower limit. \n"
                          << "Lower limit: " << Arg1 << "\n"
                          << "Upper limit: " << Arg2); 

    }
}

#endif //__fesystem_axis_base_h__
