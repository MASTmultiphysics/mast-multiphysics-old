//
//  AxisBase.cpp
//  FESystem
//
//  Created by Manav Bhatia on 4/14/11.
//  Copyright 2011 . All rights reserved.
//

// FESystem includes
#include "Plotting/AxisBase.h"
#include "Plotting/RealAxis.h"
#include "Plotting/LogAxis.h"


template <typename ValType>        
FESystem::Plotting::AxisBase<ValType>::AxisBase():
label(),
if_automatic_range(true),
low(0),
up(1)
{
    
}
            

template <typename ValType>        
FESystem::Plotting::AxisBase<ValType>::~AxisBase()
{
    
}
            

template <typename ValType>        
void
FESystem::Plotting::AxisBase<ValType>::setLabel(const std::string& s)
{
    this->label = s;
}
            


template <typename ValType>        
void
FESystem::Plotting::AxisBase<ValType>::setRange(ValType v_low, ValType v_up)
{
    FESystemAssert2(v_up > v_low, 
                    FESystem::Plotting::InvalidAxisRange,
                    v_low, v_up);
    
    this->low = v_low;
    this->up = v_up;
    
    this->if_automatic_range = false;
}



template <typename ValType>        
std::pair<ValType, ValType>
FESystem::Plotting::AxisBase<ValType>::getRange() const
{
    std::pair<ValType, ValType> p(this->low, this->up);
    return p;
}



template <typename ValType>        
FESystemBoolean
FESystem::Plotting::AxisBase<ValType>::ifAutomaticRange() const
{
    return this->if_automatic_range;
}



template <typename ValType>        
const std::string&
FESystem::Plotting::AxisBase<ValType>::getLabel() const
{
    return this->label;
}


template <typename ValType>
std::auto_ptr<FESystem::Plotting::AxisBase<ValType> >
FESystem::Plotting::AxisCreate(FESystem::Plotting::AxisKind a)
{
    std::auto_ptr<FESystem::Plotting::AxisBase<ValType> > r;
    
    switch (a)
    {
        case FESystem::Plotting::REAL_AXIS:
            r.reset(new FESystem::Plotting::RealAxis<ValType>() );
            break;
            
        case FESystem::Plotting::LOG_AXIS:
            r.reset(new FESystem::Plotting::LogAxis<ValType>() );
            break;
            
    }
    
    return r;
}


/***************************************************************************************/
// Template instantiations for some generic classes

template std::auto_ptr<FESystem::Plotting::AxisBase<FESystemDouble> > 
FESystem::Plotting::AxisCreate(FESystem::Plotting::AxisKind a);

template class FESystem::Plotting::AxisBase<FESystemDouble>;


/***************************************************************************************/

