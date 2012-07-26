//
//  LogAxis.cpp
//  FESystem
//
//  Created by Manav Bhatia on 4/18/11.
//  Copyright 2011 . All rights reserved.
//

#include "Plotting/LogAxis.h"


template <typename ValType>        
FESystem::Plotting::LogAxis<ValType>::LogAxis():
FESystem::Plotting::AxisBase<ValType>()
{
    
}


template <typename ValType>        
FESystem::Plotting::LogAxis<ValType>::~LogAxis()
{
    
}



template <typename ValType>        
FESystem::Plotting::AxisKind
FESystem::Plotting::LogAxis<ValType>::getAxisKind()
{
    return FESystem::Plotting::LOG_AXIS;
}

/***************************************************************************************/
// Template instantiations for some generic classes

template class FESystem::Plotting::LogAxis<FESystemDouble>;


/***************************************************************************************/

