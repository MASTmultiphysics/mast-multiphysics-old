//
//  RealAxis.cpp
//  FESystem
//
//  Created by Manav Bhatia on 4/14/11.
//  Copyright 2011 . All rights reserved.
//

#include "Plotting/RealAxis.h"

template <typename ValType>        
FESystem::Plotting::RealAxis<ValType>::RealAxis():
FESystem::Plotting::AxisBase<ValType>()
{
    
}


template <typename ValType>        
FESystem::Plotting::RealAxis<ValType>::~RealAxis()
{
    
}



template <typename ValType>        
FESystem::Plotting::AxisKind
FESystem::Plotting::RealAxis<ValType>::getAxisKind()
{
    return FESystem::Plotting::REAL_AXIS;
}

/***************************************************************************************/
// Template instantiations for some generic classes

template class FESystem::Plotting::RealAxis<FESystemDouble>;


/***************************************************************************************/

