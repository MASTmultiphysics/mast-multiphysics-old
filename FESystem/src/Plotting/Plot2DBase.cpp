//
//  Plot2DBase.cpp
//  FESystem
//
//  Created by Manav Bhatia on 4/14/11.
//  Copyright 2011 . All rights reserved.
//

// FESystem includes
#include "Plotting/Plot2DBase.h"



template <typename ValType>
FESystem::Plotting::Plot2DBase<ValType>::Plot2DBase(FESystem::Plotting::AxisKind& a1,
                                                    FESystem::Plotting::AxisKind& a2)
{
    this->axes.resize(this->getDimension());
    this->axes[0] = FESystem::Plotting::AxisCreate<ValType>(a1).release();
    this->axes[1] = FESystem::Plotting::AxisCreate<ValType>(a2).release();
}



template <typename ValType>
FESystem::Plotting::Plot2DBase<ValType>::~Plot2DBase()
{
    for (FESystemUInt i=0; i<this->getDimension(); i++)
        if (this->axes[i] != NULL)
            delete this->axes[i];
}



template <typename ValType>
FESystem::Plotting::AxisBase<ValType>& 
FESystem::Plotting::Plot2DBase<ValType>::getAxis(FESystemUInt i)
{
    FESystemAssert2(i < this->axes.size(), FESystem::Exception::IndexOutOfBound, 
                    i, this->axes.size());
    
    FESystemAssert0(this->axes[i]!=NULL, FESystem::Exception::NULLQuantity);
    
    return *(this->axes[i]);
}


template <typename ValType>
FESystemUInt 
FESystem::Plotting::Plot2DBase<ValType>::getDimension()
{
    return 2;
}


template <typename ValType>
void
FESystem::Plotting::Plot2DBase<ValType>::setTitle(const std::string& s)
{
    this->title = s;
}



/***************************************************************************************/
// Template instantiations for some generic classes

template class FESystem::Plotting::Plot2DBase<FESystemDouble>;


/***************************************************************************************/




