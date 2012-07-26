//
//  PolarFunctionMapping.cpp
//  FESystem
//
//  Created by Manav Bhatia on 4/13/12.
//  Copyright (c) 2012. All rights reserved.
//

// FESystem includes
#include "Functions/PolarFunctionMapping.h"
#include "Numerics/VectorBase.h"
#include "Numerics/DenseMatrix.h"
#include "Base/macros.h"


template <typename ValType>
FESystem::Functions::PolarFunctionMapping<ValType>::PolarFunctionMapping():
FESystem::Functions::FunctionMappingBase<ValType>()
{
    
}


template <typename ValType>
FESystem::Functions::PolarFunctionMapping<ValType>::~PolarFunctionMapping()
{
    
}


template <typename ValType>
std::pair<FESystemUInt, FESystemUInt>
FESystem::Functions::PolarFunctionMapping<ValType>::getDimensions() const
{
    return std::pair<FESystemUInt, FESystemUInt> (2,2);
}


template <typename ValType>
void
FESystem::Functions::PolarFunctionMapping<ValType>::map(const FESystem::Numerics::VectorBase<ValType>& vin, FESystem::Numerics::VectorBase<ValType>& vout) const
{
    FESystemAssert0(vin.getSize() == 2, FESystem::Exception::InvalidValue);
    FESystemAssert0(vout.getSize() == 2, FESystem::Exception::InvalidValue);
    
    ValType r = vin.getVal(0), theta = vin.getVal(2);
    vout.setVal(0, r * cos(theta)); // x
    vout.setVal(1, r * sin(theta)); // y
}


template <typename ValType>
void 
FESystem::Functions::PolarFunctionMapping<ValType>::inverseMap(const FESystem::Numerics::VectorBase<ValType>& vin, FESystem::Numerics::VectorBase<ValType>& vout) const
{
    FESystemAssert0(vin.getSize() == 2, FESystem::Exception::InvalidValue);
    FESystemAssert0(vout.getSize() == 2, FESystem::Exception::InvalidValue);
    
    ValType x = vin.getVal(0), y = vin.getVal(2);
    vout.setVal(0, sqrt(x*x + y*y));  // r
    vout.setVal(1, atan2(y, x));  // theta
}


template <typename ValType>
void
FESystem::Functions::PolarFunctionMapping<ValType>::getMappingJacobian(const FESystem::Numerics::VectorBase<ValType>& vin, FESystem::Numerics::MatrixBase<ValType>& mat) const
{
    std::pair<FESystemUInt, FESystemUInt> s=mat.getSize();
    
    FESystemAssert0(vin.getSize() == 2, FESystem::Exception::InvalidValue);
    FESystemAssert0((s.first==2) && (s.second==2), FESystem::Exception::InvalidValue);
    
    ValType r = vin.getVal(0), theta = vin.getVal(2);
    
    mat.setVal(0,0, cos(theta)); // d x / d r
    mat.setVal(0,1, -r * sin(theta)); // d x / d theta
    
    mat.setVal(1,0, sin(theta)); // d y / d r
    mat.setVal(1,1, r * cos(theta)); // d y / d theta

}


/***************************************************************************************/
// Template instantiations for some generic classes

INSTANTIATE_CLASS_FOR_ONLY_REAL_DATA_TYPES(FESystem::Functions::PolarFunctionMapping);


/***************************************************************************************/




