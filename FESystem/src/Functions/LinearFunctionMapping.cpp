//
//  LinearFunctionMapping.cpp
//  FESystem
//
//  Created by Manav Bhatia on 4/13/12.
//  Copyright (c) 2012. All rights reserved.
//

// FESystem includes
#include "Functions/LinearFunctionMapping.h"
#include "Numerics/DenseMatrix.h"
#include "Numerics/VectorBase.h"
#include "Base/macros.h"


template <typename ValType>
FESystem::Functions::LinearFunctionMapping<ValType>::LinearFunctionMapping():
FESystem::Functions::FunctionMappingBase<ValType>(),
if_initialized(false),
jac(new FESystem::Numerics::DenseMatrix<ValType>())
{
    
}

template <typename ValType>
FESystem::Functions::LinearFunctionMapping<ValType>::~LinearFunctionMapping()
{
    
}
    

template <typename ValType>
void 
FESystem::Functions::LinearFunctionMapping<ValType>::clear()
{
    this->jac->zero();
    this->if_initialized = false;
}
    

//template <typename ValType>
//void
//FESystem::Functions::LinearFunctionMapping<ValType>::reinit(const FESystem::Numerics::MatrixBase<ValType>& m1, const FESystem::Numerics::MatrixBase<ValType>& m2)
//{
//    FESystemAssert0(!this->if_initialized, FESystem::Exception::InvalidState);
//    
//    std::pair<FESystemUInt, FESystemUInt> s1 = m1.getSize(), s2 = m2.getSize();
//    FESystemAssert4(s2.second == s1.second, FESystem::Numerics::MatrixMultiplicationSizeMismatch, s2.first, s2.second, s1.second, s1.first);
//    
//    this->jac->resize(s1.first, s2.first);
//    m2.matrixRightMultiplyTranspose(1.0, m1,*jac);
//    this->if_initialized = true;
//}
    

template <typename ValType>
void
FESystem::Functions::LinearFunctionMapping<ValType>::reinit(const FESystem::Numerics::MatrixBase<ValType>& mat)
{
    FESystemAssert0(!this->if_initialized, FESystem::Exception::InvalidState);

    std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    
    this->jac->resize(s.first, s.second);
    this->jac->copyMatrix(mat);
    this->if_initialized = true;
}


template <typename ValType>
std::pair<FESystemUInt, FESystemUInt>
FESystem::Functions::LinearFunctionMapping<ValType>::getDimensions() const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    std::pair<FESystemUInt, FESystemUInt> s = this->jac->getSize();
    return std::pair<FESystemUInt, FESystemUInt> (s.second, s.first);
}


template <typename ValType>
void
FESystem::Functions::LinearFunctionMapping<ValType>::map(const FESystem::Numerics::VectorBase<ValType>& vin, FESystem::Numerics::VectorBase<ValType>& vout) const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    std::pair<FESystemUInt, FESystemUInt> s = this->jac->getSize();
    
    FESystemAssert3(s.second == vin.getSize(), FESystem::Numerics::MatrixVectorSizeMismatch, s.first, s.second, vin.getSize());
    FESystemAssert3(s.first == vout.getSize(), FESystem::Numerics::MatrixVectorSizeMismatch, s.first, s.second, vout.getSize());
    
    this->jac->rightVectorMultiply(vin, vout);
}



template <typename ValType>
void
FESystem::Functions::LinearFunctionMapping<ValType>::inverseMap(const FESystem::Numerics::VectorBase<ValType>& vin, FESystem::Numerics::VectorBase<ValType>& vout) const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    std::pair<FESystemUInt, FESystemUInt> s = this->jac->getSize();
    
    FESystemAssert3(s.second == vout.getSize(), FESystem::Numerics::MatrixVectorSizeMismatch, s.first, s.second, vout.getSize());
    FESystemAssert3(s.first == vin.getSize(), FESystem::Numerics::MatrixVectorSizeMismatch, s.first, s.second, vin.getSize());
    
    this->jac->leftVectorMultiply(vin, vout);
}




template <typename ValType>
void
FESystem::Functions::LinearFunctionMapping<ValType>::getMappingJacobian(const FESystem::Numerics::VectorBase<ValType>& vin, FESystem::Numerics::MatrixBase<ValType>& mat) const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    mat.copyMatrix(*(this->jac));
}
  


template <typename ValType>
void
FESystem::Functions::LinearFunctionMapping<ValType>::getFunctionDerivative(const std::vector<FESystemUInt>& derivative_orders, 
                                                                           const FESystem::Numerics::VectorBase<ValType>& val, 
                                                                           FESystem::Numerics::VectorBase<ValType>& deriv) const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    // this is a linear function and therefore the sensitivity does not change with respect to the value
    
    std::pair<FESystemUInt, FESystemUInt> s=this->getDimensions();
    FESystemAssert2(derivative_orders.size() == s.first, FESystem::Exception::DimensionsDoNotMatch, derivative_orders.size(), s.first);
    FESystemAssert2(deriv.getSize() == s.second, FESystem::Exception::DimensionsDoNotMatch, deriv.getSize(), s.second);
    
    // get the value of the derivative order
    FESystemUInt order=0, id=0;
    for (FESystemUInt i=0; i<s.first; i++)
        if (order < derivative_orders[i])
        {
            order = derivative_orders[i];
            id = i;
        }
    // if the derivative order is higher than 1, then the value is zero
    FESystemAssert0(order > 0, FESystem::Exception::InvalidValue);
    
    if (order == 1)
        this->jac->getColumnVals(id, 0, s.second, deriv);
    else 
        deriv.zero();
}


/***************************************************************************************/
// Template instantiations for some generic classes

INSTANTIATE_CLASS_FOR_ALL_DATA_TYPES(FESystem::Functions::LinearFunctionMapping);


/***************************************************************************************/




