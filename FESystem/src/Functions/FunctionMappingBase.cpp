//
//  FunctionMappingBase.cpp
//  FESystem
//
//  Created by Manav Bhatia on 4/14/12.
//  Copyright (c) 2012. All rights reserved.
//

// FESystem includes
#include "Functions/FunctionMappingBase.h"
#include "Numerics/DenseMatrix.h"
#include "Base/macros.h"


template <typename ValType>
FESystem::Functions::FunctionMappingBase<ValType>::FunctionMappingBase()
{

}


template <typename ValType>
FESystem::Functions::FunctionMappingBase<ValType>::~FunctionMappingBase()
{
    
}



template <typename ValType>
FESystem::Functions::DiscreteFunctionMappingBase<ValType>::DiscreteFunctionMappingBase():
FESystem::Functions::FunctionMappingBase<ValType>(),
if_discrete_function_initialized(false)
{
    this->discrete_function_values = new FESystem::Numerics::DenseMatrix<ValType>;
}
            

template <typename ValType>
FESystem::Functions::DiscreteFunctionMappingBase<ValType>::~DiscreteFunctionMappingBase()
{
    delete this->discrete_function_values;
}
            

template <typename ValType>
void
FESystem::Functions::DiscreteFunctionMappingBase<ValType>::clear()
{
    this->discrete_function_values->zero();
    this->if_discrete_function_initialized = false;
}


template <typename ValType>
void
FESystem::Functions::DiscreteFunctionMappingBase<ValType>::reinit(FESystemUInt dim, FESystemUInt n_pts)
{
    FESystemAssert0(!this->if_discrete_function_initialized, FESystem::Exception::InvalidState);

    this->discrete_function_values->resize(dim, n_pts);
}


template <typename ValType>
void
FESystem::Functions::DiscreteFunctionMappingBase<ValType>::inverseMap(const FESystem::Numerics::VectorBase<ValType>& vin, FESystem::Numerics::VectorBase<ValType>& vout) const
{
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}


template <typename ValType>
void 
FESystem::Functions::DiscreteFunctionMappingBase<ValType>::setDiscreteFunctionValues(const FESystem::Numerics::MatrixBase<ValType>& vals)
{
    FESystemAssert0(!this->if_discrete_function_initialized, FESystem::Exception::InvalidState);

    std::pair<FESystemUInt, FESystemUInt> s = this->getDimensions(), s_discrete=this->discrete_function_values->getSize(), s_in=vals.getSize();
    
    FESystemUInt n_funcs = this->getNDiscreteFunctions();
    
    // make sure that the data structures are initialized 
    FESystemAssert4((s_discrete.first == s.second) && (s_discrete.second==n_funcs), FESystem::Numerics::MatrixSizeMismatch, 
                    s_discrete.first, s_discrete.second, s.second, n_funcs);
    
    // make sure that the dimensions are appropriate
    FESystemAssert2(vals.getSize().first == s.second, FESystem::Exception::DimensionsDoNotMatch, vals.getSize().first, s.second);

    this->discrete_function_values->copyMatrixValues(vals);
    
    this->if_discrete_function_initialized = true;
}



/***************************************************************************************/
// Template instantiations for some generic classes

INSTANTIATE_CLASS_FOR_ALL_DATA_TYPES(FESystem::Functions::FunctionMappingBase);
INSTANTIATE_CLASS_FOR_ALL_DATA_TYPES(FESystem::Functions::DiscreteFunctionMappingBase);


/***************************************************************************************/

