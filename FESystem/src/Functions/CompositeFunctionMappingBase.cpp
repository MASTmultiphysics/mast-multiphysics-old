//
//  CompositeFunctionMappingBase.cpp
//  FESystem
//
//  Created by Manav Bhatia on 4/14/12.
//  Copyright (c) 2012. All rights reserved.
//

// FESystem includes
#include "Functions/CompositeFunctionMappingBase.h"
#include "Numerics/LocalVector.h"
#include "Numerics/DenseMatrix.h"
#include "Base/macros.h"



template <typename ValType>
FESystem::Functions::CompositeFunctionMappingBase<ValType>::CompositeFunctionMappingBase():
FESystem::Functions::FunctionMappingBase<ValType>(),
if_initialized(false),
original_function_mapping(NULL),
sub_function_mapping(NULL)
{
    
}
    
template <typename ValType>
FESystem::Functions::CompositeFunctionMappingBase<ValType>::~CompositeFunctionMappingBase()
{
    
}
    

template <typename ValType>
std::pair<FESystemUInt, FESystemUInt> 
FESystem::Functions::CompositeFunctionMappingBase<ValType>::getDimensions() const
{
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}
    
template <typename ValType>
void 
FESystem::Functions::CompositeFunctionMappingBase<ValType>::clear()
{
    this->if_initialized = false;
    this->original_function_mapping = NULL;
    this->sub_function_mapping = NULL;
}
    
template <typename ValType>
void 
FESystem::Functions::CompositeFunctionMappingBase<ValType>::reinit(const FESystem::Functions::FunctionMappingBase<ValType>& original,
                                                                   const FESystem::Functions::FunctionMappingBase<ValType>& sub)
{
    FESystemAssert0(!this->if_initialized, FESystem::Exception::InvalidState);
    // make sure that the dimensions of the two mappings are consistent for X1
    std::pair<FESystemUInt, FESystemUInt> s2 = sub.getDimensions(), s3 = original.getDimensions();
    FESystemAssert2(s2.first == s3.first, FESystem::Functions::InconsistentDimensions, s2.first, s3.first);
    
    this->original_function_mapping = &original;
    this->sub_function_mapping = &sub;
    this->if_initialized = true;
}
    

template <typename ValType>
const FESystem::Functions::FunctionMappingBase<ValType>& 
FESystem::Functions::CompositeFunctionMappingBase<ValType>::getOriginalFunctionMapping() const
{
    // make sure this is initialized 
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    return *(this->original_function_mapping);
}
    
template <typename ValType>
const FESystem::Functions::FunctionMappingBase<ValType>& 
FESystem::Functions::CompositeFunctionMappingBase<ValType>::getSubFunctionMapping() const
{
    // make sure this is initialized 
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    return *(this->sub_function_mapping);
}
    
template <typename ValType>
void
FESystem::Functions::CompositeFunctionMappingBase<ValType>::map(const FESystem::Numerics::VectorBase<ValType>& vin, FESystem::Numerics::VectorBase<ValType>& vout) const
{
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);    
}

template <typename ValType>
void 
FESystem::Functions::CompositeFunctionMappingBase<ValType>::getMappingJacobian(const FESystem::Numerics::VectorBase<ValType>& vin, FESystem::Numerics::MatrixBase<ValType>& mat) const
{
    // make sure this is initialized 
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    

}

template <typename ValType>
void 
FESystem::Functions::CompositeFunctionMappingBase<ValType>::getFunctionDerivative(const std::vector<FESystemUInt>& derivative_orders, 
                                                                                  const FESystem::Numerics::VectorBase<ValType>& val, 
                                                                                  FESystem::Numerics::VectorBase<ValType>& deriv) const
{
    // make sure this is initialized 
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    // make sure that the dimensions of the vector is consistent
    std::pair<FESystemUInt, FESystemUInt> s2 = this->sub_function_mapping->getDimensions(), s3 = this->original_function_mapping->getDimensions();
    FESystemAssert2(s2.second == derivative_orders.size(), FESystem::Functions::InconsistentDimensions, derivative_orders.size(), s2.second);
    FESystemAssert2(val.getSize() == s3.first, FESystem::Functions::InconsistentDimensions, val.getSize(), s3.first);
    
    // check what is the order of the derivative
    FESystemUInt order = 0, coord_num=0;
    for (FESystemUInt i=0; i<derivative_orders.size(); i++) 
    {
        order += derivative_orders[i];
        coord_num = i; // coordinate number for which this derivative is needed
    }
    
    // zero the vector before adding to it
    deriv.zero();
    
    switch (order) {
        case 1:
        {
            // d X3_i / d X2_j =  (d X3_i / d X1_k)  (d X1_k / d X2_j)
            FESystem::Numerics::LocalVector<ValType> dX3_dX1;
            FESystem::Numerics::DenseMatrix<ValType> dX2_dX1_jac, dX2_dX1_jac_inv;
            dX3_dX1.resize(s3.second);
            dX2_dX1_jac.resize(s3.second, s3.first);
            dX2_dX1_jac_inv.resize(s3.second, s3.first);
            std::vector<FESystemUInt> x3_deriv(s3.first);
            
            this->sub_function_mapping->getMappingJacobian(val, dX2_dX1_jac);
            dX2_dX1_jac.getInverse(dX2_dX1_jac_inv);
            
            for (FESystemUInt i=0; i<s3.first; i++)
            {
                for (FESystemUInt j=0; j<s3.first; j++)  x3_deriv[j] = 0;
                x3_deriv[i] = 1;
                dX3_dX1.zero();
                this->original_function_mapping->getFunctionDerivative(x3_deriv, val, dX3_dX1);
                deriv.add(dX2_dX1_jac_inv.getVal(i,coord_num), dX3_dX1);  
            }
        }
            break;
            
        case 2: // this is to be implemented 
        case 0: // this method is used specifically for derivatives, and not for function values
        default:
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
    }
}



template <typename ValType>
FESystem::Functions::DiscreteCompositeMappingFunctionBase<ValType>::DiscreteCompositeMappingFunctionBase():
FESystem::Functions::CompositeFunctionMappingBase<ValType>(),
original_discrete_function_mapping(NULL),
sub_discrete_function_mapping(NULL)
{
    
}

template <typename ValType>
FESystem::Functions::DiscreteCompositeMappingFunctionBase<ValType>::~DiscreteCompositeMappingFunctionBase()
{
    
}


template <typename ValType>
void 
FESystem::Functions::DiscreteCompositeMappingFunctionBase<ValType>::clear()
{
    this->original_discrete_function_mapping = NULL;
    this->sub_discrete_function_mapping = NULL;
    FESystem::Functions::CompositeFunctionMappingBase<ValType>::clear();
}

template <typename ValType>
void 
FESystem::Functions::DiscreteCompositeMappingFunctionBase<ValType>::reinit(const FESystem::Functions::DiscreteFunctionMappingBase<ValType>& original,
                                                                           const FESystem::Functions::DiscreteFunctionMappingBase<ValType>& sub)
{
    FESystem::Functions::CompositeFunctionMappingBase<ValType>::reinit(original, sub);
    
    this->original_discrete_function_mapping = &original;
    this->sub_discrete_function_mapping = &sub;
}


template <typename ValType>
const FESystem::Functions::DiscreteFunctionMappingBase<ValType>& 
FESystem::Functions::DiscreteCompositeMappingFunctionBase<ValType>::getOriginalFunctionMapping() const
{
    // make sure this is initialized 
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    return *(this->original_discrete_function_mapping);
}

template <typename ValType>
const FESystem::Functions::DiscreteFunctionMappingBase<ValType>& 
FESystem::Functions::DiscreteCompositeMappingFunctionBase<ValType>::getSubFunctionMapping() const
{
    // make sure this is initialized 
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    return *(this->sub_discrete_function_mapping);
}


template <typename ValType>
void
FESystem::Functions::DiscreteCompositeMappingFunctionBase<ValType>::inverseMap(const FESystem::Numerics::VectorBase<ValType>& vin, 
                                                                               FESystem::Numerics::VectorBase<ValType>& vout) const
{
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}



template <typename ValType>
void 
FESystem::Functions::DiscreteCompositeMappingFunctionBase<ValType>::getDiscreteFunction(const FESystem::Numerics::VectorBase<ValType>& vin, 
                                                                                        FESystem::Numerics::VectorBase<ValType>& vout) const
{
    // make sure this is initialized 
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    // make sure that the dimensions of the vector is consistent
    FESystemUInt n_funcs = this->original_discrete_function_mapping->getNDiscreteFunctions();
    std::pair<FESystemUInt, FESystemUInt> s2 = this->sub_function_mapping->getDimensions(), s3 = this->original_function_mapping->getDimensions();
    FESystemAssert2(vin.getSize() == s3.second, FESystem::Functions::InconsistentDimensions, vin.getSize(), s3.second);
    FESystemAssert2(vout.getSize() == n_funcs, FESystem::Functions::InconsistentDimensions, vout.getSize(), n_funcs);
    
    this->original_discrete_function_mapping->getDiscreteFunction(vin, vout);
}


template <typename ValType>
void 
FESystem::Functions::DiscreteCompositeMappingFunctionBase<ValType>::getDiscreteFunctionDerivative(const std::vector<FESystemUInt>& derivative_orders, 
                                                                                                  const FESystem::Numerics::VectorBase<ValType>& vin, 
                                                                                                  FESystem::Numerics::VectorBase<ValType>& deriv) const
{
    // make sure this is initialized 
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    // make sure that the dimensions of the vector is consistent
    FESystemUInt n_funcs = this->original_discrete_function_mapping->getNDiscreteFunctions();
    
    std::pair<FESystemUInt, FESystemUInt> s2 = this->sub_function_mapping->getDimensions(), s3 = this->original_function_mapping->getDimensions();
    FESystemAssert2(vin.getSize() == s3.second, FESystem::Functions::InconsistentDimensions, vin.getSize(), s3.second);
    FESystemAssert2(s2.second == derivative_orders.size(), FESystem::Functions::InconsistentDimensions, derivative_orders.size(), s2.second);
    FESystemAssert2(deriv.getSize() == n_funcs, FESystem::Functions::InconsistentDimensions, deriv.getSize(), n_funcs);
    
    // check what is the order of the derivative
    FESystemUInt order = 0, coord_num=0;
    for (FESystemUInt i=0; i<derivative_orders.size(); i++) 
    {
        order += derivative_orders[i];
        if (derivative_orders[i] > 0)
            coord_num = i; // coordinate number for which this derivative is needed
    }
    
    // zero the vector before adding to it
    deriv.zero();
    
    switch (order) {
        case 1:
        {
            // d X3_i / d X2_j =  (d X3_i / d X1_k)  (d X1_k / d X2_j)
            static FESystem::Numerics::LocalVector<ValType> dX3_dX1;
            static FESystem::Numerics::DenseMatrix<ValType> dX2_dX1_jac, dX2_dX1_jac_inv;
            dX3_dX1.resize(n_funcs);
            dX2_dX1_jac.resize(s3.second, s3.first);
            dX2_dX1_jac_inv.resize(s3.second, s3.first);
            static std::vector<FESystemUInt> x3_deriv;
            if (x3_deriv.size() != s3.first)
                x3_deriv.resize(s3.first);
            
            this->sub_function_mapping->getMappingJacobian(vin, dX2_dX1_jac);
            dX2_dX1_jac.getInverse(dX2_dX1_jac_inv);
            
            for (FESystemUInt i=0; i<s3.first; i++)
            {
                for (FESystemUInt j=0; j<s3.first; j++)  x3_deriv[j] = 0;
                x3_deriv[i] = 1;
                dX3_dX1.zero();
                this->original_discrete_function_mapping->getDiscreteFunctionDerivative(x3_deriv, vin, dX3_dX1);
                deriv.add(dX2_dX1_jac_inv.getVal(i,coord_num), dX3_dX1);  
            }
        }
            break;
            
        case 2: // this is to be implemented 
        case 0: // this method is used specifically for derivatives, and not for function values
        default:
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
    }
}


/***************************************************************************************/
// Template instantiations for some generic classes

INSTANTIATE_CLASS_FOR_ALL_DATA_TYPES(FESystem::Functions::CompositeFunctionMappingBase);
INSTANTIATE_CLASS_FOR_ALL_DATA_TYPES(FESystem::Functions::DiscreteCompositeMappingFunctionBase);


/***************************************************************************************/


