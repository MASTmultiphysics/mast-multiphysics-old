//
//  LegendreFunctionMapping.cpp
//  FESystem
//
//  Created by Manav Bhatia on 2/12/13.
//
//

// FESystem includes
#include "Functions/LegendreFunctionMapping.h"
#include "Functions/LegendreFunction.h"
#include "Numerics/LocalVector.h"
#include "Numerics/MatrixBase.h"
#include "Base/macros.h"


template <typename ValType>
FESystem::Functions::LegendreFunctionMapping<ValType>::LegendreFunctionMapping():
FESystem::Functions::DiscreteFunctionMappingBase<ValType>(),
dimension(0),
order(0),
if_initialized(false)
{
    
}

template <typename ValType>
FESystem::Functions::LegendreFunctionMapping<ValType>::~LegendreFunctionMapping()
{
}

template <typename ValType>
void
FESystem::Functions::LegendreFunctionMapping<ValType>::clear()
{
    this->legendre_function=NULL;
    this->dimension=0;
    this->order=0;
    this->if_initialized = false;
    FESystem::Functions::DiscreteFunctionMappingBase<ValType>::clear();
}

template <typename ValType>
void
FESystem::Functions::LegendreFunctionMapping<ValType>::reinit(const FESystemUInt dim, const FESystemUInt o, const FESystem::Functions::LegendreFunction<ValType>& func)
{
    FESystemAssert0(!this->if_initialized, FESystem::Exception::InvalidState);
    
    // make sure that the provided values are correct
    FESystemAssert0(dim > 0, FESystem::Exception::InvalidValue);
    
    // the dimension
    this->dimension = dim;
    this->order = order;
    this->legendre_function = &func;
    
    // the number of shape functions depends on the order and maximum dimension
    FESystem::Functions::DiscreteFunctionMappingBase<ValType>::reinit(dim, pow(this->order+1, dim));
    
    this->if_initialized = true;
}


template <typename ValType>
std::pair<FESystemUInt, FESystemUInt>
FESystem::Functions::LegendreFunctionMapping<ValType>::getDimensions() const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    return std::pair<FESystemUInt, FESystemUInt> (this->dimension, this->dimension);
}


template <typename ValType>
FESystemUInt
FESystem::Functions::LegendreFunctionMapping<ValType>::getNDiscreteFunctions() const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    return pow(this->order+1, this->dimension);
}


template <typename ValType>
void
FESystem::Functions::LegendreFunctionMapping<ValType>::map(const FESystem::Numerics::VectorBase<ValType>& vin, FESystem::Numerics::VectorBase<ValType>& vout) const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    FESystemAssert0(this->if_discrete_function_initialized, FESystem::Exception::InvalidState);
    
    // make sure that the dimensions are consistent
    std::pair<FESystemUInt, FESystemUInt> s=this->getDimensions(), s_discrete=this->discrete_function_values->getSize();
    FESystemAssert2(vin.getSize() == s.first, FESystem::Exception::DimensionsDoNotMatch, vin.getSize(), s.first);
    FESystemAssert2(vout.getSize() == s.second, FESystem::Exception::DimensionsDoNotMatch, vout.getSize(), s.second);
    
    FESystemUInt n_pts = this->getNDiscreteFunctions();
    
    FESystemAssert4(((s_discrete.first==s.second)&&(s_discrete.second == n_pts)), FESystem::Numerics::MatrixSizeMismatch, s_discrete.first, s_discrete.second, s.second, n_pts);
    
    FESystem::Numerics::LocalVector<ValType> vec;
    vec.resize(n_pts);
    
    this->getDiscreteFunction(vin, vec);
    
    // now that the function values are available, calculate the interpolated quantity
    ValType out_val = 0.0;
    for (FESystemUInt i_dim=0; i_dim<s.second; i_dim++)
    {
        out_val = 0.0;
        for (FESystemUInt i_pts=0; i_pts<n_pts; i_pts++)
            out_val += vec.getVal(i_pts)*this->discrete_function_values->getVal(i_dim,i_pts);
        vout.setVal(i_dim, out_val);
    }
}

template <typename ValType>
void
FESystem::Functions::LegendreFunctionMapping<ValType>::getMappingJacobian(const FESystem::Numerics::VectorBase<ValType>& vin, FESystem::Numerics::MatrixBase<ValType>& mat) const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    FESystemAssert0(this->if_discrete_function_initialized, FESystem::Exception::InvalidState);
    
    // make sure that the dimensions are consistent
    std::pair<FESystemUInt, FESystemUInt> s=this->getDimensions(), smat = mat.getSize(), s_discrete=this->discrete_function_values->getSize();
    FESystemAssert2(vin.getSize() == s.first, FESystem::Exception::DimensionsDoNotMatch, vin.getSize(), s.first);
    FESystemAssert2(s.first == smat.second, FESystem::Exception::DimensionsDoNotMatch, smat.second, s.first);
    FESystemAssert2(s.second == smat.first, FESystem::Exception::DimensionsDoNotMatch, smat.first, s.second);
    
    FESystemUInt n_pts = this->getNDiscreteFunctions();
    
    FESystemAssert4(((s_discrete.first==s.second)&&(s_discrete.second == n_pts)), FESystem::Numerics::MatrixSizeMismatch, s_discrete.first, s_discrete.second, s.second, n_pts);
    
    FESystem::Numerics::LocalVector<ValType> deriv;
    deriv.resize(s.second);
    
    std::vector<FESystemUInt> deriv_order;
    if (deriv_order.size() != s.first)
        deriv_order.resize(s.first);
    
    for (FESystemUInt i=0; i<s.first; i++)
    {
        deriv.zero();
        for (FESystemUInt j=0; j<s.first; j++)  deriv_order[j]=0; // zero out all derivative orders before setting the value
        deriv_order[i]=1; // first derivative of the i^th dimension
        this->getFunctionDerivative(deriv_order, vin, deriv);
        mat.setColumnVals(i,0,s.second-1,deriv);
    }
}


template <typename ValType>
void
FESystem::Functions::LegendreFunctionMapping<ValType>::getFunctionDerivative(const std::vector<FESystemUInt>& derivative_orders, const FESystem::Numerics::VectorBase<ValType>& vin, FESystem::Numerics::VectorBase<ValType>& deriv) const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    FESystemAssert0(this->if_discrete_function_initialized, FESystem::Exception::InvalidState);
    
    // make sure that the dimensions are consistent
    std::pair<FESystemUInt, FESystemUInt> s=this->getDimensions(), s_discrete=this->discrete_function_values->getSize();
    FESystemAssert2(vin.getSize() == s.first, FESystem::Exception::DimensionsDoNotMatch, vin.getSize(), s.first);
    FESystemAssert2(deriv.getSize() == s.second, FESystem::Exception::DimensionsDoNotMatch, deriv.getSize(), s.second);
    
    FESystemUInt n_pts = this->getNDiscreteFunctions();
    
    FESystemAssert4(((s_discrete.first==s.second)&&(s_discrete.second == n_pts)), FESystem::Numerics::MatrixSizeMismatch, s_discrete.first, s_discrete.second, s.second, n_pts);
    
    FESystem::Numerics::LocalVector<ValType> vec;
    vec.resize(n_pts);
    
    this->getDiscreteFunctionDerivative(derivative_orders, vin, vec);
    
    // now that the function values are available, calculate the interpolated quantity
    ValType out_val = 0.0;
    for (FESystemUInt i_dim=0; i_dim<s.second; i_dim++)
    {
        out_val = 0.0;
        for (FESystemUInt i_pts=0; i_pts<n_pts; i_pts++)
            out_val += vec.getVal(i_pts)*this->discrete_function_values->getVal(i_dim,i_pts);
        deriv.setVal(i_dim, out_val);
    }
}




template <typename ValType>
void
FESystem::Functions::LegendreFunctionMapping<ValType>::getDiscreteFunction(const FESystem::Numerics::VectorBase<ValType>& vin, FESystem::Numerics::VectorBase<ValType>& vout) const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    FESystemUInt n_pts = this->getNDiscreteFunctions();
    
    // make sure that the dimensions are consistent
    std::pair<FESystemUInt, FESystemUInt> s=this->getDimensions();
    FESystemAssert2(vin.getSize() == s.first, FESystem::Exception::DimensionsDoNotMatch, vin.getSize(), s.first);
    FESystemAssert2(vout.getSize() == n_pts, FESystem::Exception::DimensionsDoNotMatch, vout.getSize(), n_pts);
    
    // loops to calculate the shape functions at each point
    std::vector<FESystemUInt> dim_order_to_multiply(s.second);
    std::fill(dim_order_to_multiply.begin(), dim_order_to_multiply.end(), 0);
    
    FESystem::Functions::LegendreFunctionMapping<ValType>::OrderCounter order_for_dim(this->dimension, this->order);
    
    ValType out_val = 0.0;
    for (FESystemUInt i_pts=0; i_pts<n_pts; i_pts++)
    {
        out_val = 1.0;
        for (FESystemUInt i_dim=0; i_dim<s.first; i_dim++)
            out_val *= this->legendre_function->getFunctionValue(vin.getVal(i_dim), order_for_dim.getOrderForDim(i_dim));
        order_for_dim.increment();
        vout.setVal(i_pts, out_val);
    }
}


template <typename ValType>
void
FESystem::Functions::LegendreFunctionMapping<ValType>::getDiscreteFunctionDerivative(const std::vector<FESystemUInt>& derivative_orders, const FESystem::Numerics::VectorBase<ValType>& vin, FESystem::Numerics::VectorBase<ValType>& deriv) const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    FESystemUInt n_pts = this->getNDiscreteFunctions();
    
    // make sure that the dimensions are consistent
    std::pair<FESystemUInt, FESystemUInt> s=this->getDimensions();
    FESystemAssert2(vin.getSize() == s.first, FESystem::Exception::DimensionsDoNotMatch, vin.getSize(), s.first);
    FESystemAssert2(deriv.getSize() == n_pts, FESystem::Exception::DimensionsDoNotMatch, deriv.getSize(), n_pts);
    
    // iterate over each dimension and shape function and calculate the mapped function value
    std::vector<std::vector<ValType> > func_val_deriv;
    if (func_val_deriv.size() != s.first)
        func_val_deriv.resize(s.first);
    for (FESystemUInt i=0; i<s.first; i++)
    {
        func_val_deriv[i].resize(n_pts);
        for (FESystemUInt j=0; j<n_pts; j++) func_val_deriv[i][j]=1.0;
    }
    
    FESystem::Functions::LegendreFunctionMapping<ValType>::OrderCounter order_for_dim(this->dimension, this->order);
    
    ValType out_val = 0.0;
    for (FESystemUInt i_pts=0; i_pts<n_pts; i_pts++)
    {
        out_val = 1.0;
        for (FESystemUInt i_dim=0; i_dim<s.first; i_dim++)
            out_val *= this->legendre_function->getFunctionDerivative(vin.getVal(i_dim), order_for_dim.getOrderForDim(i_dim), derivative_orders[i_dim]);
        order_for_dim.increment();
        deriv.setVal(i_pts, out_val);
    }
}



template <typename ValType>
FESystem::Functions::LegendreFunctionMapping<ValType>::OrderCounter::OrderCounter(const FESystemUInt dim, const FESystemUInt o):
dimension(dim),
order(o)
{
    this->order_per_dim.resize(dim);
    std::fill(this->order_per_dim.begin(), this->order_per_dim.end(), 0);
}


template <typename ValType>
void
FESystem::Functions::LegendreFunctionMapping<ValType>::OrderCounter::increment()
{
    // increment the last counter, and then check for limits
    FESystemInt i_dim=this->dimension-1;
    while (i_dim >= 0)
    {
        this->order_per_dim[i_dim]++;
        if (this->order_per_dim[i_dim] <= this->order)
            break;
        else
            this->order_per_dim[i_dim] = 0;
        i_dim--;
    }
}


template <typename ValType>
void
FESystem::Functions::LegendreFunctionMapping<ValType>::OrderCounter::write(std::ostream& out) const
{
    for (FESystemUInt i_dim=0; i_dim<this->dimension; i_dim++)
        out << this->order_per_dim[i_dim] << "  ";
    out << std::endl;
}


template <typename ValType>
FESystemUInt
FESystem::Functions::LegendreFunctionMapping<ValType>::OrderCounter::getOrderForDim(const FESystemUInt dim) const
{
    FESystemAssert2(dim < this->dimension, FESystem::Exception::IndexOutOfBound, dim, this->dimension);
    return this->order_per_dim[dim];
}




/***************************************************************************************/
// Template instantiations for some generic classes

INSTANTIATE_CLASS_FOR_ONLY_REAL_DATA_TYPES(FESystem::Functions::LegendreFunctionMapping);


/***************************************************************************************/

