//
//  LagrangeFunctionMapping.cpp
//  FESystem
//
//  Created by Manav Bhatia on 4/14/12.
//  Copyright (c) 2012. All rights reserved.
//

// FESystem includes
#include "Functions/LagrangeFunctionMapping.h"
#include "Functions/LagrangeFunction.h"
#include "Numerics/LocalVector.h"
#include "Numerics/MatrixBase.h"
#include "Utils/Table.h"
#include "Base/macros.h"


template <typename ValType>
FESystem::Functions::LagrangeFunctionMapping<ValType>::LagrangeFunctionMapping():
FESystem::Functions::DiscreteFunctionMappingBase<ValType>(),
dimension(0),
if_initialized(false)
{
    
}

template <typename ValType>
FESystem::Functions::LagrangeFunctionMapping<ValType>::~LagrangeFunctionMapping()
{
    
}

template <typename ValType>
void 
FESystem::Functions::LagrangeFunctionMapping<ValType>::clear()
{
    for (FESystemUInt i=0; i<this->lagrange_functions.size(); i++)
        this->lagrange_functions[i]=NULL;
    this->dimension=0;
    this->point_id_table = NULL;
    this->if_initialized = false;
    FESystem::Functions::DiscreteFunctionMappingBase<ValType>::clear();
}

template <typename ValType>
void
FESystem::Functions::LagrangeFunctionMapping<ValType>::reinit(FESystemUInt dim, const std::vector<FESystem::Functions::LagrangeFunction<ValType>*>& funcs,
                                                              const FESystem::Utility::Table<FESystemUInt>& pt_table)
{
    FESystemAssert0(!this->if_initialized, FESystem::Exception::InvalidState);
    
    // make sure that the provided values are correct
    FESystemAssert0(dim > 0, FESystem::Exception::InvalidValue);
    FESystemAssert0(funcs.size() > 0, FESystem::Exception::InvalidValue);

    // the dimension
    this->dimension = dim;
    if (this->lagrange_functions.size() != dim)
        this->lagrange_functions.resize(dim);
    for (FESystemUInt i=0; i<dim; i++)
        this->lagrange_functions[i] = funcs[i];

    // calculate the number of points
    FESystemUInt n=1;
    for (FESystemUInt i=0; i<dim; i++)
        n *= this->lagrange_functions[i]->getNPoints();

    // make sure that the dimension of the table is accurate
    FESystemAssert2(pt_table.getDimension() == 2, FESystem::Exception::DimensionsDoNotMatch, pt_table.getDimension(), 2);
    FESystemAssert2(pt_table.getNElementsInDimension(0) == n, FESystem::Exception::DimensionsDoNotMatch, pt_table.getNElementsInDimension(0), n);
    FESystemAssert2(pt_table.getNElementsInDimension(1) == dim, FESystem::Exception::DimensionsDoNotMatch, pt_table.getNElementsInDimension(1), dim);
    this->point_id_table = &pt_table;
    
    FESystem::Functions::DiscreteFunctionMappingBase<ValType>::reinit(dim, n);
    
    this->if_initialized = true;
}
            

template <typename ValType>
std::pair<FESystemUInt, FESystemUInt>
FESystem::Functions::LagrangeFunctionMapping<ValType>::getDimensions() const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    return std::pair<FESystemUInt, FESystemUInt> (this->lagrange_functions.size(), this->dimension);
}


template <typename ValType>
FESystemUInt 
FESystem::Functions::LagrangeFunctionMapping<ValType>::getNDiscreteFunctions() const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);

    FESystemUInt n=1;
    for (FESystemUInt i=0; i<this->lagrange_functions.size(); i++)
        n *= this->lagrange_functions[i]->getNPoints();
    return n;
}
            

template <typename ValType>
void 
FESystem::Functions::LagrangeFunctionMapping<ValType>::map(const FESystem::Numerics::VectorBase<ValType>& vin, FESystem::Numerics::VectorBase<ValType>& vout) const
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
FESystem::Functions::LagrangeFunctionMapping<ValType>::getMappingJacobian(const FESystem::Numerics::VectorBase<ValType>& vin, FESystem::Numerics::MatrixBase<ValType>& mat) const
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
FESystem::Functions::LagrangeFunctionMapping<ValType>::getFunctionDerivative(const std::vector<FESystemUInt>& derivative_orders, const FESystem::Numerics::VectorBase<ValType>& vin, FESystem::Numerics::VectorBase<ValType>& deriv) const
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
FESystem::Functions::LagrangeFunctionMapping<ValType>::getDiscreteFunction(const FESystem::Numerics::VectorBase<ValType>& vin, FESystem::Numerics::VectorBase<ValType>& vout) const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    FESystemUInt n_pts = this->getNDiscreteFunctions();

    // make sure that the dimensions are consistent
    std::pair<FESystemUInt, FESystemUInt> s=this->getDimensions();
    FESystemAssert2(vin.getSize() == s.first, FESystem::Exception::DimensionsDoNotMatch, vin.getSize(), s.first);
    FESystemAssert2(vout.getSize() == n_pts, FESystem::Exception::DimensionsDoNotMatch, vout.getSize(), n_pts);
    
    // iterate over each dimension and shape function and calculate the mapped function value
    std::vector<ValType> func_vals;
    if (func_vals.size() != n_pts)
        func_vals.resize(n_pts);
    for (FESystemUInt i=0; i<n_pts; i++) func_vals[i]=1.0;
    
    // loops to calculate the shape functions at each point
    ValType out_val = 0.0;
    for (FESystemUInt i_pts=0; i_pts<n_pts; i_pts++)
    {
        out_val = 1.0;
        for (FESystemUInt i_dim=0; i_dim<s.first; i_dim++)
            out_val *= this->lagrange_functions[i_dim]->getFunctionValue(vin.getVal(i_dim), this->getDiscretePointNumberAlongLocalDimension(i_pts, i_dim));
        vout.setVal(i_pts, out_val);
    }
}


template <typename ValType>
void 
FESystem::Functions::LagrangeFunctionMapping<ValType>::getDiscreteFunctionDerivative(const std::vector<FESystemUInt>& derivative_orders, const FESystem::Numerics::VectorBase<ValType>& vin, FESystem::Numerics::VectorBase<ValType>& deriv) const
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
    
    //  calculate the function derivative
    ValType out_val = 0.0;
    for (FESystemUInt i_pts=0; i_pts<n_pts; i_pts++)
    {
        out_val = 1.0;
        for (FESystemUInt i_dim=0; i_dim<s.first; i_dim++)
            out_val *= this->lagrange_functions[i_dim]->getFunctionDerivative(vin.getVal(i_dim), this->getDiscretePointNumberAlongLocalDimension(i_pts, i_dim), derivative_orders[i_dim]); 
        deriv.setVal(i_pts, out_val);
    }
}



template <typename ValType>
FESystemUInt 
FESystem::Functions::LagrangeFunctionMapping<ValType>::getDiscretePointNumberAlongLocalDimension(FESystemUInt p_id, FESystemUInt dim_id) const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    std::vector<FESystemUInt> el(2);
    el[0] = p_id;
    el[1] = dim_id;
    return this->point_id_table->getVal(el);
}


/***************************************************************************************/
// Template instantiations for some generic classes

INSTANTIATE_CLASS_FOR_ALL_DATA_TYPES(FESystem::Functions::LagrangeFunctionMapping);


/***************************************************************************************/

