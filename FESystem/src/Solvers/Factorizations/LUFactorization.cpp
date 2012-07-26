//
//  LUFactorization.cpp
//  FESystem
//
//  Created by Manav Bhatia on 6/17/12.
//  Copyright (c) 2012. All rights reserved.
//


// FESystem includes
#include "Solvers/Factorizations/LUFactorization.h"
#include "Numerics/MatrixBase.h"
#include "Numerics/VectorBase.h"
#include "Numerics/SparseMatrix.h"
#include "Numerics/SparsityPattern.h"
#include "Base/FESystemExceptions.h"
#include "Base/macros.h"



template <typename ValType>
FESystem::Solvers::LUFactorization<ValType>::LUFactorization():
factorization_complete(false),
mat(NULL)
{
    
}


template <typename ValType>
FESystem::Solvers::LUFactorization<ValType>::~LUFactorization()
{
    
}



template <typename ValType>
void
FESystem::Solvers::LUFactorization<ValType>::clear()
{
    this->mat = NULL;
    
    this->l_mat.reset();
    this->u_mat.reset();
    this->l_sparsity_pattern.reset();
    this->u_sparsity_pattern.reset();
    
    this->factorization_complete = false;
}



template <typename ValType>
void
FESystem::Solvers::LUFactorization<ValType>::initializeMatrices()
{
    FESystemAssert0(this->mat != NULL, FESystem::Exception::InvalidState);
    
    this->l_mat.reset(FESystem::Numerics::MatrixCreate<ValType>(this->mat->getType()).release());
    this->u_mat.reset(FESystem::Numerics::MatrixCreate<ValType>(this->mat->getType()).release());

    std::pair<FESystemUInt, FESystemUInt> s = this->mat->getSize();

    switch (this->mat->getType()) {
        case FESystem::Numerics::LOCAL_DENSE_MATRIX:
        {            
            this->l_mat->resize(s.first, s.second);
            this->u_mat->resize(s.first, s.second);    
        }
            break;

        case FESystem::Numerics::LOCAL_SPARSE_MATRIX:
        {
            // initialize the sparsity pattern
            this->l_sparsity_pattern.reset(new FESystem::Numerics::SparsityPattern);
            this->u_sparsity_pattern.reset(new FESystem::Numerics::SparsityPattern);

            // ask the sparsity pattern to initialize the sparsity pattern for the LU matrices
            dynamic_cast<const FESystem::Numerics::SparseMatrix<ValType>* >(this->mat)->getSparsityPattern().initSparsityPatternForLUFactorization(*(this->l_sparsity_pattern), *(this->u_sparsity_pattern));
            
            dynamic_cast<FESystem::Numerics::SparseMatrix<ValType>* >(this->l_mat.get())->resize(*this->l_sparsity_pattern);
            dynamic_cast<FESystem::Numerics::SparseMatrix<ValType>* >(this->u_mat.get())->resize(*this->u_sparsity_pattern);
        }
            break;
            
        default:
            FESystemAssert1(false, FESystem::Exception::EnumerationNotHandled, this->mat->getType());
            break;
    }

    this->mat->initializeLUFactoredMatrices(*(this->l_mat), *(this->u_mat));
    
    this->factorization_complete = true;
}




template <typename ValType>
void
FESystem::Solvers::LUFactorization<ValType>::setMatrix(const FESystem::Numerics::MatrixBase<ValType>* m)
{
    FESystemAssert0(m != NULL, FESystem::Exception::NULLQuantity);
    FESystemAssert0(this->mat == NULL, FESystem::Exception::InvalidState);
    FESystemAssert0(!this->factorization_complete, FESystem::Exception::InvalidState);
    
    this->mat = m;
    
    this->initializeMatrices();
}


template <typename ValType>
const FESystem::Numerics::MatrixBase<ValType>&
FESystem::Solvers::LUFactorization<ValType>::getMatrix() const
{
    FESystemAssert0(this->mat != NULL, FESystem::Exception::NULLQuantity);
    
    return *(this->mat);
}


template <typename ValType>
const FESystem::Numerics::MatrixBase<ValType>&
FESystem::Solvers::LUFactorization<ValType>::getLMatrix() const
{
    FESystemAssert0(this->factorization_complete, FESystem::Exception::InvalidState);
    FESystemAssert0(this->l_mat.get() != NULL, FESystem::Exception::NULLQuantity);
    
    return *(this->l_mat);
}


template <typename ValType>
const FESystem::Numerics::MatrixBase<ValType>&
FESystem::Solvers::LUFactorization<ValType>::getUMatrix() const
{
    FESystemAssert0(this->factorization_complete, FESystem::Exception::InvalidState);
    FESystemAssert0(this->u_mat.get() != NULL, FESystem::Exception::NULLQuantity);
    
    return *(this->u_mat);
}


/***************************************************************************************/
// Template instantiations for some generic classes

INSTANTIATE_CLASS_FOR_ALL_DATA_TYPES(FESystem::Solvers::LUFactorization);


/***************************************************************************************/


