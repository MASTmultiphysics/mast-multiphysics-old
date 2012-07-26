//
//  MatrixQRFactorizationBase.cpp
//  FESystem
//
//  Created by Manav Bhatia on 4/8/11.
//  Copyright 2011 . All rights reserved.
//


// FESystem includes
#include "Solvers/Factorizations/MatrixQRFactorizationBase.h"
#include "Numerics/MatrixBase.h"
#include "Numerics/VectorBase.h"
#include "Base/macros.h"



template <typename ValType>
FESystem::Solvers::MatrixQRFactorizationBase<ValType>::MatrixQRFactorizationBase():
factorization_complete(false),
mat(NULL)
{
    
}


template <typename ValType>
FESystem::Solvers::MatrixQRFactorizationBase<ValType>::~MatrixQRFactorizationBase()
{
    
}



template <typename ValType>
void
FESystem::Solvers::MatrixQRFactorizationBase<ValType>::clear()
{
    this->mat = NULL;
    
    this->Q_mat.reset();
    this->R_mat.reset();
    
    this->factorization_complete = false;
}



template <typename ValType>
void
FESystem::Solvers::MatrixQRFactorizationBase<ValType>::initializeMatrices()
{
    FESystemAssert0(this->mat != NULL, FESystem::Solvers::MatrixNotSetBeforeFactorization);
    
    this->Q_mat.reset(FESystem::Numerics::MatrixCreate<ValType>(this->mat->getType()).release());
    this->R_mat.reset(FESystem::Numerics::MatrixCreate<ValType>(this->mat->getType()).release());
    
    std::pair<FESystemUInt, FESystemUInt> s = this->mat->getSize();
    
    this->Q_mat->resize(s.first, s.second);
    this->R_mat->resize(s.first, s.second);    
}




template <typename ValType>
void
FESystem::Solvers::MatrixQRFactorizationBase<ValType>::setMatrix(const FESystem::Numerics::MatrixBase<ValType>* m)
{
    FESystemAssert0(m != NULL, FESystem::Exception::NULLQuantity);
    FESystemAssert0(this->mat == NULL, FESystem::Exception::InvalidState);
    FESystemAssert0(!this->factorization_complete, FESystem::Exception::InvalidState);
    
    this->mat = m;
}


template <typename ValType>
const FESystem::Numerics::MatrixBase<ValType>&
FESystem::Solvers::MatrixQRFactorizationBase<ValType>::getMatrix() const
{
    FESystemAssert0(this->mat != NULL, FESystem::Exception::NULLQuantity);
    
    return *(this->mat);
}


template <typename ValType>
FESystem::Numerics::MatrixBase<ValType>&
FESystem::Solvers::MatrixQRFactorizationBase<ValType>::getQMatrix()
{
    FESystemAssert0(this->factorization_complete, FESystem::Solvers::FactorizationNotComplete);
    FESystemAssert0(this->Q_mat.get() != NULL, FESystem::Exception::NULLQuantity);
    
    return *(this->Q_mat);
}


template <typename ValType>
FESystem::Numerics::MatrixBase<ValType>&
FESystem::Solvers::MatrixQRFactorizationBase<ValType>::getRMatrix()
{
    FESystemAssert0(this->factorization_complete, FESystem::Solvers::FactorizationNotComplete);
    FESystemAssert0(this->R_mat.get() != NULL, FESystem::Exception::NULLQuantity);
    
    return *(this->R_mat);
}


/***************************************************************************************/
// Template instantiations for some generic classes

INSTANTIATE_CLASS_FOR_ALL_DATA_TYPES(FESystem::Solvers::MatrixQRFactorizationBase);


/***************************************************************************************/


