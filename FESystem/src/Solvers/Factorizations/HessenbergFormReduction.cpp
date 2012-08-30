//
//  HessenbergFormReduction.cpp
//  FESystem
//
//  Created by Manav Bhatia on 5/8/11.
//  Copyright 2011 . All rights reserved.
//




// FESystem includes
#include "Solvers/Factorizations/HessenbergFormReduction.h"
#include "Base/FESystemExceptions.h"
#include "Numerics/MatrixBase.h"
#include "Numerics/VectorBase.h"
#include "Solvers/Factorizations/HouseHolderTriangulation.h"
#include "Base/macros.h"


template <typename ValType>
FESystem::FactorizationSolvers::HessenbergFormReduction<ValType>::HessenbergFormReduction():
mat(NULL),
factorization_complete(false)
{
    
}

template <typename ValType>
FESystem::FactorizationSolvers::HessenbergFormReduction<ValType>::~HessenbergFormReduction()
{
    
}



template <typename ValType>
void
FESystem::FactorizationSolvers::HessenbergFormReduction<ValType>::initializeMatrices()
{
    
    FESystemAssert0(this->mat != NULL, FESystem::Exception::NULLQuantity);
    
    this->Q_mat.reset(FESystem::Numerics::MatrixCreate<ValType>(this->mat->getType()).release());
    this->H_mat.reset(FESystem::Numerics::MatrixCreate<ValType>(this->mat->getType()).release());
    
    std::pair<FESystemUInt, FESystemUInt> s = this->mat->getSize();
    
    this->Q_mat->resize(s.first, s.second);
    this->H_mat->resize(s.first, s.second);    
}




template <typename ValType>
void
FESystem::FactorizationSolvers::HessenbergFormReduction<ValType>::setMatrix(FESystem::Numerics::MatrixBase<ValType>* m)
{
    FESystemAssert0(m != NULL, FESystem::Exception::NULLQuantity);
    FESystemAssert0(this->mat == NULL, FESystem::Exception::InvalidState);
    
    this->mat = m;    
}



template <typename ValType>
FESystem::Numerics::MatrixBase<ValType>& 
FESystem::FactorizationSolvers::HessenbergFormReduction<ValType>::getMatrix()
{
    FESystemAssert0(this->mat != NULL, FESystem::Exception::NULLQuantity);
    
    return *(this->mat);
}


template <typename ValType>
FESystem::Numerics::MatrixBase<ValType>& 
FESystem::FactorizationSolvers::HessenbergFormReduction<ValType>::getQMatrix()
{
    FESystemAssert0(this->Q_mat.get() != NULL, FESystem::Exception::NULLQuantity);
    FESystemAssert0(this->factorization_complete, FESystem::Exception::InvalidState);
    
    return *(this->Q_mat);
}



template <typename ValType>
FESystem::Numerics::MatrixBase<ValType>& 
FESystem::FactorizationSolvers::HessenbergFormReduction<ValType>::getHMatrix()
{
    FESystemAssert0(this->H_mat.get() != NULL, FESystem::Exception::NULLQuantity);
    FESystemAssert0(this->factorization_complete, FESystem::Exception::InvalidState);
    
    return *(this->H_mat);
}



template <typename ValType>
void 
FESystem::FactorizationSolvers::HessenbergFormReduction<ValType>::factorize()
{
    FESystemAssert0(!this->factorization_complete, FESystem::Exception::InvalidState);
        
    std::pair<FESystemUInt, FESystemUInt> s = this->mat->getSize();
    
    this->initializeMatrices();
    
    std::auto_ptr<FESystem::Numerics::MatrixBase<ValType> > 
    tmat1(FESystem::Numerics::MatrixCreate<ValType>(this->mat->getType())),
    tmat2(FESystem::Numerics::MatrixCreate<ValType>(this->mat->getType())),
    tmat3(FESystem::Numerics::MatrixCreate<ValType>(this->mat->getType()));
    
    tmat2->resize(s.first, s.second);
    tmat3->resize(s.first, s.second);
        
    std::auto_ptr<FESystem::FactorizationSolvers::MatrixQRFactorizationBase<ValType> > 
    qr_factorization(new FESystem::FactorizationSolvers::HouseholderTriangulation<ValType>());

    this->Q_mat->setToIdentity();
    this->H_mat->copyMatrix(*(this->mat));
    
    for (FESystemUInt i=0; i<s.first-1; i++)
    {
        tmat1->resize(s.first-1-i,s.second-1-i);
        
        this->H_mat->getSubMatrixVals(i+1, s.first-1, // first and last row id of the block
                                      i, s.second-2, // first and last column id of the block
                                      0,s.first-2-i, // all rows of tmat1
                                      0,s.second-2-i, // all columns of tmat1
                                      *(tmat1));
        
        qr_factorization->clear();
        qr_factorization->setMatrix(tmat1.get());
        qr_factorization->factorize();
        
        tmat2->setToIdentity();
        tmat2->setSubMatrixVals(i+1, s.first-1, // first and last row id of the block
                                i+1, s.second-1, // first and last column id of the block
                                0,s.first-2-i, // all rows of tmat1
                                0,s.second-2-i, // all columns of tmat1
                                qr_factorization->getQMatrix());
        
        tmat2->matrixTransposeRightMultiply(1.0, *(this->H_mat), *tmat3); // T3 = Q^* H
        tmat3->matrixRightMultiply(1.0, *tmat2, *(this->H_mat));  // H = T3  Q
        
        tmat3->zero();
        this->Q_mat->matrixRightMultiply(1.0, *tmat2, *tmat3);
        this->Q_mat->copyMatrix(*tmat3);
    }
    
    this->factorization_complete = true;

}



/***************************************************************************************/
// Template instantiations for some generic classes

INSTANTIATE_CLASS_FOR_ALL_DATA_TYPES(FESystem::FactorizationSolvers::HessenbergFormReduction);

/***************************************************************************************/


