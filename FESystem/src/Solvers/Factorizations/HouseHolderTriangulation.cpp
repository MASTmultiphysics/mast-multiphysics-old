//
//  HouseholderTriangulation.cpp
//  FESystem
//
//  Created by Manav Bhatia on 3/27/11.
//  Copyright 2011 . All rights reserved.
//

// FESystem includes
#include "Solvers/Factorizations/HouseHolderTriangulation.h"
#include "Solvers/Factorizations/TriangularBacksubstitution.h"
#include "Numerics/MatrixBase.h"
#include "Numerics/VectorBase.h"
#include "Base/macros.h"



template <typename ValType>
FESystem::FactorizationSolvers::HouseholderTriangulation<ValType>::HouseholderTriangulation():
FESystem::FactorizationSolvers::MatrixQRFactorizationBase<ValType>()
{
    
}


template <typename ValType>
FESystem::FactorizationSolvers::HouseholderTriangulation<ValType>::~HouseholderTriangulation()
{
    
}



template <typename ValType>
void 
FESystem::FactorizationSolvers::HouseholderTriangulation<ValType>::factorize()
{
    FESystemAssert0(!this->factorization_complete, FESystem::Exception::InvalidState);
    
    std::pair<FESystemUInt, FESystemUInt> s = this->mat->getSize();
    
    this->initializeMatrices();
    
    std::auto_ptr<FESystem::Numerics::MatrixBase<ValType> > 
    tmat1(FESystem::Numerics::MatrixCreate<ValType>(this->mat->getType())),
    tmat2(FESystem::Numerics::MatrixCreate<ValType>(this->mat->getType())),
    tmat3(FESystem::Numerics::MatrixCreate<ValType>(this->mat->getType()));
    tmat1->resize(s.first, s.second);
    tmat2->resize(s.first, s.second);
    
    std::auto_ptr<FESystem::Numerics::VectorBase<ValType> > 
    vec1(FESystem::Numerics::VectorCreate<ValType>(FESystem::Numerics::LOCAL_VECTOR));
        
    tmat2->setToIdentity();
    this->R_mat->copyMatrix(this->getMatrix());
    
    for (FESystemUInt i=0; i<s.second; i++)
    {
        // create orthogonal projector using the first column of the submatrix
        vec1->resize(s.first-i);
        tmat3->resize(s.first-i, s.first-i);
        tmat1->resize(s.first-i, s.second);
        this->R_mat->getColumnVals(i, i, s.first-1, *vec1);
        
        tmat3->createOrthogonalReflectorForVector(*vec1);
                
        // multiply this matrix to the A matrix and update the product of the R matrix
        tmat3->matrixRightMultiplySubMatrix(i, s.first-1, 0, s.second-1,
                                            1.0, *(this->R_mat), *tmat1);
        
        this->R_mat->setSubMatrixVals(i, s.first-1, 0, s.second-1, 
                                      0, s.first-i-1, 0, s.second-1,
                                      *tmat1);
        tmat3->matrixRightMultiplySubMatrix(i, s.first-1, 0, s.second-1,
                                            1.0, *tmat2, *tmat1);
        tmat2->setSubMatrixVals(i, s.first-1, 0, s.second-1, 
                                0, s.first-i-1, 0, s.second-1,
                                *tmat1);
        
        tmat1->zero();
    }

    // Q = Q_{HT}^H 
    this->Q_mat->copyMatrixTranspose(*tmat2);

    
    this->factorization_complete = true;
}



/***************************************************************************************/
// Template instantiations for some generic classes

INSTANTIATE_CLASS_FOR_ALL_DATA_TYPES(FESystem::FactorizationSolvers::HouseholderTriangulation);


/***************************************************************************************/


