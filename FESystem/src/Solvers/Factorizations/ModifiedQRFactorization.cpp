//
//  ModifiedQRFactorization.cpp
//  FESystem
//
//  Created by Manav Bhatia on 4/8/11.
//  Copyright 2011 . All rights reserved.
//

// C++ includes
#include <vector>

// FESystem includes
#include "Solvers/Factorizations/ModifiedQRFactorization.h"
#include "Solvers/Factorizations/TriangularBacksubstitution.h"
#include "Numerics/MatrixBase.h"
#include "Numerics/VectorBase.h"
#include "Base/macros.h"



template <typename ValType>
FESystem::FactorizationSolvers::ModifiedQRFactorization<ValType>::ModifiedQRFactorization():
FESystem::FactorizationSolvers::MatrixQRFactorizationBase<ValType>()
{
    
}


template <typename ValType>
FESystem::FactorizationSolvers::ModifiedQRFactorization<ValType>::~ModifiedQRFactorization()
{
    
}



template <typename ValType>
void 
FESystem::FactorizationSolvers::ModifiedQRFactorization<ValType>::factorize()
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
    tmat3->resize(s.first, s.second);
    
    std::vector<FESystemUInt> null_vector_indices;
    
    std::auto_ptr<FESystem::Numerics::VectorBase<ValType> > 
    vec1(FESystem::Numerics::VectorCreate<ValType>(FESystem::Numerics::LOCAL_VECTOR)),
    vec2(FESystem::Numerics::VectorCreate<ValType>(FESystem::Numerics::LOCAL_VECTOR));
    
    vec1->resize(s.first);
    vec2->resize(s.first);
    
    this->R_mat->setToIdentity();
    this->Q_mat->copyMatrix(this->getMatrix());
    
    FESystemDouble v;
    ValType val;
    
    for (FESystemUInt i=0; i<s.second; i++)
    {
        // this vector is used the base for this operation
        this->Q_mat->getColumnVals(i, 0, s.first-1, *vec1);
        
        // set the multiplying matrix 
        tmat1->setToIdentity();
        v = vec1->getL2Norm();
        
        if (v > MACHINE_EPSILON)
        {
            vec1->scale(1.0/v);
            this->Q_mat->setColumnVals(i, 0, s.first-1, *vec1);
            
            for (FESystemUInt j=i+1; j<s.second; j++)
                tmat1->setVal(i, j, -(this->Q_mat->dotProductWithColumn(j, *vec1)));

            // multiply this matrix to the A matrix and update the product of the R matrix
            this->Q_mat->matrixRightMultiply(1.0, *tmat1, *tmat2);
            this->Q_mat->copyMatrix(*tmat2);
            
            tmat1->scaleRow(i, 1.0/v);
            
            this->R_mat->matrixRightMultiply(1.0, *tmat1, *tmat2);
            this->R_mat->copyMatrix(*tmat2);
            tmat1->zero();
        }
        else
        {
            null_vector_indices.push_back(i);
            
            // create a random vector and orthogonalize it with respect to the previous bases
            vec1->initializeToRandomUnitVector(-1.0, 1.0);
            
            // orthogonalize with respect to previous vectors
            for (FESystemUInt j=0; j<i; j++)
            {
                this->Q_mat->getColumnVals(j, 0, s.first-1, *vec2);
                val = -vec1->dotProduct(*vec2);
                val /= vec2->getL2Norm();
                vec1->add(val, *vec2);
            }
            // orthogonalize with respect to the null space vectors
            for (FESystemUInt j=0; j<null_vector_indices.size()-1; j++)
            {
                tmat3->getColumnVals(j, 0, s.first-1, *vec2);
                val = -vec1->dotProduct(*vec2);
                val /= vec2->getL2Norm();
                vec1->add(val, *vec2);
            }
            
            vec1->scaleToUnitLength();
            
            // store in the matrix
            tmat3->setColumnVals(null_vector_indices.size()-1, 0, s.first-1, *vec1);            
        }
        
    }
    
    // invert the R factor matrix to get the R matrix for QR factorization
    FESystem::FactorizationSolvers::TriangularBacksubstitution<ValType> tri_bs;
    tri_bs.setTriangularMatrixType(FESystem::FactorizationSolvers::UPPER_TRIANGULAR);
    tri_bs.setMatrix(*(this->R_mat));
    
    tmat1->setToIdentity();
    tri_bs.backSubstitute(*tmat1,*tmat2);
    this->R_mat->copyMatrix(*tmat2);
    
    FESystemUInt id;
    for (FESystemUInt i=0; i < null_vector_indices.size(); i++)
    {
        id = null_vector_indices[i];
        this->R_mat->setVal(id, id, 0.0);
        tmat3->getColumnVals(i, 0, s.first-1, *vec1);
        this->Q_mat->setColumnVals(id, 0, s.first-1, *vec1);
    }

    this->factorization_complete = true;
}



/***************************************************************************************/
// Template instantiations for some generic classes

INSTANTIATE_CLASS_FOR_ALL_DATA_TYPES(FESystem::FactorizationSolvers::ModifiedQRFactorization);

/***************************************************************************************/


