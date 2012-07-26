/*
 *  MatrixBase.cpp
 *  FESystem
 *
 *  Created by Manav Bhatia on 1/9/11.
 *  Copyright 2011 . All rights reserved.
 *
 */


// FESystem includes
#include "Numerics/MatrixBase.h"
#include "Numerics/DenseMatrix.h"
#include "Numerics/SparseMatrix.h"
#include "Base/macros.h"


template <typename ValType>
std::auto_ptr<FESystem::Numerics::MatrixBase<ValType> >
FESystem::Numerics::MatrixCreate(FESystem::Numerics::MatrixType mat_type)
{
    std::auto_ptr<FESystem::Numerics::MatrixBase<ValType> > mat;
    
    switch (mat_type) {
        case FESystem::Numerics::LOCAL_DENSE_MATRIX:
            mat.reset(new FESystem::Numerics::DenseMatrix<ValType>());
            break;
            
        case FESystem::Numerics::LOCAL_SPARSE_MATRIX:
            mat.reset(new FESystem::Numerics::SparseMatrix<ValType>());
            break;

        default:
            FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
    }
    
    return mat;
}



/***************************************************************************************/
// Template instantiations for some generic classes

template std::auto_ptr<FESystem::Numerics::MatrixBase<FESystemFloat> > FESystem::Numerics::MatrixCreate(FESystem::Numerics::MatrixType mat_type);
template std::auto_ptr<FESystem::Numerics::MatrixBase<FESystemDouble> > FESystem::Numerics::MatrixCreate(FESystem::Numerics::MatrixType mat_type);
template std::auto_ptr<FESystem::Numerics::MatrixBase<FESystemComplexFloat> > FESystem::Numerics::MatrixCreate(FESystem::Numerics::MatrixType mat_type);
template std::auto_ptr<FESystem::Numerics::MatrixBase<FESystemComplexDouble> > FESystem::Numerics::MatrixCreate(FESystem::Numerics::MatrixType mat_type);

/***************************************************************************************/

