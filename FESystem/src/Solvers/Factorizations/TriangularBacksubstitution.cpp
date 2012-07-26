//
//  TriangularBacksubstitution.cpp
//  FESystem
//
//  Created by Manav Bhatia on 4/6/11.
//  Copyright 2011 . All rights reserved.
//


// FESystem includes
#include "Base/FESystemExceptions.h"
#include "Solvers/Factorizations/TriangularBacksubstitution.h"
#include "Numerics/MatrixBase.h"
#include "Numerics/VectorBase.h"
#include "Base/macros.h"


template <typename ValType> 
FESystem::Solvers::TriangularBacksubstitution<ValType>::TriangularBacksubstitution():
tri_type(FESystem::Solvers::INVALID_TRIANGULAR_PART),
mat(NULL)
{
    
}



template <typename ValType> 
FESystem::Solvers::TriangularBacksubstitution<ValType>::~TriangularBacksubstitution()
{
    
}



template <typename ValType> 
void
FESystem::Solvers::TriangularBacksubstitution<ValType>::clear()
{
    this->mat = NULL;
    this->tri_type = FESystem::Solvers::INVALID_TRIANGULAR_PART;
}




template <typename ValType> 
FESystem::Solvers::TriangularType
FESystem::Solvers::TriangularBacksubstitution<ValType>::getTriangularMatrixType()
{
    return this->tri_type;
}



template <typename ValType> 
void 
FESystem::Solvers::TriangularBacksubstitution<ValType>::setTriangularMatrixType(FESystem::Solvers::TriangularType t)
{
    FESystemAssert0(t != FESystem::Solvers::INVALID_TRIANGULAR_PART,
                    FESystem::Exception::InvalidValue);
    
    this->tri_type = t;
}
    


template <typename ValType> 
void 
FESystem::Solvers::TriangularBacksubstitution<ValType>::setMatrix(const FESystem::Numerics::MatrixBase<ValType>& m)
{
    this->mat = &m;
}


template <typename ValType> 
const FESystem::Numerics::MatrixBase<ValType>&
FESystem::Solvers::TriangularBacksubstitution<ValType>::getMatrix() const
{
    FESystemAssert0(this->mat != NULL, FESystem::Exception::NULLQuantity);
    return *(this->mat);
}

    

template <typename ValType> 
void 
FESystem::Solvers::TriangularBacksubstitution<ValType>::backSubstitute(const FESystem::Numerics::VectorBase<ValType>& v1,
                                                                        FESystem::Numerics::VectorBase<ValType>& res) 
{
    FESystemAssert0(this->mat != NULL, FESystem::Exception::NULLQuantity);
    FESystemAssert0(this->tri_type != FESystem::Solvers::INVALID_TRIANGULAR_PART, FESystem::Exception::InvalidValue);
    
    std::pair<FESystemUInt, FESystemUInt> s = this->getMatrix().getSize();
    res.zero();
    ValType v = 0.0;
    FESystemUInt row;
    
    switch (this->getTriangularMatrixType())
    {
        case FESystem::Solvers::UPPER_TRIANGULAR:
        {
            for (FESystemUInt i=0; i < s.first; i++)
            {
                row = s.first - i - 1;
                v = v1.getVal(row);
                if (row < s.first-1)
                    v -= this->mat->multiplySubVectorWithSubRow(row, row+1, s.second-1, row+1, s.second-1, res);
                res.setVal(row, v/this->mat->getVal(row, row));
            }
        }
            break;
            
        case FESystem::Solvers::LOWER_TRIANGULAR:
        {
            for (FESystemUInt row=0; row < s.first; row++)
            {
                v = v1.getVal(row);
                if (row > 0)
                    v -= this->mat->multiplySubVectorWithSubRow(row, 0, row-1, 0, row-1, res);
                res.setVal(row, v/this->mat->getVal(row, row));
            }
        }
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
    }
}
    


template <typename ValType> 
void 
FESystem::Solvers::TriangularBacksubstitution<ValType>::backSubstitute(const FESystem::Numerics::MatrixBase<ValType>& m1,
                                                                        FESystem::Numerics::MatrixBase<ValType>& res) 
{
    FESystemAssert0(this->mat != NULL, FESystem::Exception::NULLQuantity);
    FESystemAssert0(this->tri_type != FESystem::Solvers::INVALID_TRIANGULAR_PART, FESystem::Exception::InvalidValue);
    
    std::pair<FESystemUInt, FESystemUInt> s = this->getMatrix().getSize();
    res.zero();

    std::auto_ptr<FESystem::Numerics::VectorBase<ValType> > vec1(FESystem::Numerics::VectorCreate<ValType>(FESystem::Numerics::LOCAL_VECTOR).release()), 
    vec2(FESystem::Numerics::VectorCreate<ValType>(FESystem::Numerics::LOCAL_VECTOR).release());
    vec1->resize(s.first); vec2->resize(s.first);
    
    for (FESystemUInt i=0; i<s.second; i++)
    {
        vec1->zero(); 
        m1.getColumnVals(i, 0, s.first-1, *vec1);
        this->backSubstitute(*vec1, *vec2);
        res.setColumnVals(i, 0, s.first-1, *vec2);
    }
}
    
    
/***************************************************************************************/
// Template instantiations for some generic classes

INSTANTIATE_CLASS_FOR_ALL_DATA_TYPES(FESystem::Solvers::TriangularBacksubstitution);

/***************************************************************************************/

