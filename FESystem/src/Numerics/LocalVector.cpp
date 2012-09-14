
/*
 *  LocalVector.cpp
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/27/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

// C++ includes
#include <cmath>
#include <stdlib.h>


// FESystem includes
#include "Numerics/LocalVector.h"
#include "Base/macros.h"

template <typename ValType>
FESystem::Numerics::LocalVector<ValType>::LocalVector():
vec_vals(NULL),
n_vals(0),
if_owns_pointer(false)
{
    
}


template <typename ValType>
FESystem::Numerics::LocalVector<ValType>::LocalVector(FESystemUInt n):
vec_vals(NULL),
n_vals(0),
if_owns_pointer(false)
{
    this->resize(n);
}


template <typename ValType>
FESystem::Numerics::LocalVector<ValType>::~LocalVector()
{
    this->clear();
}


template <typename ValType>
void
FESystem::Numerics::LocalVector<ValType>::clear()
{
    if (this->if_owns_pointer)
        if (vec_vals != NULL)
            delete[] this->vec_vals;
    this->vec_vals = NULL;
    this->n_vals = 0;
    this->if_owns_pointer = false;
}


template <typename ValType>
FESystem::Numerics::VectorType
FESystem::Numerics::LocalVector<ValType>::getType() const
{
    return FESystem::Numerics::LOCAL_VECTOR;
}


template <typename ValType>
FESystemUInt
FESystem::Numerics::LocalVector<ValType>::getSize() const
{
    return this->n_vals;
} 


template <typename ValType>
void 
FESystem::Numerics::LocalVector<ValType>::resize(FESystemUInt i)
{
    if (this->if_owns_pointer)
    {
        if (this->getSize() != i)
        {
            delete[] this->vec_vals;
            this->vec_vals = new ValType[i];
        }
    }
    else // if it does not own the pointer, then reset the pointer to a data created locally
    {
        this->vec_vals = NULL;
        this->vec_vals = new ValType[i];
    }
    
    this->n_vals = i;
    this->if_owns_pointer = true;
    
    this->zero();
} 


template <typename ValType>
void 
FESystem::Numerics::LocalVector<ValType>::resize(FESystemUInt i, ValType* vals)
{
    this->clear();
    this->vec_vals = vals;
    this->n_vals = i;
    this->if_owns_pointer = false;
} 



template <typename ValType>
void 
FESystem::Numerics::LocalVector<ValType>::copyVector(const VectorBase<ValType>& t)
{
    this->resize(t.getSize());
    
    const ValType* t_val = t.getVectorValues();
    
    for (unsigned int i=0; i< this->getSize(); i++)
        this->vec_vals[i] = t_val[i];
}




template <typename ValType>
void 
FESystem::Numerics::LocalVector<ValType>::copyVectorVals(const VectorBase<ValType>& t)
{
    FESystemAssert2(t.getSize() == this->getSize(), 
                    FESystem::Numerics::VectorSizeMismatch,
                    this->getSize(), t.getSize()); 
    
    const ValType* t_val = t.getVectorValues();
    
    for (unsigned int i=0; i< this->getSize(); i++)
        this->vec_vals[i] = t_val[i];
}




template <typename ValType>
void 
FESystem::Numerics::LocalVector<ValType>::setSubVectorVals(FESystemUInt row1, FESystemUInt row2, 
                                                           FESystemUInt v_row1, FESystemUInt v_row2, 
                                                           const FESystem::Numerics::VectorBase<ValType>& v)
{
    FESystemAssert2( row1 < this->getSize(),
                    FESystem::Numerics::ElementLocationExceedesVectorSize,
                    this->getSize(), row1);
    FESystemAssert2( row2 < this->getSize(),
                    FESystem::Numerics::ElementLocationExceedesVectorSize,
                    this->getSize(), row2);
    
    FESystemAssert2( v_row1 < v.getSize(),
                    FESystem::Numerics::ElementLocationExceedesVectorSize,
                    v.getSize(), v_row1);
    FESystemAssert2( v_row2 < v.getSize(),
                    FESystem::Numerics::ElementLocationExceedesVectorSize,
                    v.getSize(), v_row2);
    
    FESystemAssert2((row2-row1) >= 0, 
                    FESystem::Exception::PositiveDifferenceNeeded,
                    row1, row2);
    FESystemAssert2((v_row2-v_row1) >= 0, 
                    FESystem::Exception::PositiveDifferenceNeeded,
                    v_row1, v_row2);
    FESystemAssert2((row2-row1) == (v_row2-v_row1),
                    FESystem::Numerics::VectorSizeMismatch,
                    (row2-row1), (v_row2-v_row1));
    
    const ValType* val = v.getVectorValues();
    
    for (FESystemUInt i=0; i<=(row2-row1); i++)
        this->vec_vals[row1+i] =  val[v_row1+i];

}



template <typename ValType>
void 
FESystem::Numerics::LocalVector<ValType>::addSubVectorVals(FESystemUInt row1, FESystemUInt row2, 
                                                           FESystemUInt v_row1, FESystemUInt v_row2, 
                                                           const ValType f, const FESystem::Numerics::VectorBase<ValType>& v)
{
    FESystemAssert2( row1 < this->getSize(),
                    FESystem::Numerics::ElementLocationExceedesVectorSize,
                    this->getSize(), row1);
    FESystemAssert2( row2 < this->getSize(),
                    FESystem::Numerics::ElementLocationExceedesVectorSize,
                    this->getSize(), row2);
    
    FESystemAssert2( v_row1 < v.getSize(),
                    FESystem::Numerics::ElementLocationExceedesVectorSize,
                    v.getSize(), v_row1);
    FESystemAssert2( v_row2 < v.getSize(),
                    FESystem::Numerics::ElementLocationExceedesVectorSize,
                    v.getSize(), v_row2);
    
    FESystemAssert2((row2-row1) >= 0, 
                    FESystem::Exception::PositiveDifferenceNeeded,
                    row1, row2);
    FESystemAssert2((v_row2-v_row1) >= 0, 
                    FESystem::Exception::PositiveDifferenceNeeded,
                    v_row1, v_row2);
    FESystemAssert2((row2-row1) == (v_row2-v_row1),
                    FESystem::Numerics::VectorSizeMismatch,
                    (row2-row1), (v_row2-v_row1));
    
    const ValType* val = v.getVectorValues();
    
    for (FESystemUInt i=0; i<=(row2-row1); i++)
        this->vec_vals[row1+i] += f* val[v_row1+i];
    
}



template <typename ValType>
void
FESystem::Numerics::LocalVector<ValType>::getSubVectorVals(FESystemUInt row1, FESystemUInt row2, 
                                                           FESystemUInt v_row1, FESystemUInt v_row2, 
                                                           FESystem::Numerics::VectorBase<ValType>& v) const
{
    FESystemAssert2( row1 < this->getSize(),
                    FESystem::Numerics::ElementLocationExceedesVectorSize,
                    this->getSize(), row1);
    FESystemAssert2( row2 < this->getSize(),
                    FESystem::Numerics::ElementLocationExceedesVectorSize,
                    this->getSize(), row2);
    
    FESystemAssert2( v_row1 < v.getSize(),
                    FESystem::Numerics::ElementLocationExceedesVectorSize,
                    v.getSize(), v_row1);
    FESystemAssert2( v_row2 < v.getSize(),
                    FESystem::Numerics::ElementLocationExceedesVectorSize,
                    v.getSize(), v_row2);
    
    FESystemAssert2((row2-row1) >= 0, 
                    FESystem::Exception::PositiveDifferenceNeeded,
                    row1, row2);
    FESystemAssert2((v_row2-v_row1) >= 0, 
                    FESystem::Exception::PositiveDifferenceNeeded,
                    v_row1, v_row2);
    FESystemAssert2((row2-row1) == (v_row2-v_row1),
                    FESystem::Numerics::VectorSizeMismatch,
                    (row2-row1), (v_row2-v_row1));
    
    ValType* val = v.getVectorValues();
    
    for (FESystemUInt i=0; i<=(row2-row1); i++)
        val[v_row1+i] = this->vec_vals[row1+i];
    
}



template <typename ValType>
void
FESystem::Numerics::LocalVector<ValType>::setSubVectorValsFromIndices(const std::vector<FESystemUInt>& rows, const FESystem::Numerics::VectorBase<ValType>& vec) 
{
    FESystemAssert2(rows.size() == vec.getSize(), FESystem::Exception::DimensionsDoNotMatch, vec.getSize(), rows.size());
    FESystemAssert2(rows[0] < this->getSize(), FESystem::Exception::IndexOutOfBound, rows[0], this->getSize());
    FESystemAssert2(*(rows.rbegin()) < this->getSize(), FESystem::Exception::IndexOutOfBound, *(rows.rbegin()), this->getSize());    
    
    const ValType* vals = vec.getVectorValues();
    
    // copy the values from the vector to the present subvector
    for (FESystemUInt i=0; i<rows.size(); i++)
        this->vec_vals[rows[i]] = vals[i];
}



template <typename ValType>
void
FESystem::Numerics::LocalVector<ValType>::getSubVectorValsFromIndices(const std::vector<FESystemUInt>& rows, FESystem::Numerics::VectorBase<ValType>& vec) const
{
    FESystemAssert2(rows.size() == vec.getSize(), FESystem::Exception::DimensionsDoNotMatch, vec.getSize(), rows.size());
    FESystemAssert2(rows[0] < this->getSize(), FESystem::Exception::IndexOutOfBound, rows[0], this->getSize());
    FESystemAssert2(*(rows.rbegin()) < this->getSize(), FESystem::Exception::IndexOutOfBound, *(rows.rbegin()), this->getSize());    
    
    ValType* vals = vec.getVectorValues();
    
    // copy the values from the vector to the present subvector
    for (FESystemUInt i=0; i<rows.size(); i++)
        vals[i] = this->vec_vals[rows[i]];
}



template <typename ValType>
ValType
FESystem::Numerics::LocalVector<ValType>::getMinVal() const
{
    ValType v_min;
    typename RealOperationType(ValType) v = FESystem::Base::getMachineMax<typename RealOperationType(ValType)>();
    for (FESystemUInt i=0; i<this->getSize(); i++)
        if ( FESystem::Base::comparisonValue<ValType, typename RealOperationType(ValType)>(this->vec_vals[i]) < v)
        {
            v_min = this->vec_vals[i];
            v = FESystem::Base::comparisonValue<ValType, typename RealOperationType(ValType)>(v_min);
        }
    
    return v_min;
}



template <typename ValType>
ValType
FESystem::Numerics::LocalVector<ValType>::getMaxVal() const
{
    ValType v_max;
    typename RealOperationType(ValType) v = FESystem::Base::getMachineMin<typename RealOperationType(ValType)>();
    for (FESystemUInt i=0; i<this->getSize(); i++)
        if ( FESystem::Base::comparisonValue<ValType, typename RealOperationType(ValType)>(this->vec_vals[i]) > v)
        {
            v_max = this->vec_vals[i];
            v = FESystem::Base::comparisonValue<ValType, typename RealOperationType(ValType)>(v_max);
        }
    
    return v_max;
}



template <typename ValType>
ValType
FESystem::Numerics::LocalVector<ValType>::getSum() const
{
    ValType v=0.0;
    for (FESystemUInt i=0; i<this->getSize(); i++)
        v += this->vec_vals[i];
    
    return v;
}



template <typename ValType>
const ValType& 
FESystem::Numerics::LocalVector<ValType>::getVal(const FESystemUInt i) const
{
    FESystemAssert2( i < this->getSize(),
                    FESystem::Numerics::ElementLocationExceedesVectorSize,
                    this->getSize(), i);
    
    return this->vec_vals[i];
} 			


template <typename ValType>
void 
FESystem::Numerics::LocalVector<ValType>::zero()
{
    for (unsigned int i=0; i< this->getSize(); i++)
        this->vec_vals[i] = 0.0;
}



template <typename ValType>
void 
FESystem::Numerics::LocalVector<ValType>::setVal(const FESystemUInt i, const ValType& val)
{
    FESystemAssert2( i < this->getSize(),
                    FESystem::Numerics::ElementLocationExceedesVectorSize,
                    this->getSize(), i);
    
    this->vec_vals[i] = val; 
} 


template <typename ValType>
void 
FESystem::Numerics::LocalVector<ValType>::addVal(const FESystemUInt i, const ValType& val)
{
    FESystemAssert2( i < this->getSize(),
                    FESystem::Numerics::ElementLocationExceedesVectorSize,
                    this->getSize(), i);
    
    this->vec_vals[i] += val; 
} 


template <typename ValType>
void 
FESystem::Numerics::LocalVector<ValType>::setAllVals(const ValType& val)
{    
    FESystemAssert0(this->getSize() > 0, FESystem::Exception::InvalidState);
    
    for (FESystemUInt i=0; i<this->getSize(); i++)
        this->vec_vals[i] = val; 
} 


template <typename ValType>
void 
FESystem::Numerics::LocalVector<ValType>::scale(const ValType& t)
{
    for (unsigned int i=0; i< this->getSize(); i++)
        this->vec_vals[i] *= t;
}


template <typename ValType>
void 
FESystem::Numerics::LocalVector<ValType>::scaleSubVector(FESystemUInt n1, FESystemUInt n2, const ValType &t)
{
    FESystemAssert2( n1 < this->getSize(),
                    FESystem::Numerics::ElementLocationExceedesVectorSize,
                    this->getSize(), n1);
    FESystemAssert2( n2 < this->getSize(),
                    FESystem::Numerics::ElementLocationExceedesVectorSize,
                    this->getSize(), n2);
    
    for (unsigned int i=n1; i<=n2; i++)
        this->vec_vals[i] *= t;
}


template <typename ValType>
void 
FESystem::Numerics::LocalVector<ValType>::scaleToUnitLength()
{
    this->scale(1.0 / this->getL2Norm());
}



template <typename ValType>
void
FESystem::Numerics::LocalVector<ValType>::initializeToRandomUnitVector(ValType low, ValType up)
{
    ValType val;
    for (FESystemUInt i=0; i < this->getSize(); i++)
    {
        val = (1.0 * rand()) / (1.0 * RAND_MAX);
        this->vec_vals[i] = val * (up - low) - low ;
    }
     
    this->scaleToUnitLength();
}




template <typename ValType>
typename RealOperationType(ValType)
FESystem::Numerics::LocalVector<ValType>::getL1Norm() const
{
    const ValType* vec_vals = this->getVectorValues();

    typename RealOperationType(ValType) val1 = 0.0, val2 = 0.0;
    
    for (FESystemUInt i=0; i<this->getSize(); i++)
    {
        val2 = FESystem::Base::magnitude<ValType, typename RealOperationType(ValType)> (vec_vals[i]);
        if (val2 > val1)
            val2 = val1;
    }

    return val1;
}



template <typename ValType>
typename RealOperationType(ValType)
FESystem::Numerics::LocalVector<ValType>::getL2Norm() const
{
    FESystemDouble v=0.0;
    const ValType* vec_vals = this->getVectorValues();
    
    for (FESystemUInt i=0; i<this->getSize(); i++)
        v += pow(FESystem::Base::magnitude<ValType, typename RealOperationType(ValType)>(vec_vals[i]), 2);
    
    return pow(v, 0.5);
}





template <typename ValType>
void 
FESystem::Numerics::LocalVector<ValType>::elementwiseMultiply(const FESystem::Numerics::VectorBase<ValType>& t)
{
    FESystemAssert2(t.getSize() == this->getSize(), 
                    FESystem::Numerics::VectorSizeMismatch,
                    this->getSize(), t.getSize()); 
    
    const ValType* val = t.getVectorValues();
    for (unsigned int i=0; i < this->getSize(); i++) 
        this->vec_vals[i] *= val[i];
}



template <typename ValType>
void 
FESystem::Numerics::LocalVector<ValType>::elementwiseMultiply
(const FESystem::Numerics::VectorBase<ValType>& t,
 FESystem::Numerics::VectorBase<ValType>& res)
{
    FESystemAssert2(t.getSize() == this->getSize(), 
                    FESystem::Numerics::VectorSizeMismatch,
                    this->getSize(), t.getSize()); 
    
    FESystemAssert2(res.getSize() == this->getSize(), 
                    FESystem::Numerics::VectorSizeMismatch,
                    this->getSize(), res.getSize()); 
    
    const ValType* val = t.getVectorValues();
    ValType* rval = res.getVectorValues();
    
    for (unsigned int i=0; i < this->getSize(); i++) 
        rval[i] =  this->vec_vals[i] * val[i];
}


template <typename ValType>
ValType
FESystem::Numerics::LocalVector<ValType>::dotProduct
(const FESystem::Numerics::VectorBase<ValType>& t) const
{
    FESystemAssert2(t.getSize() == this->getSize(), 
                    FESystem::Numerics::VectorSizeMismatch,
                    this->getSize(), t.getSize()); 
    
    const ValType* val = t.getVectorValues();
    
    ValType rval = 0.0;
    
    for (unsigned int i=0; i < this->getSize(); i++) 
        rval +=  this->vec_vals[i] * val[i];
    
    return rval;
}




template <typename ValType>
void
FESystem::Numerics::LocalVector<ValType>::crossProduct(const FESystem::Numerics::VectorBase<ValType>& v1,
                                                       FESystem::Numerics::VectorBase<ValType>& v2) const
{
    // this is allowed only for 3-D vectors
    FESystemAssert0(this->getSize() == 3, FESystem::Exception::InvalidFunctionCall); 

    FESystemAssert2(v1.getSize() == this->getSize(), FESystem::Numerics::VectorSizeMismatch, this->getSize(), v1.getSize()); 
    FESystemAssert2(v2.getSize() == this->getSize(), FESystem::Numerics::VectorSizeMismatch, this->getSize(), v2.getSize()); 
    
    const ValType *v0_val = this->getVectorValues(); const ValType *v1_val = v1.getVectorValues(); ValType *v2_val = v2.getVectorValues();
        
    v2_val[0] =  v0_val[1] * v1_val[2] - v0_val[2] * v1_val[1]; // x-value
    v2_val[1] = -v0_val[0] * v1_val[2] + v0_val[2] * v1_val[0]; // y-value
    v2_val[2] =  v0_val[0] * v1_val[1] - v0_val[1] * v1_val[0]; // z-value
}



template <typename ValType>
void
FESystem::Numerics::LocalVector<ValType>::add(ValType f, const FESystem::Numerics::VectorBase<ValType>& t)
{
    FESystemAssert2(t.getSize() == this->getSize(),
                    FESystem::Numerics::VectorSizeMismatch,
                    this->getSize(), t.getSize());
    
    const ValType* val = t.getVectorValues();
    
    for (unsigned int i=0; i < this->getSize(); i++)
        this->vec_vals[i] += f * val[i];
}



template <typename ValType>
void 
FESystem::Numerics::LocalVector<ValType>::addVal(const std::vector<FESystemUInt>& indices, const FESystem::Numerics::VectorBase<ValType>& t)
{
    FESystemAssert2(t.getSize() == indices.size(), FESystem::Numerics::VectorSizeMismatch, indices.size(), t.getSize());
    
    const ValType* val = t.getVectorValues();
    FESystemUInt n = this->getSize();
    
    for (unsigned int i=0; i < indices.size(); i++)
    {
        FESystemAssert2(indices[i] < n, FESystem::Exception::IndexOutOfBound, indices[i], n);
        this->vec_vals[indices[i]] += val[i];
    }
}


template <typename ValType>
const ValType* 
FESystem::Numerics::LocalVector<ValType>::getVectorValues() const
{
    return &(this->vec_vals[0]);
}



template <typename ValType>
ValType* 
FESystem::Numerics::LocalVector<ValType>::getVectorValues()
{
    return &(this->vec_vals[0]);
}



template <typename ValType>
void
FESystem::Numerics::LocalVector<ValType>::write(std::ostream& out) const
{
    unsigned int m = this->getSize();
    
    out << "Local Vector" << std::endl;
    out << "Size: " << m << std::endl;
    for (unsigned int i=0; i<m; i++)
        out << "#" << i << "  :   " << this->getVal(i) << std::endl; 
}



template <>
void
FESystem::Numerics::LocalVector<FESystemDouble>::complexConjugate()
{
    // no such operation exists for complex
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}


template <>
void
FESystem::Numerics::LocalVector<FESystemFloat>::complexConjugate()
{
    // no such operation exists for complex
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}


template <>
void
FESystem::Numerics::LocalVector<FESystemComplexDouble>::complexConjugate()
{
    for (unsigned int i=0; i< this->getSize(); i++)
        this->vec_vals[i] = std::conj(this->vec_vals[i]);
}



template <>
void
FESystem::Numerics::LocalVector<FESystemComplexFloat>::complexConjugate()
{
    for (unsigned int i=0; i< this->getSize(); i++)
        this->vec_vals[i] = std::conj(this->vec_vals[i]);
}



/***************************************************************************************/
// Template instantiations for some generic classes

INSTANTIATE_CLASS_FOR_ALL_DATA_TYPES(FESystem::Numerics::LocalVector);


/***************************************************************************************/

