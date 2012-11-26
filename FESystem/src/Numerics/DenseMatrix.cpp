/*
 *  DenseMatrix.cpp
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/5/10.
 *  Copyright 2010 . All rights reserved.
 *
 */

// C++ includes
#include <iostream>
#include <iomanip>
#include <float.h>

// FESystem includes
#include "Numerics/DenseMatrix.h"
#include "Numerics/VectorBase.h"
#include "Base/macros.h"
#include "Solvers/Factorizations/HouseHolderTriangulation.h"
#include "Solvers/Factorizations/TriangularBacksubstitution.h"


template <typename ValType>
FESystem::Numerics::DenseMatrix<ValType>::DenseMatrix():
FESystem::Numerics::LocalMatrixBase<ValType>()
{ 
    this->dims.first = 0;
    this->dims.second = 0; 
} 



template <typename ValType>
FESystem::Numerics::DenseMatrix<ValType>::~DenseMatrix()
{ 

} 



template <typename ValType>
FESystem::Numerics::MatrixType 
FESystem::Numerics::DenseMatrix<ValType>::getType() const
{
    return FESystem::Numerics::LOCAL_DENSE_MATRIX;
}


template <typename ValType>
const std::pair<FESystemUInt, FESystemUInt>
FESystem::Numerics::DenseMatrix<ValType>::getSize() const
{
    return this->dims;
} 



template <typename ValType>
void 
FESystem::Numerics::DenseMatrix<ValType>::copyMatrixValues(const FESystem::Numerics::MatrixBase<ValType>& m)
{
    std::pair<FESystemUInt, FESystemUInt> s =this->getSize(), s_m = m.getSize();
    FESystemAssert4(((s.first == s_m.first)&&(s.second==s_m.second)), FESystem::Numerics::MatrixSizeMismatch, s.first, s.second, s_m.first, s_m.second);
    
    
    switch(m.getType())
    {
        case FESystem::Numerics::LOCAL_DENSE_MATRIX:
        {
            const ValType* v = m.getMatrixValues();
            
            for (FESystemUInt i=0; i<this->n_vals; i++)
                this->mat_vals[i] = v[i];
        }
            break;
            
        default:
            for (FESystemUInt i=0; i<s.first; i++)
                for (FESystemUInt j=0; j<s.second; j++)
                    this->mat_vals[j*s.first + i] = m.getVal(i,j);
    }
}


template <typename ValType>
void 
FESystem::Numerics::DenseMatrix<ValType>::copyMatrix(const FESystem::Numerics::MatrixBase<ValType>& m)
{
    std::pair<FESystemUInt, FESystemUInt> s = m.getSize();
    this->resize(s.first, s.second);

    switch(m.getType())
    {
        case FESystem::Numerics::LOCAL_DENSE_MATRIX:
        {
            const ValType* v = m.getMatrixValues();
            
            for (FESystemUInt i=0; i<s.first*s.second; i++)
                this->mat_vals[i] = v[i];
        }
            break;
            
        default:
            for (FESystemUInt i=0; i<s.first; i++)
                for (FESystemUInt j=0; j<s.second; j++)
                    this->mat_vals[j*s.first + i] = m.getVal(i,j);
    }
}



template <typename ValType>
void
FESystem::Numerics::DenseMatrix<ValType>::copyRealMatrix(const FESystem::Numerics::MatrixBase<typename RealOperationType(ValType)>& m)
{
    std::pair<FESystemUInt, FESystemUInt> s = m.getSize();
    this->resize(s.first, s.second);
    
    switch(m.getType())
    {
        case FESystem::Numerics::LOCAL_DENSE_MATRIX:
        {
            const typename RealOperationType(ValType)* v = m.getMatrixValues();
            
            for (FESystemUInt i=0; i<s.first*s.second; i++)
                this->mat_vals[i] = ValType(v[i]);
        }
            break;
            
        default:
            for (FESystemUInt i=0; i<s.first; i++)
                for (FESystemUInt j=0; j<s.second; j++)
                    this->mat_vals[j*s.first + i] = m.getVal(i,j);
    }
}



template <typename ValType>
void 
FESystem::Numerics::DenseMatrix<ValType>::copyMatrixTranspose(const FESystem::Numerics::MatrixBase<ValType>& m)
{
    std::pair<FESystemUInt, FESystemUInt> s = m.getSize();
    this->resize(s.second, s.first);
    
    switch(m.getType())
    {
        case FESystem::Numerics::LOCAL_DENSE_MATRIX:
        {
            const ValType* v = m.getMatrixValues();
            
            for (FESystemUInt i=0; i<s.second; i++) // row counter for transpose matrix
                for (FESystemUInt j=0; j<s.first; j++) // column counter for transpose matrix
                    this->mat_vals[j*s.second + i] = v[i*s.first + j];
        }
            break;
            
        default:
            for (FESystemUInt i=0; i<s.second; i++) // row counter for transpose matrix
                for (FESystemUInt j=0; j<s.first; j++) // column counter for transpose matrix
                    this->mat_vals[j*s.second + i] = m.getVal(j,i);
    }
}




template <typename ValType>
void
FESystem::Numerics::DenseMatrix<ValType>::resize(FESystemUInt i, FESystemUInt j) 
{
    // if the size matches, simply change the dimensions and zero the matrix
    this->resizeValPtr(i*j);
    
    // else create the matrix
    this->dims.first = i;
    this->dims.second = j;
} 



template <typename ValType>
void
FESystem::Numerics::DenseMatrix<ValType>::resize(ValType* v, FESystemUInt i, FESystemUInt j) 
{
    // if the size matches, simply change the dimensions and zero the matrix
    this->resizeValPtr(v, i*j);
    
    // else create the matrix
    this->dims.first = i;
    this->dims.second = j;
} 



template <typename ValType>
void 
FESystem::Numerics::DenseMatrix<ValType>::zero() 
{ 
    for (FESystemUInt i=0; i < this->n_vals; i++)
        this->mat_vals[i] = 0.0;
} 



template <typename ValType>
void 
FESystem::Numerics::DenseMatrix<ValType>::setToIdentity() 
{ 
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();

    FESystemAssert0(s.first == s.second, FESystem::Numerics::MatrixMustBeSquareForOperation);

    this->zero();
    for (FESystemUInt i=0; i < s.first; i++)
        this->mat_vals[s.first*i + i] = 1.0;
} 



template <typename ValType>
ValType
FESystem::Numerics::DenseMatrix<ValType>::getVal(const FESystemUInt i, const FESystemUInt j) const
{ 
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    
    FESystemAssert4( ((i < s.first) && (j < s.second)),
                    FESystem::Numerics::ElementLocationExceedesMatrixSize,
                    s.first, s.second, i, j);
    
    // stored in a column major format: column1 vals -> column2 vals -> ... -> columnN vals
    return this->mat_vals[j*s.first + i];
}  



template <typename ValType>
void 
FESystem::Numerics::DenseMatrix<ValType>::setVal(const FESystemUInt i, const FESystemUInt j, 
                                                 const ValType& val) 
{
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    
    FESystemAssert4( ((i < s.first) && (j < s.second)),
                    FESystem::Numerics::ElementLocationExceedesMatrixSize,
                    s.first, s.second, i, j);
    
    // stored in a column major format: column1 vals -> column2 vals -> ... -> columnN vals
    this->mat_vals[j*s.first + i] = val;    
}  



template <typename ValType>
void 
FESystem::Numerics::DenseMatrix<ValType>::addVal(const FESystemUInt i, const FESystemUInt j, 
                                                 const ValType& val) 
{
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    
    FESystemAssert4( ((i < s.first) && (j < s.second)),
                    FESystem::Numerics::ElementLocationExceedesMatrixSize,
                    s.first, s.second, i, j);
    
    // stored in a column major format: column1 vals -> column2 vals -> ... -> columnN vals
    this->mat_vals[j*s.first + i] += val;    
}  


template <typename ValType>
void
FESystem::Numerics::DenseMatrix<ValType>::addVal(const std::vector<FESystemUInt>& i, 
                                                 const std::vector<FESystemUInt>& j, 
                                                 const FESystem::Numerics::MatrixBase<ValType>& m)
{
    const std::pair<FESystemUInt, FESystemUInt> s=this->getSize(), s_m =m.getSize();
    FESystemAssert0(m.getType()==FESystem::Numerics::LOCAL_DENSE_MATRIX, FESystem::Exception::InvalidValue); // currently, only for dense input matrix
    
    const ValType* m_vals = dynamic_cast<const FESystem::Numerics::DenseMatrix<ValType>&>(m).getMatrixValues(); 

    FESystemUInt row_n=0, col_n=0;
    std::vector<FESystemUInt>::const_iterator row_it=i.begin(), row_end=i.end(), col_it=j.begin(), col_end=j.end();    
    
    for (; row_it!=row_end; row_it++) 
    {
        col_it=j.begin();
        col_n=0;
        for (; col_it!=col_end; col_it++) 
        {
            this->mat_vals[(*col_it)*s.first+(*row_it)] += m_vals[col_n*s_m.first+row_n];
            col_n++;
        }
        row_n++;
    }
}


template <typename ValType>
void 
FESystem::Numerics::DenseMatrix<ValType>::getDiagonal(FESystem::Numerics::VectorBase<ValType>& vec) const
{
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    FESystemUInt vec_size = vec.getSize();
    
    FESystemAssert0(s.first == s.second, FESystem::Numerics::MatrixMustBeSquareForOperation);
    
    FESystemAssert3((vec_size == s.second), FESystem::Numerics::MatrixVectorSizeMismatch,
                    s.first, s.second, vec_size); 
    
    ValType* val = vec.getVectorValues();
    
    // stored in a column major format: column1 vals -> column2 vals -> ... -> columnN vals
    for (FESystemUInt i=0; i<s.first; i++)
        val[i] = this->mat_vals[i*s.first + i];    
}  



template <typename ValType>
void 
FESystem::Numerics::DenseMatrix<ValType>::setDiagonal(const FESystem::Numerics::VectorBase<ValType>& vec) 
{
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    FESystemUInt vec_size = vec.getSize();
    
    FESystemAssert0(s.first == s.second, FESystem::Numerics::MatrixMustBeSquareForOperation);
    
    FESystemAssert3((vec_size == s.second), FESystem::Numerics::MatrixVectorSizeMismatch,
                    s.first, s.second, vec_size); 
    
    const ValType* val = vec.getVectorValues();
    
    // stored in a column major format: column1 vals -> column2 vals -> ... -> columnN vals
    for (FESystemUInt i=0; i<s.first; i++)
        this->mat_vals[i*s.first + i] = val[i];    
}  




template <typename ValType>
void 
FESystem::Numerics::DenseMatrix<ValType>::setRowVals(const FESystemUInt row_num, 
                                                     const FESystemUInt col1,
                                                     const FESystemUInt col2,
                                                     const FESystem::Numerics::VectorBase<ValType>& vec) 
{
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    
    FESystemAssert2(row_num < s.first,
                    FESystem::Numerics::RowNumberSizeMismatch,
                    s.first, row_num);
    FESystemAssert2(col1 < s.second,
                    FESystem::Numerics::ColumnNumberSizeMismatch,
                    s.second, col1);
    FESystemAssert2(col2 < s.second,
                    FESystem::Numerics::ColumnNumberSizeMismatch,
                    s.second, col2);
    FESystemAssert2(vec.getSize() == col2-col1+1,
                    FESystem::Exception::DimensionsDoNotMatch,
                    vec.getSize(), col2-col1+1);
    
    // stored in a column major format: column1 vals -> column2 vals -> ... -> columnN vals
    for (FESystemUInt j=col1; j<=col2; j++)
        this->mat_vals[j*s.first + row_num] = vec.getVal(j-col1); 
}



template <typename ValType>
void 
FESystem::Numerics::DenseMatrix<ValType>::getRowVals(const FESystemUInt row_num, 
                                                     const FESystemUInt col1,
                                                     const FESystemUInt col2,
                                                     FESystem::Numerics::VectorBase<ValType>& vec) const
{
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    
    FESystemAssert2(row_num < s.first,
                    FESystem::Numerics::RowNumberSizeMismatch,
                    s.first, row_num);
    FESystemAssert2(col1 < s.second,
                    FESystem::Numerics::ColumnNumberSizeMismatch,
                    s.second, col1);
    FESystemAssert2(col2 < s.second,
                    FESystem::Numerics::ColumnNumberSizeMismatch,
                    s.second, col2);
    FESystemAssert2(vec.getSize() == col2-col1+1,
                    FESystem::Exception::DimensionsDoNotMatch,
                    vec.getSize(), col2-col1+1);
    
    // stored in a column major format: column1 vals -> column2 vals -> ... -> columnN vals
    for (FESystemUInt j=col1; j<=col2; j++)
        vec.setVal(j-col1, this->mat_vals[j*s.first + row_num]);
}


template <typename ValType>
void 
FESystem::Numerics::DenseMatrix<ValType>::getColumnVals(const FESystemUInt col_num, 
                                                        const FESystemUInt row1,
                                                        const FESystemUInt row2,
                                                        FESystem::Numerics::VectorBase<ValType>& vec) const
{
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    
    FESystemAssert2(col_num < s.second,
                    FESystem::Numerics::ColumnNumberSizeMismatch,
                    s.second, col_num);
    FESystemAssert2(row1 < s.first,
                    FESystem::Numerics::RowNumberSizeMismatch,
                    s.first, row1);
    FESystemAssert2(row2 < s.first,
                    FESystem::Numerics::RowNumberSizeMismatch,
                    s.first, row2);
    FESystemAssert2(vec.getSize() == row2-row1+1,
                    FESystem::Exception::DimensionsDoNotMatch,
                    vec.getSize(), row2-row1+1);
    
    // stored in a column major format: column1 vals -> column2 vals -> ... -> columnN vals
    for (FESystemUInt j=row1; j<=row2; j++)
        vec.setVal(j-row1, this->mat_vals[col_num*s.first + j]);
}



template <typename ValType>
void 
FESystem::Numerics::DenseMatrix<ValType>::setColumnVals(const FESystemUInt col_num, 
                                                        const FESystemUInt row1,
                                                        const FESystemUInt row2,
                                                        const FESystem::Numerics::VectorBase<ValType>& vec) 
{
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    
    FESystemAssert2(col_num < s.second,
                    FESystem::Numerics::ColumnNumberSizeMismatch,
                    s.second, col_num);
    FESystemAssert2(row1 < s.first,
                    FESystem::Numerics::RowNumberSizeMismatch,
                    s.first, row1);
    FESystemAssert2(row2 < s.first,
                    FESystem::Numerics::RowNumberSizeMismatch,
                    s.first, row2);
    FESystemAssert2(vec.getSize() == row2-row1+1,
                    FESystem::Exception::DimensionsDoNotMatch,
                    vec.getSize(), row2-row1+1);
    
    // stored in a column major format: column1 vals -> column2 vals -> ... -> columnN vals
    for (FESystemUInt j=row1; j<=row2; j++)
        this->mat_vals[col_num*s.first + j] = vec.getVal(j-row1);
}



template <typename ValType>
void 
FESystem::Numerics::DenseMatrix<ValType>::getSubMatrixVals(const FESystemUInt row1,
                                                           const FESystemUInt row2,
                                                           const FESystemUInt col1,
                                                           const FESystemUInt col2,
                                                           const FESystemUInt m_row1,
                                                           const FESystemUInt m_row2,
                                                           const FESystemUInt m_col1,
                                                           const FESystemUInt m_col2,
                                                           FESystem::Numerics::MatrixBase<ValType>& m) const
{
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    std::pair<FESystemUInt, FESystemUInt> s_sub = m.getSize();
    
    FESystemAssert2(row1 < s.first,
                    FESystem::Numerics::RowNumberSizeMismatch,
                    s.first, row1);
    FESystemAssert2(row2 < s.first,
                    FESystem::Numerics::RowNumberSizeMismatch,
                    s.first, row2);
    FESystemAssert2(col1 < s.second,
                    FESystem::Numerics::ColumnNumberSizeMismatch,
                    s.second, col1);
    FESystemAssert2(col2 < s.second,
                    FESystem::Numerics::ColumnNumberSizeMismatch,
                    s.second, col2);

    FESystemAssert2(m_row1 < s.first,
                    FESystem::Numerics::RowNumberSizeMismatch,
                    s.first, m_row1);
    FESystemAssert2(m_row2 < s_sub.first,
                    FESystem::Numerics::RowNumberSizeMismatch,
                    s_sub.first, m_row2);
    FESystemAssert2(m_col1 < s_sub.second,
                    FESystem::Numerics::ColumnNumberSizeMismatch,
                    s_sub.second, m_col1);
    FESystemAssert2(m_col2 < s_sub.second,
                    FESystem::Numerics::ColumnNumberSizeMismatch,
                    s_sub.second, m_col2);

    FESystemAssert2((row2-row1+1) == (m_row2-m_row1+1),
                    FESystem::Numerics::RowNumberSizeMismatch,
                    (row2-row1+1), (m_row2-m_row1+1));
    FESystemAssert2((col2-col1+1) == (m_col2-m_col1+1),
                    FESystem::Numerics::ColumnNumberSizeMismatch,
                    (col2-col1+1), (m_col2-m_col1+1));
    
    ValType* m_mat = m.getMatrixValues();
    
    // copy the values from the submatrix to the present matrix
    for (FESystemUInt i=0; i<(row2-row1+1); i++)
        for (FESystemUInt j=0; j<(col2-col1+1); j++)
            m_mat[(m_col1+j)*s_sub.first + (m_row1+i)] = this->mat_vals[(col1+j)*s.first + (row1+i)];
}


template <typename ValType>
void 
FESystem::Numerics::DenseMatrix<ValType>::zeroSubMatrixVals(const FESystemUInt row1, const FESystemUInt row2,
                                                            const FESystemUInt col1, const FESystemUInt col2)
{
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    
    FESystemAssert2(row1 < s.first,
                    FESystem::Numerics::RowNumberSizeMismatch,
                    s.first, row1);
    FESystemAssert2(row2 < s.first,
                    FESystem::Numerics::RowNumberSizeMismatch,
                    s.first, row2);
    FESystemAssert2(col1 < s.second,
                    FESystem::Numerics::ColumnNumberSizeMismatch,
                    s.second, col1);
    FESystemAssert2(col2 < s.second,
                    FESystem::Numerics::ColumnNumberSizeMismatch,
                    s.second, col2);
        
    // copy the values from the submatrix to the present matrix
    for (FESystemUInt i=row1; i<=row2; i++)
        for (FESystemUInt j=col1; j<=col2; j++)
            this->mat_vals[j*s.first + i] = 0.0;
}



template <typename ValType>
void 
FESystem::Numerics::DenseMatrix<ValType>::setSubMatrixVals(const FESystemUInt row1,
                                                           const FESystemUInt row2,
                                                           const FESystemUInt col1,
                                                           const FESystemUInt col2,
                                                           const FESystemUInt m_row1,
                                                           const FESystemUInt m_row2,
                                                           const FESystemUInt m_col1,
                                                           const FESystemUInt m_col2,
                                                           const FESystem::Numerics::MatrixBase<ValType>& m)
{
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    std::pair<FESystemUInt, FESystemUInt> s_sub = m.getSize();
    
    FESystemAssert2(row1 < s.first,
                    FESystem::Numerics::RowNumberSizeMismatch,
                    s.first, row1);
    FESystemAssert2(row2 < s.first,
                    FESystem::Numerics::RowNumberSizeMismatch,
                    s.first, row2);
    FESystemAssert2(col1 < s.second,
                    FESystem::Numerics::ColumnNumberSizeMismatch,
                    s.second, col1);
    FESystemAssert2(col2 < s.second,
                    FESystem::Numerics::ColumnNumberSizeMismatch,
                    s.second, col2);
    
    FESystemAssert2(m_row1 < s.first,
                    FESystem::Numerics::RowNumberSizeMismatch,
                    s.first, m_row1);
    FESystemAssert2(m_row2 < s_sub.first,
                    FESystem::Numerics::RowNumberSizeMismatch,
                    s_sub.first, m_row2);
    FESystemAssert2(m_col1 < s_sub.second,
                    FESystem::Numerics::ColumnNumberSizeMismatch,
                    s_sub.second, m_col1);
    FESystemAssert2(m_col2 < s_sub.second,
                    FESystem::Numerics::ColumnNumberSizeMismatch,
                    s_sub.second, m_col2);
    
    FESystemAssert2((row2-row1+1) == (m_row2-m_row1+1),
                    FESystem::Numerics::RowNumberSizeMismatch,
                    (row2-row1+1), (m_row2-m_row1+1));
    FESystemAssert2((col2-col1+1) == (m_col2-m_col1+1),
                    FESystem::Numerics::ColumnNumberSizeMismatch,
                    (col2-col1+1), (m_col2-m_col1+1));
    
    const ValType* m_mat = m.getMatrixValues();
    
    // copy the values from the submatrix to the present matrix
    for (FESystemUInt i=0; i<(row2-row1+1); i++)
        for (FESystemUInt j=0; j<(col2-col1+1); j++)
            this->mat_vals[(col1+j)*s.first + (row1+i)] = m_mat[(m_col1+j)*s_sub.first + (m_row1+i)];
}



template <typename ValType>
void 
FESystem::Numerics::DenseMatrix<ValType>::addSubMatrixVals(const FESystemUInt row1,
                                                           const FESystemUInt row2,
                                                           const FESystemUInt col1,
                                                           const FESystemUInt col2,
                                                           const FESystemUInt m_row1,
                                                           const FESystemUInt m_row2,
                                                           const FESystemUInt m_col1,
                                                           const FESystemUInt m_col2,
                                                           const ValType f, const FESystem::Numerics::MatrixBase<ValType>& m)
{
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    std::pair<FESystemUInt, FESystemUInt> s_sub = m.getSize();
    
    FESystemAssert2(row1 < s.first,
                    FESystem::Numerics::RowNumberSizeMismatch,
                    s.first, row1);
    FESystemAssert2(row2 < s.first,
                    FESystem::Numerics::RowNumberSizeMismatch,
                    s.first, row2);
    FESystemAssert2(col1 < s.second,
                    FESystem::Numerics::ColumnNumberSizeMismatch,
                    s.second, col1);
    FESystemAssert2(col2 < s.second,
                    FESystem::Numerics::ColumnNumberSizeMismatch,
                    s.second, col2);
    
    FESystemAssert2(m_row1 < s.first,
                    FESystem::Numerics::RowNumberSizeMismatch,
                    s.first, m_row1);
    FESystemAssert2(m_row2 < s_sub.first,
                    FESystem::Numerics::RowNumberSizeMismatch,
                    s_sub.first, m_row2);
    FESystemAssert2(m_col1 < s_sub.second,
                    FESystem::Numerics::ColumnNumberSizeMismatch,
                    s_sub.second, m_col1);
    FESystemAssert2(m_col2 < s_sub.second,
                    FESystem::Numerics::ColumnNumberSizeMismatch,
                    s_sub.second, m_col2);
    
    FESystemAssert2((row2-row1+1) == (m_row2-m_row1+1),
                    FESystem::Numerics::RowNumberSizeMismatch,
                    (row2-row1+1), (m_row2-m_row1+1));
    FESystemAssert2((col2-col1+1) == (m_col2-m_col1+1),
                    FESystem::Numerics::ColumnNumberSizeMismatch,
                    (col2-col1+1), (m_col2-m_col1+1));
    
    const ValType* m_mat = m.getMatrixValues();
    
    // copy the values from the submatrix to the present matrix
    for (FESystemUInt i=0; i<(row2-row1+1); i++)
        for (FESystemUInt j=0; j<(col2-col1+1); j++)
            this->mat_vals[(col1+j)*s.first + (row1+i)] += f*m_mat[(m_col1+j)*s_sub.first + (m_row1+i)];
}



template <typename ValType>
void
FESystem::Numerics::DenseMatrix<ValType>::getSubMatrixValsFromRowAndColumnIndices(const std::vector<FESystemUInt>& rows, const std::vector<FESystemUInt>& cols,
                                                                                  const std::map<FESystemUInt, FESystemUInt>& old_to_new_id_map,
                                                                                  FESystem::Numerics::MatrixBase<ValType>& mat) const
{
    mat.zero();
    
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    std::pair<FESystemUInt, FESystemUInt> s_sub = mat.getSize();
    
    FESystemAssert2(rows.size() == s_sub.first, FESystem::Exception::DimensionsDoNotMatch, s_sub.first, rows.size());
    FESystemAssert2(cols.size() == s_sub.second, FESystem::Exception::DimensionsDoNotMatch, s_sub.second, cols.size());
    FESystemAssert2(rows[0] < s.first, FESystem::Exception::IndexOutOfBound, rows[0], s.first);
    FESystemAssert2(*(rows.rbegin()) < s.first, FESystem::Exception::IndexOutOfBound, *(rows.rbegin()), s.first);
    FESystemAssert2(cols[0] < s.second, FESystem::Exception::IndexOutOfBound, cols[0], s.second);
    FESystemAssert2(*(cols.rbegin()) < s.second, FESystem::Exception::IndexOutOfBound, *(cols.rbegin()), s.second);
    
    
    // iterators for row and colums
    std::vector<FESystemUInt>::const_iterator row_it=rows.begin(), row_end=rows.end(), col_it=cols.begin(), col_end=cols.end();
    std::map<FESystemUInt, FESystemUInt>::const_iterator map_row_it, map_col_it, map_end = old_to_new_id_map.end();
    
    ValType* m_mat = mat.getMatrixValues();
    
    // copy the values from the submatrix to the present matrix
    for ( ; row_it!=row_end; row_it++)
    {
        map_row_it = old_to_new_id_map.find(*row_it);
        FESystemAssert1(map_row_it != map_end, FESystem::Exception::InvalidID, *row_it);
        
        col_it = cols.begin();
        for ( ; col_it != col_end; col_it++)
        {
            map_col_it = old_to_new_id_map.find(*col_it);
            FESystemAssert1(map_col_it != map_end, FESystem::Exception::InvalidID, *col_it);
            m_mat[map_col_it->second*s_sub.first+map_row_it->second] = this->mat_vals[(*col_it)*s.first+(*row_it)];
        }
    }
}




template <typename ValType>
void
FESystem::Numerics::DenseMatrix<ValType>::scaleRow(const FESystemInt row_num, const ValType& t)
{
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    
    FESystemAssert2(row_num < s.first,
                    FESystem::Numerics::RowNumberSizeMismatch,
                    s.first, row_num);
    
    // stored in a column major format: column1 vals -> column2 vals -> ... -> columnN vals
    for (FESystemUInt j=0; j<=s.second-1; j++)
        this->mat_vals[j*s.first + row_num] *= t;
}



template <typename ValType>
void
FESystem::Numerics::DenseMatrix<ValType>::scaleColumn(const FESystemInt col_num, const ValType& t)
{
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    
    FESystemAssert2(col_num < s.second,
                    FESystem::Numerics::ColumnNumberSizeMismatch,
                    s.second, col_num);
    
    // stored in a column major format: column1 vals -> column2 vals -> ... -> columnN vals
    for (FESystemUInt j=0; j<=s.first-1; j++)
        this->mat_vals[col_num*s.first + j] *= t;
}



template <typename ValType>
void
FESystem::Numerics::DenseMatrix<ValType>::addScaledSubRow(FESystemUInt row1, ValType f1, FESystemUInt row2, ValType f2, FESystemUInt col1, FESystemUInt col2)
{
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    
    FESystemAssert2(row1 < s.first, FESystem::Exception::DimensionsDoNotMatch, row1, s.first); 
    FESystemAssert2(row2 < s.first, FESystem::Exception::DimensionsDoNotMatch, row2, s.first); 
    
    FESystemAssert2(col1 < s.second, FESystem::Exception::DimensionsDoNotMatch, col1, s.second); 
    FESystemAssert2(col2 < s.second, FESystem::Exception::DimensionsDoNotMatch, col2, s.second); 
    FESystemAssert2(col1 <= col2, FESystem::Exception::PositiveDifferenceNeeded, col1, col2); 
    
    // get the columns for the two specified rows
    for (FESystemUInt i=col1; i<=col2; i++) 
        this->mat_vals[i*s.first+row1] = f1*this->mat_vals[i*s.first+row1] + f2*this->mat_vals[i*s.first+row2];
}



template <typename ValType>
void
FESystem::Numerics::DenseMatrix<ValType>::shiftDiagonal(const ValType& v)
{
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    
    FESystemAssert0(s.first == s.second, FESystem::Numerics::MatrixMustBeSquareForOperation);

    for (FESystemUInt i=0; i < s.first; i++)
        this->mat_vals[i*s.first + i] += v; 
}



template <typename ValType>
void 
FESystem::Numerics::DenseMatrix<ValType>::createOrthogonalProjectorOntoVector(const FESystem::Numerics::VectorBase<ValType>& v)
{
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    
    FESystemAssert2(v.getSize() == s.first,
                    FESystem::Numerics::RowNumberSizeMismatch,
                    s.first, v.getSize());
    FESystemAssert2(v.getSize() == s.second,
                    FESystem::Numerics::ColumnNumberSizeMismatch,
                    s.second, v.getSize());
    
    const ValType* v_val = v.getVectorValues();
    
    for (FESystemUInt i=0; i < s.first; i++)
        for (FESystemUInt j=0; j < s.second; j++)
            this->mat_vals[j*s.first + i] = v_val[i] * v_val[j];
    
    this->scale(1.0/pow(v.getL2Norm(),2));
}


template <typename ValType>
void 
FESystem::Numerics::DenseMatrix<ValType>::createOrthogonalProjectorAlongVector(const FESystem::Numerics::VectorBase<ValType>& v)
{
    // first create the projector onto the vector
    this->createOrthogonalProjectorOntoVector(v);
    
    // now, scale it by -1, and shift the diagonal by 1.0
    this->scale(-1.0);
    this->shiftDiagonal(1.0);
}


template <typename ValType>
void 
FESystem::Numerics::DenseMatrix<ValType>::createOrthogonalReflectorForVector(const FESystem::Numerics::VectorBase<ValType>& v)
{
    std::auto_ptr<FESystem::Numerics::VectorBase<ValType> > 
    v1(FESystem::Numerics::VectorCreate<ValType>(v.getType()).release()), 
    v2(FESystem::Numerics::VectorCreate<ValType>(v.getType()).release());
    
    v1->resize(v.getSize());
    v2->resize(v.getSize());

    // for numerical stability, choose the reflector with the longer reflection distance
    v1->setVal(0, -1.0 * v.getL2Norm());
    v1->add(1.0, v);

    v2->setVal(0, 1.0 * v.getL2Norm());
    v2->add(1.0, v);

    // create the projector onto the longer vector
    if (v1->getL2Norm() > v2->getL2Norm())
        this->createOrthogonalProjectorOntoVector(*v1);
    else 
        this->createOrthogonalProjectorOntoVector(*v2);
    
    // now, scale it by -2, and shift the diagonal by 1.0
    // the scaling of 2.0 is done because the projector along vector (|v|e_1 - v) would bisect the angle between 
    // (v, 0 and |v|e_1), but we need the reflector to project the vector v onto the e_1 axis. 
    this->scale(-2.0);
    this->shiftDiagonal(1.0);
}



template <typename ValType>
ValType
FESystem::Numerics::DenseMatrix<ValType>::dotProductWithColumn(const FESystemUInt col_num, FESystem::Numerics::VectorBase<ValType>& vec) const
{
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    
    FESystemAssert2(col_num < s.second,
                    FESystem::Numerics::ColumnNumberSizeMismatch,
                    s.second, col_num);
    FESystemAssert3(vec.getSize() == s.first,
                    FESystem::Numerics::MatrixVectorSizeMismatch,
                    s.first, s.second, vec.getSize());

    const ValType* vec_vals = vec.getVectorValues();

    ValType v=0;
    for (FESystemUInt i=0; i<s.first; i++)
        v += vec_vals[i] * this->mat_vals[s.first*col_num +i];
    
    return v;
}



template <typename ValType>
void 
FESystem::Numerics::DenseMatrix<ValType>::rightVectorMultiply(const FESystem::Numerics::VectorBase<ValType>& t,
                                                              FESystem::Numerics::VectorBase<ValType>& res) const
{
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    FESystemUInt vec_size = t.getSize();
    FESystemUInt rvec_size = res.getSize();
    
    FESystemAssert3((vec_size == s.second), FESystem::Numerics::MatrixVectorSizeMismatch,
                    s.first, s.second, vec_size); 
    FESystemAssert3((rvec_size == s.first), FESystem::Numerics::MatrixVectorSizeMismatch,
                    s.first, s.second, rvec_size); 
    
    res.zero();
    ValType* rvec_vals = res.getVectorValues();

    // get the vector components
    const ValType* vec_vals = t.getVectorValues();
    
    for (FESystemUInt i=0; i < s.first; i++) 
        for (FESystemUInt j=0; j < s.second; j++)     
            rvec_vals[i] += this->mat_vals[j*s.first + i] * vec_vals[j];
}




template <typename ValType>
void 
FESystem::Numerics::DenseMatrix<ValType>::rightSubVectorMultiply(FESystemUInt row1, FESystemUInt row2, FESystemUInt col1, FESystemUInt col2, 
                                                                  FESystemUInt v_row1, FESystemUInt v_row2, FESystemUInt res_row1, FESystemUInt res_row2,
                                                                  const FESystem::Numerics::VectorBase<ValType>& t, FESystem::Numerics::VectorBase<ValType>& res) const
{
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    FESystemUInt vec_size = t.getSize();
    FESystemUInt rvec_size = res.getSize();
    
    FESystemAssert2(row2>=row1, FESystem::Exception::PositiveDifferenceNeeded, row2, row1);
    FESystemAssert2(col2>=col1, FESystem::Exception::PositiveDifferenceNeeded, col2, col1);
    FESystemAssert2((row2-row1) == (res_row2-res_row1), FESystem::Exception::DimensionsDoNotMatch, (row2-row1), (res_row2-res_row1));
    FESystemAssert2((col2-col1) == (v_row2-v_row1), FESystem::Exception::DimensionsDoNotMatch, (col2-col1), (v_row2-v_row1));
    FESystemAssert2((v_row2-v_row1+1) <= vec_size, FESystem::Exception::DimensionsDoNotMatch, (v_row2-v_row1+1), vec_size);
    FESystemAssert2((res_row2-res_row1+1) <= rvec_size, FESystem::Exception::DimensionsDoNotMatch, (res_row2-res_row1+1), rvec_size);
    
    ValType* rvec_vals = res.getVectorValues();
    
    // get the vector components
    const ValType* vec_vals = t.getVectorValues();
    std::vector<std::pair<FESystemUInt, FESystemUInt> >::const_iterator it, end;
    
    for (FESystemUInt i=row1; i <= row2; i++)
    {
        rvec_vals[res_row1+(i-row1)] = 0.0;
        for (FESystemUInt j=col1; j<=col2; j++)
            rvec_vals[res_row1+(i-row1)] += this->mat_vals[j*s.first+i] * vec_vals[v_row1+(j-col1)];
    }
}



template <typename ValType>
void 
FESystem::Numerics::DenseMatrix<ValType>::leftVectorMultiply(const FESystem::Numerics::VectorBase<ValType>& t,
                                                             FESystem::Numerics::VectorBase<ValType>& res) const
{
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    FESystemUInt vec_size = t.getSize();
    FESystemUInt rvec_size = res.getSize();
    
    FESystemAssert3((vec_size == s.first), FESystem::Numerics::MatrixVectorSizeMismatch,
                    s.first, s.second, vec_size); 
    FESystemAssert3((rvec_size == s.second), FESystem::Numerics::MatrixVectorSizeMismatch,
                    s.first, s.second, rvec_size); 
    
    // get the vector components 
    const ValType* vec_vals = t.getVectorValues();
    res.zero();
    ValType* rvec_vals = res.getVectorValues();
    
    for (FESystemUInt i=0; i < s.second; i++) 
        for (FESystemUInt j=0; j < s.first; j++)     
            rvec_vals[i] += this->mat_vals[i*s.first + j] * vec_vals[j];
}


template <typename ValType>
ValType 
FESystem::Numerics::DenseMatrix<ValType>::multiplySubVectorWithSubRow(const FESystemUInt row, const FESystemUInt col1, const FESystemUInt col2,
                                                                       const FESystemUInt vec_el1, const FESystemUInt vec_el2, const FESystem::Numerics::VectorBase<ValType>& vec) const
{
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    FESystemUInt vec_size = vec.getSize();
    
    FESystemAssert2(row < s.first, FESystem::Exception::DimensionsDoNotMatch, row, s.first); 
    FESystemAssert2(col1 < s.second, FESystem::Exception::DimensionsDoNotMatch, col1, s.second); 
    FESystemAssert2(col2 < s.second, FESystem::Exception::DimensionsDoNotMatch, col2, s.second); 
    FESystemAssert2(col1 <= col2, FESystem::Exception::PositiveDifferenceNeeded, col1, col2); 
    
    FESystemAssert2(vec_el1 < vec_size, FESystem::Exception::DimensionsDoNotMatch, vec_el1, vec_size); 
    FESystemAssert2(vec_el2 < vec_size, FESystem::Exception::DimensionsDoNotMatch, vec_el2, vec_size); 
    FESystemAssert2(vec_el1 <= vec_el2, FESystem::Exception::PositiveDifferenceNeeded, vec_el1, vec_el2); 
    
    FESystemAssert2((col2-col1) == (vec_el2-vec_el1), FESystem::Exception::DimensionsDoNotMatch, (col2-col1), (vec_el2-vec_el1)); 
    
    // get the vector components 
    const ValType* vec_vals = vec.getVectorValues();
    ValType res = 0.0;
        
    for (FESystemUInt i=0; i<=(col2-col1); i++)
        res += this->mat_vals[(col1+i)*s.first + row] * vec_vals[vec_el1 + i]; // the vector location needs to be offset
    
    return res;
}



template <typename ValType>
ValType 
FESystem::Numerics::DenseMatrix<ValType>::multiplySubVectorWithSubColumn(const FESystemUInt col, const FESystemUInt row1, const FESystemUInt row2,
                                                                          const FESystemUInt vec_el1, const FESystemUInt vec_el2, const FESystem::Numerics::VectorBase<ValType>& vec) const
{
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    FESystemUInt vec_size = vec.getSize();
    
    FESystemAssert2(col < s.second, FESystem::Exception::DimensionsDoNotMatch, col, s.second); 
    FESystemAssert2(row1 < s.first, FESystem::Exception::DimensionsDoNotMatch, row1, s.first); 
    FESystemAssert2(row2 < s.first, FESystem::Exception::DimensionsDoNotMatch, row2, s.first); 
    FESystemAssert2(row1 <= row2, FESystem::Exception::PositiveDifferenceNeeded, row1, row2); 
    
    FESystemAssert2(vec_el1 < vec_size, FESystem::Exception::DimensionsDoNotMatch, vec_el1, vec_size); 
    FESystemAssert2(vec_el2 < vec_size, FESystem::Exception::DimensionsDoNotMatch, vec_el2, vec_size); 
    FESystemAssert2(vec_el1 <= vec_el2, FESystem::Exception::PositiveDifferenceNeeded, vec_el1, vec_el2); 
    
    FESystemAssert2((row2-row1) == (vec_el2-vec_el1), FESystem::Exception::DimensionsDoNotMatch, (row2-row1), (vec_el2-vec_el1)); 
    
    // get the vector components 
    const ValType* vec_vals = vec.getVectorValues();
    ValType res = 0.0;
        
    for (FESystemUInt i=0; i<=(row2-row1); i++)
        res += this->mat_vals[col*s.first + i+row1] * vec_vals[vec_el1 + i]; // the vector location needs to be offset
    
    return res;
}



template <typename ValType>
void 
FESystem::Numerics::DenseMatrix<ValType>::matrixRightMultiplySubMatrix(FESystemUInt row1, FESystemUInt row2, FESystemUInt col1, FESystemUInt col2,
                                                                       ValType f, const MatrixBase<ValType>& t, MatrixBase<ValType>& res) const
{
    std::pair<FESystemUInt, FESystemUInt> m1_s = this->getSize();
    std::pair<FESystemUInt, FESystemUInt> m2_s = t.getSize();
    std::pair<FESystemUInt, FESystemUInt> res_s = res.getSize();
    
    FESystemAssert2(row1 <= row2, FESystem::Exception::PositiveDifferenceNeeded, row1, row2);
    FESystemAssert2(col1 <= col2, FESystem::Exception::PositiveDifferenceNeeded, col1, col2);
    FESystemAssert4(m1_s.second == row2-row1+1, FESystem::Numerics::MatrixMultiplicationSizeMismatch,
                    m1_s.first, m1_s.second, row2-row1+1, col2-col1+1); 
    FESystemAssert4(res_s.first == m1_s.first, FESystem::Numerics::MatrixMultiplicationSizeMismatch,
                    m1_s.first, m1_s.second, res_s.first, res_s.second); 
    FESystemAssert4(res_s.second == col2-col1+1, FESystem::Numerics::MatrixMultiplicationSizeMismatch,
                    m2_s.first, m2_s.second, row2-row1+1, col2-col1+1); 
    
    // get the vector components 
    const ValType* m2_vals = t.getMatrixValues();
    res.zero();
    ValType* r_vals = res.getMatrixValues();
    
    for (FESystemUInt i=0; i< row2-row1+1; i++)
        for (FESystemUInt j=0; j < col2-col1+1; j++) 
        {
            for (FESystemUInt k=0; k < m1_s.second; k++) 
                r_vals[j*res_s.first + i] += this->mat_vals[k*m1_s.first + i] * m2_vals[(j+col1)*m2_s.first + k+row1];
            r_vals[j*res_s.first + i] *= f;
        }
}



template <typename ValType>
void
FESystem::Numerics::DenseMatrix<ValType>::matrixRightMultiply (ValType f, const MatrixBase<ValType>& t,
                                                               MatrixBase<ValType>& res) const
{
    std::pair<FESystemUInt, FESystemUInt> m1_s = this->getSize();
    std::pair<FESystemUInt, FESystemUInt> m2_s = t.getSize();
    std::pair<FESystemUInt, FESystemUInt> res_s = res.getSize();
    
    FESystemAssert4((m1_s.second == m2_s.first),
                    FESystem::Numerics::MatrixMultiplicationSizeMismatch,
                    m1_s.first, m1_s.second, m2_s.first, m2_s.second);
    FESystemAssert4((res_s.first == m1_s.first),
                    FESystem::Numerics::MatrixMultiplicationSizeMismatch,
                    m1_s.first, m1_s.second, res_s.first, res_s.second);
    FESystemAssert4((res_s.second == m2_s.second),
                    FESystem::Numerics::MatrixMultiplicationSizeMismatch,
                    m2_s.first, m2_s.second, res_s.first, res_s.second);
    
    // get the vector components
    const ValType* m2_vals = t.getMatrixValues();
    res.zero();
    ValType* r_vals = res.getMatrixValues();
    
    for (FESystemUInt i=0; i< m1_s.first; i++)
        for (FESystemUInt j=0; j < m2_s.second; j++)
        {
            for (FESystemUInt k=0; k < m1_s.second; k++)
                r_vals[j*res_s.first + i] += this->mat_vals[k*m1_s.first + i] * m2_vals[j*m2_s.first + k];
            r_vals[j*res_s.first + i] *= f;
        }
}




template <typename ValType>
void 
FESystem::Numerics::DenseMatrix<ValType>::matrixTransposeRightMultiply (ValType f, const MatrixBase<ValType>& t,
                                                                        MatrixBase<ValType>& res) const
{
    std::pair<FESystemUInt, FESystemUInt> m1_s = this->getSize();
    std::pair<FESystemUInt, FESystemUInt> m2_s = t.getSize();
    std::pair<FESystemUInt, FESystemUInt> res_s = res.getSize();
    
    FESystemAssert4((m1_s.first == m2_s.first), 
                    FESystem::Numerics::MatrixMultiplicationSizeMismatch,
                    m1_s.first, m1_s.second, m2_s.first, m2_s.second); 
    FESystemAssert4((res_s.first == m1_s.second), 
                    FESystem::Numerics::MatrixMultiplicationSizeMismatch,
                    m1_s.first, m1_s.second, res_s.first, res_s.second); 
    FESystemAssert4((res_s.second == m2_s.second), 
                    FESystem::Numerics::MatrixMultiplicationSizeMismatch,
                    m2_s.first, m2_s.second, res_s.first, res_s.second); 
    
    // get the vector components 
    const ValType* m2_vals = t.getMatrixValues();
    res.zero();
    ValType* r_vals = res.getMatrixValues();
    
    for (FESystemUInt i=0; i< m1_s.second; i++)
        for (FESystemUInt j=0; j < m2_s.second; j++) 
        {
            for (FESystemUInt k=0; k < m1_s.first; k++) 
                r_vals[j*res_s.first + i] += this->mat_vals[i*m1_s.first + k] *  m2_vals[j*m2_s.first + k];
            r_vals[j*res_s.first + i] *= f;
        }
}




template <typename ValType>
void 
FESystem::Numerics::DenseMatrix<ValType>::matrixRightMultiplyTranspose (ValType f, const MatrixBase<ValType>& t,
                                                                        MatrixBase<ValType>& res) const
{
    std::pair<FESystemUInt, FESystemUInt> m1_s = this->getSize();
    std::pair<FESystemUInt, FESystemUInt> m2_s = t.getSize();
    std::pair<FESystemUInt, FESystemUInt> res_s = res.getSize();
    
    FESystemAssert4((m1_s.second == m2_s.second), 
                    FESystem::Numerics::MatrixMultiplicationSizeMismatch,
                    m1_s.first, m1_s.second, m2_s.second, m2_s.first); 
    FESystemAssert4((res_s.first == m1_s.first), 
                    FESystem::Numerics::MatrixMultiplicationSizeMismatch,
                    m1_s.first, m1_s.second, res_s.first, res_s.second); 
    FESystemAssert4((res_s.second == m2_s.first), 
                    FESystem::Numerics::MatrixMultiplicationSizeMismatch,
                    m2_s.first, m2_s.second, res_s.first, res_s.second); 
    
    // get the vector components 
    const ValType* m2_vals = t.getMatrixValues();
    res.zero();
    ValType* r_vals = res.getMatrixValues();
    
    for (FESystemUInt i=0; i< m1_s.first; i++)
        for (FESystemUInt j=0; j < m2_s.first; j++) 
        {
            for (FESystemUInt k=0; k < m1_s.second; k++) 
                r_vals[j*res_s.first + i] += this->mat_vals[k*m1_s.first + i] * m2_vals[k*m2_s.first + j];
            r_vals[j*res_s.first + i] *= f;
        }
}



template <typename ValType>
void 
FESystem::Numerics::DenseMatrix<ValType>::add(ValType f, const MatrixBase<ValType>& t)
{
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize(),
    s_in = t.getSize();
    
    FESystemAssert4(s.first == s_in.first,
                    FESystem::Numerics::MatrixSizeMismatch,
                    s.first, s.second,
                    s_in.first, s_in.second);
    
    // get values of the input matrix
    const ValType* v_in = t.getMatrixValues();
    
    // stored in a column major format: column1 vals -> column2 vals -> ... -> columnN vals
    for (FESystemUInt i=0; i<s.first*s.second; i++)
        this->mat_vals[i] += f * v_in[i]; 

}



template <typename ValType>
void 
FESystem::Numerics::DenseMatrix<ValType>::initializeLUFactoredMatrices(FESystem::Numerics::MatrixBase<ValType>& l_mat, FESystem::Numerics::MatrixBase<ValType>& u_mat) const
{
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
    
//    const std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
//    const std::pair<FESystemUInt, FESystemUInt> s_l = l_mat.getSize();
//    const std::pair<FESystemUInt, FESystemUInt> s_u = u_mat.getSize();
//    
//    FESystemAssert0(l_mat.getType() == FESystem::Numerics::LOCAL_SPARSE_MATRIX, FESystem::Exception::InvalidValue);
//    FESystemAssert0(u_mat.getType() == FESystem::Numerics::LOCAL_SPARSE_MATRIX, FESystem::Exception::InvalidValue);
//    FESystemAssert4((s_l.first == s.first) & (s_l.second == s.second), FESystem::Numerics::MatrixSizeMismatch, s.first, s.second, s_l.first, s_l.second);
//    FESystemAssert4((s_u.first == s.first) & (s_u.second == s.second), FESystem::Numerics::MatrixSizeMismatch, s.first, s.second, s_u.first, s_u.second);
//    
//    
//    std::set<FESystemUInt>::const_iterator it_row, end_row, it_col, end_col;
//    ValType v=0.0, v2=0.0;
//    l_mat.zero();
//    u_mat.zero();
//    
//    // copy the upper triangular matrix to u_mat
//    for (FESystemUInt i=1; i<s.first; i++) // 
//    {
//        std::set<FESystemUInt> cols = this->sparsity_pattern->getAllNonzeroColumnsInRow(i);
//        it_col = cols.begin(); end_col = cols.end();
//        for ( ; it_col != end_col; it_col++)
//            if (*it_col >= i)
//                u_mat.setVal(i, *it_col, this->getVal(i, *it_col));
//    }
//    
//    // it is assumed that the matrix and its sparsity pattern are initialized
//    for (FESystemUInt i=1; i<s.first; i++) // zero everything in the i^th column for all rows > i
//    {
//        std::set<FESystemUInt> cols = this->sparsity_pattern->getAllNonzeroColumnsInRow(i);
//        it_col = cols.begin(); end_col = cols.end();
//        
//        std::set<FESystemUInt> rows = this->sparsity_pattern->getAllNonzeroRowsInColumn(i);
//        it_row = rows.begin(); end_row = rows.end();
//        
//        for ( ; it_row != end_row; it_row++) // these are all the nonzero rows in this column
//            if (*it_row > i) // factors are based on the lower triangular portion of the matrix
//            {
//                v = -u_mat.getVal(*it_row, i);
//                if (FESystem::Base::magnitude<ValType, typename RealOperationType(ValType)>(v) > 0.0) // perform only if the value left of diagonal is nonzero
//                {
//                    // calculate the factor for this row
//                    v /= this->getVal(i, i);
//                    // set the factor in l_mat
//                    l_mat.setVal(*it_row, i, v);
//                    
//                    // now iterate over all the columns in this row and set the values
//                    it_col = cols.begin();
//                    for ( ; it_col != end_col; it_col++)
//                        if (*it_col >= i) // only for upper triangular entries
//                        {
//                            v2 = u_mat.getVal(*it_row, *it_col) + v*u_mat.getVal(i,*it_col); 
//                            u_mat.setVal(*it_row, *it_col, v2);
//                        }
//                }
//            }
//    }
//    
//    l_mat.shiftDiagonal(1.0); // the diagonal is always 1.0 for the L matrix
}


/***************************************************************************************/
// Template instantiations for some generic classes

INSTANTIATE_CLASS_FOR_ALL_DATA_TYPES(FESystem::Numerics::DenseMatrix);


/***************************************************************************************/
