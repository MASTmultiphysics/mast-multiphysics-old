//
//  SparseMatrix.cpp
//  FESystem
//
//  Created by Manav Bhatia on 6/13/12.
//  Copyright (c) 2012 Virginia Tech. All rights reserved.
//

// C++ includes
#include <iostream>
#include <iomanip>
#include <algorithm>

// FESystem includes
#include "Numerics/SparseMatrix.h"
#include "Numerics/DenseMatrix.h"
#include "Numerics/SparsityPattern.h"
#include "Numerics/VectorBase.h"
#include "Base/macros.h"
#include "Solvers/Factorizations/HouseHolderTriangulation.h"
#include "Solvers/Factorizations/TriangularBacksubstitution.h"


template <typename ValType>
FESystem::Numerics::SparseMatrix<ValType>::SparseMatrix():
FESystem::Numerics::LocalMatrixBase<ValType>(),
sparsity_pattern(NULL)
{
    
} 



template <typename ValType>
FESystem::Numerics::SparseMatrix<ValType>::~SparseMatrix()
{ 
    
} 



template <typename ValType>
FESystem::Numerics::MatrixType 
FESystem::Numerics::SparseMatrix<ValType>::getType() const
{
    return FESystem::Numerics::LOCAL_SPARSE_MATRIX;
}



template <typename ValType>
const FESystem::Numerics::SparsityPattern& 
FESystem::Numerics::SparseMatrix<ValType>::getSparsityPattern() const
{
    FESystemAssert0(this->sparsity_pattern != NULL, FESystem::Exception::InvalidState);
    return *(this->sparsity_pattern);
}




template <typename ValType>
void 
FESystem::Numerics::SparseMatrix<ValType>::copyMatrixValues(const FESystem::Numerics::MatrixBase<ValType>& m)
{
    const std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    const std::pair<FESystemUInt, FESystemUInt> s_m = m.getSize();
    
    FESystemAssert4((s.first == s_m.first) & (s.second == s_m.second), FESystem::Numerics::MatrixSizeMismatch, s.first, s.second, s_m.first, s_m.second);

    this->zero();

    switch(m.getType())
    {
        case FESystem::Numerics::LOCAL_SPARSE_MATRIX:
        {
            // iterate over the rows and copy the nonzero values
            const FESystem::Numerics::SparsityPattern& m_sparse = *(dynamic_cast<const FESystem::Numerics::SparseMatrix<ValType>&>(m).sparsity_pattern);
            const ValType* m_vals = dynamic_cast<const FESystem::Numerics::SparseMatrix<ValType>&>(m).getMatrixValues();
            std::vector<std::pair<FESystemUInt, FESystemUInt> >::const_iterator source_it, source_end, sink_it, sink_end;
            for (FESystemUInt i=0; i<s_m.first; i++)
            {
                const std::vector<std::pair<FESystemUInt, FESystemUInt> >& source_cols = m_sparse.getAllNonzeroColumnsInRow(i);
                const std::vector<std::pair<FESystemUInt, FESystemUInt> >& sink_cols = this->sparsity_pattern->getAllNonzeroColumnsInRow(i);
                source_it = source_cols.begin(); source_end = source_cols.end();
                sink_it = sink_cols.begin(); sink_end = sink_cols.end();
                this->sparsity_pattern->find(sink_it, sink_end, source_it->first); // find the column iterator
                for ( ; source_it!=source_end; source_it++)
                {
                    if (sink_it->first != source_it->first) // this might happen if the two sparsity patterns are different. 
                    {
                        this->sparsity_pattern->find(sink_it, sink_end, source_it->first); // find the column iterator
                        FESystemAssert1(sink_it!=sink_end, FESystem::Exception::InvalidID, source_it->first);
                    }
                    this->mat_vals[sink_it->second] = m_vals[source_it->second];
                    sink_it++;
                }
            }
        }
            break;
            
            // only sparse matrix is handled for now
        case FESystem::Numerics::LOCAL_DENSE_MATRIX:
        default:
            FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
    }
}



template <typename ValType>
void 
FESystem::Numerics::SparseMatrix<ValType>::copyMatrix(const FESystem::Numerics::MatrixBase<ValType>& m)
{
    switch(m.getType())
    {
        case FESystem::Numerics::LOCAL_SPARSE_MATRIX:
        {
            // make sure that the sparsity pattern for the two are same
            const FESystem::Numerics::SparseMatrix<ValType>& m_sparse = dynamic_cast<const FESystem::Numerics::SparseMatrix<ValType>&>(m);
            if (&(this->sparsity_pattern) != &(m_sparse.sparsity_pattern))
                this->resize(*m_sparse.sparsity_pattern);
            for (FESystemUInt i=0; i<this->n_vals; i++)
                this->mat_vals[i] = m_sparse.mat_vals[i];
        }
            break;
            
            // only sparse matrix is handled for now
        case FESystem::Numerics::LOCAL_DENSE_MATRIX:
        default:
            FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
    }
}



template <typename ValType>
void 
FESystem::Numerics::SparseMatrix<ValType>::copyMatrixTranspose(const FESystem::Numerics::MatrixBase<ValType>& m)
{
    // this is not yet implemented for sparse matrices
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}



template <typename ValType>
void
FESystem::Numerics::SparseMatrix<ValType>::resize(FESystemUInt i, FESystemUInt j) 
{
    // this is not allowed for sparse matrices. The sparse matrix should only be sized through the sparsity pattern
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}



template <typename ValType>
void
FESystem::Numerics::SparseMatrix<ValType>::resize(const FESystem::Numerics::SparsityPattern& s_pattern) 
{
    this->sparsity_pattern = &s_pattern;
    // matrix is square
    this->dims.first = s_pattern.getNDOFs();
    this->dims.second = s_pattern.getNDOFs();
    
    // if the size matches, simply change the dimensions and zero the matrix
    this->resizeValPtr(s_pattern.getNNonzeroValues());
}



template <typename ValType>
void 
FESystem::Numerics::SparseMatrix<ValType>::setToIdentity()
{ 
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    
    FESystemAssert0(s.first == s.second, FESystem::Numerics::MatrixMustBeSquareForOperation);
    
    FESystemUInt id_val=0; FESystemBoolean if_exists = false;
    
    this->zero();
    for (FESystemUInt i=0; i < s.first; i++)
    {
        this->sparsity_pattern->getIDForRowAndColumn(i,i,if_exists,id_val);
        FESystemAssert0(if_exists, FESystem::Exception::InvalidValue);
        this->mat_vals[id_val] = 1.0;
    }
}



template <typename ValType>
ValType
FESystem::Numerics::SparseMatrix<ValType>::getVal(const FESystemUInt i, const FESystemUInt j) const
{ 
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    
    FESystemAssert4( ((i < s.first) && (j < s.second)),
                    FESystem::Numerics::ElementLocationExceedesMatrixSize,
                    s.first, s.second, i, j);
    
    // location ID is obtained from the sparsity pattern
    FESystemUInt id_val=0; FESystemBoolean if_exists = false;
    this->sparsity_pattern->getIDForRowAndColumn(i,j,if_exists,id_val);

    if (if_exists)
        return this->mat_vals[id_val];
    else
        return 0.0;
}  



template <typename ValType>
void 
FESystem::Numerics::SparseMatrix<ValType>::setVal(const FESystemUInt i, const FESystemUInt j, 
                                                  const ValType& val) 
{
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    
    FESystemAssert4( ((i < s.first) && (j < s.second)),
                    FESystem::Numerics::ElementLocationExceedesMatrixSize,
                    s.first, s.second, i, j);
    
    // location ID is obtained from the sparsity pattern
    FESystemUInt id_val=0; FESystemBoolean if_exists = false;
    this->sparsity_pattern->getIDForRowAndColumn(i,j,if_exists,id_val);
    FESystemAssert0(if_exists, FESystem::Exception::InvalidValue);
    
    this->mat_vals[id_val] = val;
}  



template <typename ValType>
void 
FESystem::Numerics::SparseMatrix<ValType>::addVal(const FESystemUInt i, const FESystemUInt j, 
                                                  const ValType& val) 
{
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    
    FESystemAssert4( ((i < s.first) && (j < s.second)),
                    FESystem::Numerics::ElementLocationExceedesMatrixSize,
                    s.first, s.second, i, j);
    
    FESystemUInt id_val=0; FESystemBoolean if_exists = false;
    this->sparsity_pattern->getIDForRowAndColumn(i,j,if_exists,id_val);
    FESystemAssert0(if_exists, FESystem::Exception::InvalidValue);

    // location ID is obtained from the sparsity pattern
    this->mat_vals[id_val] += val;
}  



template <typename ValType>
void
FESystem::Numerics::SparseMatrix<ValType>::addVal(const std::vector<FESystemUInt>& i, 
                                                  const std::vector<FESystemUInt>& j, 
                                                  const FESystem::Numerics::MatrixBase<ValType>& m)
{
    FESystemAssert0(m.getType()==FESystem::Numerics::LOCAL_DENSE_MATRIX, FESystem::Exception::InvalidValue); // currently, only for dense input matrix
    
    FESystemUInt id_val=0, row_n=0, col_n=0; FESystemBoolean if_exists = false;
    std::vector<FESystemUInt>::const_iterator row_it=i.begin(), row_end=i.end(), col_it=j.begin(), col_end=j.end();

    std::pair<FESystemUInt, FESystemUInt> s_m = m.getSize();
    const ValType* m_vals = dynamic_cast<const FESystem::Numerics::DenseMatrix<ValType>&>(m).getMatrixValues(); 
    
    for (; row_it!=row_end; row_it++) 
    {
        col_it=j.begin();
        col_n=0;
        for (; col_it!=col_end; col_it++) 
        {
            this->sparsity_pattern->getIDForRowAndColumn(*row_it,*col_it,if_exists,id_val);
            FESystemAssert0(if_exists, FESystem::Exception::InvalidValue);
            this->mat_vals[id_val] += m_vals[col_n*s_m.first+row_n];
            col_n++;
        }
        row_n++;
    }
}



template <typename ValType>
void 
FESystem::Numerics::SparseMatrix<ValType>::getDiagonal(FESystem::Numerics::VectorBase<ValType>& vec) const
{
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    FESystemUInt vec_size = vec.getSize();
    
    FESystemAssert0(s.first == s.second, FESystem::Numerics::MatrixMustBeSquareForOperation);
    
    FESystemAssert3((vec_size == s.second), FESystem::Numerics::MatrixVectorSizeMismatch,
                    s.first, s.second, vec_size); 
    
    // first zero the vector
    vec.zero();
    ValType* val = vec.getVectorValues();

    FESystemUInt id_val=0; FESystemBoolean if_exists = false;

    // location ID is obtained from the sparsity pattern
    for (FESystemUInt i=0; i < s.first; i++)           
    {
        this->sparsity_pattern->getIDForRowAndColumn(i,i,if_exists,id_val);
        if (if_exists)
            val[i] = this->mat_vals[id_val];
    }
}  



template <typename ValType>
void 
FESystem::Numerics::SparseMatrix<ValType>::setDiagonal(const FESystem::Numerics::VectorBase<ValType>& vec) 
{
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    FESystemUInt vec_size = vec.getSize();
    
    FESystemAssert0(s.first == s.second, FESystem::Numerics::MatrixMustBeSquareForOperation);
    
    FESystemAssert3((vec_size == s.second), FESystem::Numerics::MatrixVectorSizeMismatch,
                    s.first, s.second, vec_size); 
    
    const ValType* val = vec.getVectorValues();
    
    FESystemUInt id_val=0; FESystemBoolean if_exists = false;

    // location ID is obtained from the sparsity pattern
    for (FESystemUInt i=0; i < s.first; i++)
    {
        this->sparsity_pattern->getIDForRowAndColumn(i,i,if_exists,id_val);
        FESystemAssert0(if_exists, FESystem::Exception::InvalidValue);
        this->mat_vals[id_val] = val[i];
    }
}




template <typename ValType>
void 
FESystem::Numerics::SparseMatrix<ValType>::setRowVals(const FESystemUInt row_num, 
                                                      const FESystemUInt col1,
                                                      const FESystemUInt col2,
                                                      const FESystem::Numerics::VectorBase<ValType>& vec) 
{
    // this cannot be used for the sparse matrix
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}



template <typename ValType>
void 
FESystem::Numerics::SparseMatrix<ValType>::getRowVals(const FESystemUInt row_num, 
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
    
    // first zero the vector
    vec.zero();
    
    FESystemUInt id_val=0; FESystemBoolean if_exists = false;

    // location ID is obtained from the sparsity pattern
    for (FESystemUInt j=col1; j<=col2; j++)
    {
        this->sparsity_pattern->getIDForRowAndColumn(row_num,j,if_exists,id_val);
        if (if_exists)
            vec.setVal(j-col1, this->mat_vals[id_val]);
    }
}


template <typename ValType>
void 
FESystem::Numerics::SparseMatrix<ValType>::getColumnVals(const FESystemUInt col_num, 
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
    
    // first zero the vector
    vec.zero();
    
    FESystemUInt id_val=0; FESystemBoolean if_exists = false;

    // location ID is obtained from the sparsity pattern
    for (FESystemUInt j=row1; j<=row2; j++)
    {
        this->sparsity_pattern->getIDForRowAndColumn(j,col_num,if_exists,id_val);
        if (if_exists)
            vec.setVal(j-row1, this->mat_vals[id_val]);
    }
}



template <typename ValType>
void 
FESystem::Numerics::SparseMatrix<ValType>::setColumnVals(const FESystemUInt col_num, 
                                                         const FESystemUInt row1,
                                                         const FESystemUInt row2,
                                                         const FESystem::Numerics::VectorBase<ValType>& vec) 
{
    // this cannot be used for the sparse matrix
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}



template <typename ValType>
void 
FESystem::Numerics::SparseMatrix<ValType>::getSubMatrixVals(const FESystemUInt row1,
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
    
    // first zero the matrix
    m.zero();
    ValType* m_mat = m.getMatrixValues();

    FESystemUInt id_val=0; FESystemBoolean if_exists = false;

    // copy the values from the submatrix to the present matrix
    for (FESystemUInt i=0; i<(row2-row1+1); i++)
        for (FESystemUInt j=0; j<(col2-col1+1); j++)
        {
            this->sparsity_pattern->getIDForRowAndColumn(row1+i, col1+j,if_exists,id_val);
            if (if_exists)
                m_mat[(m_col1+j)*s_sub.first + (m_row1+i)] = this->mat_vals[id_val];
        }
}





template <typename ValType>
void 
FESystem::Numerics::SparseMatrix<ValType>::zeroSubMatrixVals(const FESystemUInt row1, const FESystemUInt row2,
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
    
    std::vector<std::pair<FESystemUInt, FESystemUInt> >::const_iterator it, end;
    
    // copy the values from the submatrix to the present matrix
    for (FESystemUInt i=row1; i<=row2; i++)
    {
        it = this->sparsity_pattern->getAllNonzeroColumnsInRow(i).begin();
        end = this->sparsity_pattern->getAllNonzeroColumnsInRow(i).end();
        
        this->sparsity_pattern->lowerBound(it, end, col1); 
        
        for ( ; (it!=end) && (it->first <= col2); it++) 
            this->mat_vals[it->second] = 0.0;
    }
}



template <typename ValType>
void 
FESystem::Numerics::SparseMatrix<ValType>::setSubMatrixVals(const FESystemUInt row1,
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
    
    // if the matrix to be copied is dense, then copy the values as is, and the user should ensure that 
    // the current matrix has the same nonzeros as the values being copied. Else, an error is thrown.
    // else, if the matrix is sparse, make use of the sparsity patten to reduce the number of values being copied. 
    switch (m.getType())
    {
        case FESystem::Numerics::LOCAL_DENSE_MATRIX:
        {
            const ValType* m_mat = m.getMatrixValues();
            
            // copy the values from the submatrix to the present matrix
            for (FESystemUInt i=0; i<(row2-row1+1); i++)
                for (FESystemUInt j=0; j<(col2-col1+1); j++)
                    this->setVal(row1+i, col1+j, m_mat[(m_col1+j)*s_sub.first + (m_row1+i)]);
        }
            break;
            
        case FESystem::Numerics::LOCAL_SPARSE_MATRIX:
        {
            this->zeroSubMatrixVals(row1, row2, col1, col2);
            
            const FESystem::Numerics::SparseMatrix<ValType>& m_s = dynamic_cast<const FESystem::Numerics::SparseMatrix<ValType>&>(m);
            
            std::vector<std::pair<FESystemUInt, FESystemUInt> >::const_iterator it, end, it_m, end_m;
            
            // copy the values from the submatrix to the present matrix
            for (FESystemUInt i=0; i<=(row2-row1); i++)
            {
                it = this->sparsity_pattern->getAllNonzeroColumnsInRow(row1+i).begin();
                end = this->sparsity_pattern->getAllNonzeroColumnsInRow(row1+i).end();
                it_m = m_s.sparsity_pattern->getAllNonzeroColumnsInRow(m_row1+i).begin();
                end_m = m_s.sparsity_pattern->getAllNonzeroColumnsInRow(m_row1+i).end();

                m_s.sparsity_pattern->lowerBound(it_m, end_m, m_col1);
                this->sparsity_pattern->lowerBound(it, end, col1+(it_m->first-m_col1)); // the first element to be written is same as the lower bound it_m
                                
                for ( ; (it_m!=end_m) && (it_m->first <= m_col2); it_m++) // this is done till one of the two conditions becomes true: col num > col2, or the row does not have any more nonzeros 
                {
                    if ((it->first-col1) != (it_m->first-m_col1))
                        this->sparsity_pattern->find(it, end, col1+(it_m->first-m_col1)); // the column to be written to must exist
                    FESystemAssert1(it!=end, FESystem::Exception::InvalidID, col1+(it_m->first-m_col1));
                    this->mat_vals[it->second] = m_s.mat_vals[it_m->second];
                    it++;
                }
            }
        }
            break;

        default:
            break;
    }
}




template <typename ValType>
void 
FESystem::Numerics::SparseMatrix<ValType>::addSubMatrixVals(const FESystemUInt row1,
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
    
    // if the matrix to be copied is dense, then copy the values as is, and the user should ensure that 
    // the current matrix has the same nonzeros as the values being copied. Else, an error is thrown.
    // else, if the matrix is sparse, make use of the sparsity patten to reduce the number of values being copied. 
    switch (m.getType())
    {
        case FESystem::Numerics::LOCAL_DENSE_MATRIX:
        {
            const ValType* m_mat = m.getMatrixValues();
            
            // copy the values from the submatrix to the present matrix
            for (FESystemUInt i=0; i<(row2-row1+1); i++)
                for (FESystemUInt j=0; j<(col2-col1+1); j++)
                    this->addVal(row1+i, col1+j, m_mat[(m_col1+j)*s_sub.first + (m_row1+i)]);
        }
            break;
            
        case FESystem::Numerics::LOCAL_SPARSE_MATRIX:
        {            
            const FESystem::Numerics::SparseMatrix<ValType>& m_s = dynamic_cast<const FESystem::Numerics::SparseMatrix<ValType>&>(m);
            
            std::vector<std::pair<FESystemUInt, FESystemUInt> >::const_iterator it, end, it_m, end_m;
            
            // copy the values from the submatrix to the present matrix
            for (FESystemUInt i=0; i<=(row2-row1); i++)
            {
                it = this->sparsity_pattern->getAllNonzeroColumnsInRow(row1+i).begin();
                end = this->sparsity_pattern->getAllNonzeroColumnsInRow(row1+i).end();
                it_m = m_s.sparsity_pattern->getAllNonzeroColumnsInRow(m_row1+i).begin();
                end_m = m_s.sparsity_pattern->getAllNonzeroColumnsInRow(m_row1+i).end();
                
                m_s.sparsity_pattern->lowerBound(it_m, end_m, m_col1);
                this->sparsity_pattern->lowerBound(it, end, col1+(it_m->first-m_col1)); // the first element to be written is same as the lower bound it_m
                
                for ( ; (it_m!=end_m) && (it_m->first <= m_col2); it_m++) // this is done till one of the two conditions becomes true: col num > col2, or the row does not have any more nonzeros 
                {
                    if ((it->first-col1) != (it_m->first-m_col1))
                        this->sparsity_pattern->find(it, end, col1+(it_m->first-m_col1)); // the column to be written to must exist
                    FESystemAssert1(it!=end, FESystem::Exception::InvalidID, col1+(it_m->first-m_col1));
                    this->mat_vals[it->second] += f*m_s.mat_vals[it_m->second];
                    it++;
                }
            }
        }
            break;
            
        default:
            break;
    }
}




template <typename ValType>
void
FESystem::Numerics::SparseMatrix<ValType>::getSubMatrixValsFromRowAndColumnIndices(const std::vector<FESystemUInt>& rows,
                                                                                   const std::vector<FESystemUInt>& cols,
                                                                                   const std::map<FESystemUInt, FESystemUInt>& old_to_new_id_map,
                                                                                   FESystem::Numerics::MatrixBase<ValType>& mat) const
{
    // this assumes that the rows and cols are give in ascending order
    
    // works only for sparse matrix
    FESystemAssert0(mat.getType() == FESystem::Numerics::LOCAL_SPARSE_MATRIX, FESystem::Exception::InvalidValue);
    
    FESystem::Numerics::SparseMatrix<ValType>& mat_sparse = dynamic_cast<FESystem::Numerics::SparseMatrix<ValType>&>(mat);
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
    std::vector<FESystemUInt>::const_iterator row_it=rows.begin(), row_end=rows.end(), col_it, col_it2, col_end;
    std::map<FESystemUInt, FESystemUInt>::const_iterator map_row_it, map_col_it, map_end = old_to_new_id_map.end();
    std::vector<std::pair<FESystemUInt, FESystemUInt> >::const_iterator source_c_it, source_c_end, sink_c_it, sink_c_end;;
        
    // copy the values from the submatrix to the present matrix
    for ( ; row_it!=row_end; row_it++)
    {
        map_row_it = old_to_new_id_map.find(*row_it);
        FESystemAssert1(map_row_it != map_end, FESystem::Exception::InvalidID, *row_it);
        const std::vector<std::pair<FESystemUInt, FESystemUInt> >& source_nonzero_cols = this->sparsity_pattern->getAllNonzeroColumnsInRow(*row_it);
        const std::vector<std::pair<FESystemUInt, FESystemUInt> >& sink_nonzero_cols = mat_sparse.sparsity_pattern->getAllNonzeroColumnsInRow(map_row_it->second);
        source_c_it = source_nonzero_cols.begin();
        source_c_end = source_nonzero_cols.end();
        sink_c_it = sink_nonzero_cols.begin();
        sink_c_end = sink_nonzero_cols.end();

        col_it = cols.begin(); 
        col_end = cols.end();
        for (; source_c_it!=source_c_end; source_c_it++)
        {
            // if the column in source sparsity is in the requested set of dofs, then write it
            col_it2 = std::find(col_it, col_end, source_c_it->first);

            if (col_it2 != col_end) // dof was found in the requested set
            {
                // get the id of this dof in the reduced set
                map_col_it = old_to_new_id_map.find(*col_it2);
                FESystemAssert1(map_col_it != map_end, FESystem::Exception::InvalidID, *col_it); // make sure that the id exists in the map
                
                // get the iteratore to the location of this element in the reduced matrix. This is to get the ID in the 1-D sparse value array
                while (sink_c_it != sink_c_end)
                    if (sink_c_it->first == map_col_it->second)
                        break;
                    else 
                        sink_c_it++;
                FESystemAssert1(sink_c_it!=sink_c_end, FESystem::Exception::InvalidID, map_col_it->first); // make sure that col exists in the sub matrix row
                mat_sparse.mat_vals[sink_c_it->second] = this->mat_vals[source_c_it->second];
                col_it = col_it2; // use this as the beginning iterator for the next search
            }
        }
    }
}




template <typename ValType>
void
FESystem::Numerics::SparseMatrix<ValType>::scaleRow(const FESystemInt row_num, const ValType& t)
{
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    
    FESystemAssert2(row_num < s.first,
                    FESystem::Numerics::RowNumberSizeMismatch,
                    s.first, row_num);
    
    const std::vector<std::pair<FESystemUInt, FESystemUInt> >& cols = this->sparsity_pattern->getAllNonzeroColumnsInRow(row_num);
    std::vector<std::pair<FESystemUInt, FESystemUInt> >::const_iterator it = cols.begin(), end=cols.end();
    for ( ; it != end; it++)
        this->mat_vals[it->second] *= t;
}



template <typename ValType>
void
FESystem::Numerics::SparseMatrix<ValType>::scaleColumn(const FESystemInt col_num, const ValType& t)
{
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    
    FESystemAssert2(col_num < s.second,
                    FESystem::Numerics::ColumnNumberSizeMismatch,
                    s.second, col_num);
    
    FESystemUInt id;
    FESystemBoolean if_exists;
    
    for (FESystemUInt i=0; i<s.first; i++)
    {
        this->sparsity_pattern->getIDForRowAndColumn(i, col_num, if_exists, id);
        if (if_exists)
            this->mat_vals[id] *= t;
    }
}




template <typename ValType>
void
FESystem::Numerics::SparseMatrix<ValType>::addScaledSubRow(FESystemUInt row1, ValType f1, FESystemUInt row2, ValType f2, FESystemUInt col1, FESystemUInt col2)
{
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    
    FESystemAssert2(row1 < s.first, FESystem::Exception::DimensionsDoNotMatch, row1, s.first); 
    FESystemAssert2(row2 < s.first, FESystem::Exception::DimensionsDoNotMatch, row2, s.first); 

    FESystemAssert2(col1 < s.second, FESystem::Exception::DimensionsDoNotMatch, col1, s.second); 
    FESystemAssert2(col2 < s.second, FESystem::Exception::DimensionsDoNotMatch, col2, s.second); 
    FESystemAssert2(col1 <= col2, FESystem::Exception::PositiveDifferenceNeeded, col1, col2); 
        
    // get the columns for the two specified rows
    const std::vector<std::pair<FESystemUInt, FESystemUInt> > &cols_row1=this->sparsity_pattern->getAllNonzeroColumnsInRow(row1),
    &cols_row2=this->sparsity_pattern->getAllNonzeroColumnsInRow(row2);
    
    std::vector<std::pair<FESystemUInt, FESystemUInt> >::const_iterator it1=cols_row1.begin(), end1=cols_row1.end(), it2=cols_row2.begin(), end2=cols_row2.end();
    this->sparsity_pattern->find(it1, end1, col1);
    this->sparsity_pattern->find(it2, end2, col1);
    
    for ( ; it2!=end2; it2++)
    {
        if (it2->first <= col2)
        {
            // make sure that the column numbers are consistent
            // if a column in row2 is nonzero, it will not influence row1, and can be skipped
            if (it1->first != it2->first)
                this->sparsity_pattern->find(it1, end1, it2->first);
            FESystemAssert0(it1!=end1, FESystem::Exception::InvalidValue); // row1 should contain the same nonzeros as the row2
            this->mat_vals[it1->second] *= f1;
            this->mat_vals[it1->second] += f2*this->mat_vals[it2->second];
        }
        it1++;
    }
}



template <typename ValType>
void
FESystem::Numerics::SparseMatrix<ValType>::shiftDiagonal(const ValType& v)
{
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    
    FESystemAssert0(s.first == s.second, FESystem::Numerics::MatrixMustBeSquareForOperation);

    FESystemUInt id_val=0; FESystemBoolean if_exists = false;

    // location ID is obtained from the sparsity pattern
    for (FESystemUInt i=0; i < s.first; i++)
    {
        this->sparsity_pattern->getIDForRowAndColumn(i,i,if_exists,id_val);
        FESystemAssert0(if_exists, FESystem::Exception::InvalidValue);
        this->mat_vals[id_val] += v;
    }
}



template <typename ValType>
void 
FESystem::Numerics::SparseMatrix<ValType>::createOrthogonalProjectorOntoVector(const FESystem::Numerics::VectorBase<ValType>& v)
{
    // this gives a dense matrix, and cannot be used with the sparse matrix
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}


template <typename ValType>
void 
FESystem::Numerics::SparseMatrix<ValType>::createOrthogonalProjectorAlongVector(const FESystem::Numerics::VectorBase<ValType>& v)
{
    // this gives a dense matrix, and cannot be used with the sparse matrix
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}


template <typename ValType>
void 
FESystem::Numerics::SparseMatrix<ValType>::createOrthogonalReflectorForVector(const FESystem::Numerics::VectorBase<ValType>& v)
{
    // this gives a dense matrix, and cannot be used with the sparse matrix
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}



template <typename ValType>
ValType
FESystem::Numerics::SparseMatrix<ValType>::dotProductWithColumn(const FESystemUInt col_num, FESystem::Numerics::VectorBase<ValType>& vec) const
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
    FESystemUInt id;
    FESystemBoolean if_exists;
    
    for (FESystemUInt i=0; i<s.first; i++)
    {
        this->sparsity_pattern->getIDForRowAndColumn(i, col_num, if_exists, id);
        if (if_exists)
            v += this->mat_vals[id]*vec_vals[i];
    }
    
    return v;
}




template <typename ValType>
void
FESystem::Numerics::SparseMatrix<ValType>::rightVectorMultiply(const FESystem::Numerics::VectorBase<ValType>& t, FESystem::Numerics::VectorBase<ValType>& res) const
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
    std::vector<std::pair<FESystemUInt, FESystemUInt> >::const_iterator it, end;
    
    for (FESystemUInt i=0; i < s.first; i++)
    {
        const std::vector<std::pair<FESystemUInt, FESystemUInt> >& columns = this->sparsity_pattern->getAllNonzeroColumnsInRow(i);
        it = columns.begin(); end=columns.end();
        
        for ( ; it!=end; it++)
            rvec_vals[i] += this->mat_vals[it->second] * vec_vals[it->first];
    }
}




template <typename ValType>
void 
FESystem::Numerics::SparseMatrix<ValType>::rightSubVectorMultiply(FESystemUInt row1, FESystemUInt row2, FESystemUInt col1, FESystemUInt col2, 
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
        const std::vector<std::pair<FESystemUInt, FESystemUInt> >& columns = this->sparsity_pattern->getAllNonzeroColumnsInRow(i);
        it = columns.begin(); end=columns.end();
        this->sparsity_pattern->lowerBound(it, end, col1);

        rvec_vals[res_row1+(i-row1)] = 0.0;
        for ( ; (it!=end) && (it->first <= col2); it++)
            rvec_vals[res_row1+(i-row1)] += this->mat_vals[it->second] * vec_vals[v_row1+(it->first-col1)];
    }
}




template <typename ValType>
void 
FESystem::Numerics::SparseMatrix<ValType>::leftVectorMultiply(const FESystem::Numerics::VectorBase<ValType>& t, FESystem::Numerics::VectorBase<ValType>& res) const
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
    std::vector<std::pair<FESystemUInt, FESystemUInt> >::const_iterator it, end;

    for (FESystemUInt i=0; i < s.second; i++)
    {
        const std::vector<std::pair<FESystemUInt, FESystemUInt> >& cols = this->sparsity_pattern->getAllNonzeroColumnsInRow(i);
        it = cols.begin(); end=cols.end();
        
        for ( ; it!=end; it++)
            rvec_vals[it->first] += this->mat_vals[it->second] * vec_vals[i];
    }
}


template <typename ValType>
ValType 
FESystem::Numerics::SparseMatrix<ValType>::multiplySubVectorWithSubRow(const FESystemUInt row, const FESystemUInt col1, const FESystemUInt col2,
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

    const std::vector<std::pair<FESystemUInt, FESystemUInt> >& cols = this->sparsity_pattern->getAllNonzeroColumnsInRow(row);
    std::vector<std::pair<FESystemUInt, FESystemUInt> >::const_iterator it1 = cols.begin(), it2, end=cols.end();
    
    // find the iterators for the two columns
    while (it1 != end)
        if (it1->first >= col1) break;
        else it1++;
    if (it1 == end) it1--;
    // look for col2 starting from it1
    it2 = it1;
    while (it2 < end)
    {
        if (it2->first >= col2) break;
        else it2++;
    }
    if (it2 == end) it2--;

    ValType res = 0.0;
    for ( ; it1<=it2; it1++)
        res += this->mat_vals[it1->second] * vec_vals[vec_el1 + (it1->first-col1)]; // the vector location needs to be offset
    
    return res;
}



template <typename ValType>
ValType 
FESystem::Numerics::SparseMatrix<ValType>::multiplySubVectorWithSubColumn(const FESystemUInt col, const FESystemUInt row1, const FESystemUInt row2,
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
    
    FESystemUInt id;
    FESystemBoolean if_exists;

    for (FESystemUInt i=row1; i<=row2; i++)
    {
        this->sparsity_pattern->getIDForRowAndColumn(i, col, if_exists, id);
        if (if_exists)
            res += this->mat_vals[id]*vec_vals[vec_el1+(i-row1)];
    }

    return res;
}


template <typename ValType>
void 
FESystem::Numerics::SparseMatrix<ValType>::matrixRightMultiply(ValType f, const MatrixBase<ValType>& t, MatrixBase<ValType>& res) const
{
    // not yet implemented for this class
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}


template <typename ValType>
void 
FESystem::Numerics::SparseMatrix<ValType>::matrixRightMultiplySubMatrix(FESystemUInt row1, FESystemUInt row2, FESystemUInt col1, FESystemUInt col2,
                                                                        ValType f, const MatrixBase<ValType>& t, MatrixBase<ValType>& res) const
{
    // not yet implemented for this class
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}





template <typename ValType>
void 
FESystem::Numerics::SparseMatrix<ValType>::matrixTransposeRightMultiply(ValType f, const MatrixBase<ValType>& t,
 MatrixBase<ValType>& res) const
{
    // not yet implemented for this class
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}




template <typename ValType>
void 
FESystem::Numerics::SparseMatrix<ValType>::matrixRightMultiplyTranspose(ValType f, const MatrixBase<ValType>& t,
 MatrixBase<ValType>& res) const
{
    // not yet implemented for this class
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}



template <typename ValType>
void 
FESystem::Numerics::SparseMatrix<ValType>::add(ValType f, const MatrixBase<ValType>& t)
{
    // make sure that the two have the same sparsity pattern
    FESystemAssert0(t.getType() == FESystem::Numerics::LOCAL_SPARSE_MATRIX, FESystem::Exception::InvalidValue);
    FESystemAssert0(this->sparsity_pattern == dynamic_cast<const FESystem::Numerics::SparseMatrix<ValType>&>(t).sparsity_pattern, FESystem::Exception::InvalidValue);
    
    // get values of the input matrix
    const ValType* v_in = t.getMatrixValues();
    
    // stored in a column major format: column1 vals -> column2 vals -> ... -> columnN vals
    for (FESystemUInt i=0; i<this->n_vals; i++)
        this->mat_vals[i] += f * v_in[i]; 
    
}




template <typename ValType>
void 
FESystem::Numerics::SparseMatrix<ValType>::initializeLUFactoredMatrices(FESystem::Numerics::MatrixBase<ValType>& l_mat, FESystem::Numerics::MatrixBase<ValType>& u_mat) const
{
    const std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    const std::pair<FESystemUInt, FESystemUInt> s_l = l_mat.getSize();
    const std::pair<FESystemUInt, FESystemUInt> s_u = u_mat.getSize();

    FESystemAssert0(l_mat.getType() == FESystem::Numerics::LOCAL_SPARSE_MATRIX, FESystem::Exception::InvalidValue);
    FESystemAssert0(u_mat.getType() == FESystem::Numerics::LOCAL_SPARSE_MATRIX, FESystem::Exception::InvalidValue);
    FESystemAssert4((s_l.first == s.first) & (s_l.second == s.second), FESystem::Numerics::MatrixSizeMismatch, s.first, s.second, s_l.first, s_l.second);
    FESystemAssert4((s_u.first == s.first) & (s_u.second == s.second), FESystem::Numerics::MatrixSizeMismatch, s.first, s.second, s_u.first, s_u.second);
    

    std::vector<std::pair<FESystemUInt, FESystemUInt> >::const_iterator it_col, end_col;
    
    ValType v=0.0;
    l_mat.zero();
    u_mat.zero();
    
    const FESystem::Numerics::SparsityPattern& l_sparsity = *(dynamic_cast<FESystem::Numerics::SparseMatrix<ValType>&>(l_mat).sparsity_pattern), 
    &u_sparsity = *(dynamic_cast<FESystem::Numerics::SparseMatrix<ValType>&>(u_mat).sparsity_pattern);
    
    FESystem::Numerics::SparsityPattern lu_combined_sparsity;
    lu_combined_sparsity.setNDofs(s.first);
    lu_combined_sparsity.unionWithSparsityPattern(l_sparsity);
    lu_combined_sparsity.unionWithSparsityPattern(u_sparsity);
    lu_combined_sparsity.reinit();

    FESystem::Numerics::SparseMatrix<ValType> lu_combined;
    lu_combined.resize(lu_combined_sparsity);
    
    lu_combined.copyMatrixValues(*this);
    
    // it is assumed that the matrix and its sparsity pattern are initialized
    for (FESystemUInt i=1; i<s.first; i++) // zero everything in the columns before the diagonal 
    {
        const std::vector<std::pair<FESystemUInt, FESystemUInt> >& cols = lu_combined_sparsity.getAllNonzeroColumnsInRow(i);
        it_col = cols.begin();
        end_col = cols.end();
        
        for ( ; it_col->first < i; it_col++) // this will iterate till the sub-diagonal element
        {
            v = lu_combined.mat_vals[it_col->second];
            if (FESystem::Base::magnitude<ValType, typename RealOperationType(ValType)>(v) > 0.0) // perform only if the value left of diagonal is nonzero
            {
                // calculate the factor for this row
                v /= lu_combined.getVal(it_col->first, it_col->first);
                lu_combined.addScaledSubRow(i, 1.0, it_col->first, -v, it_col->first, s.second-1);
                lu_combined.setVal(i, it_col->first, v);
            }
        }
    }
    
    
    // copy the value from the combined lu decomposition matrix
    for (FESystemUInt i=0; i<s.first; i++) // 
    {
        if (i < s.first)
            u_mat.setSubMatrixVals(i, i, i, s.second-1, i, i, i, s.second-1, lu_combined); // upper triangular
        if (i > 0)
            l_mat.setSubMatrixVals(i, i, 0, i-1, i, i, 0, i-1, lu_combined); // lower triangular
    }
    
    l_mat.shiftDiagonal(1.0); // the diagonal is always 1.0 for the L matrix

}


/***************************************************************************************/
// Template instantiations for some generic classes

INSTANTIATE_CLASS_FOR_ALL_DATA_TYPES(FESystem::Numerics::SparseMatrix);


/***************************************************************************************/


