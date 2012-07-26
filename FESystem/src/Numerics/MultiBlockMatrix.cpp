//
//  MultiBlockMatrix.cpp
//  FESystem
//
//  Created by Manav Bhatia on 4/26/12.
//  Copyright (c) 2012. All rights reserved.
//

// FESystem includes
#include "Numerics/MultiBlockMatrix.h"
#include "Base/macros.h"


template <typename ValType> 
FESystem::Numerics::MultiBlockMatrix<ValType>::MultiBlockMatrix():
if_initialized(false),
n_matrices_initialized(0),
row_dimension(0),
column_dimension(0)
{
    
}



template <typename ValType> 
FESystem::Numerics::MultiBlockMatrix<ValType>::~MultiBlockMatrix()
{
    
}



template <typename ValType> 
void
FESystem::Numerics::MultiBlockMatrix<ValType>::clear()
{
    this->if_initialized = false;
    this->n_matrices_initialized = 0;
    this->row_dimension = 0;
    this->column_dimension = 0;
    this->sub_matrices.clear();
    this->rows_per_row_wise_block.clear();
    this->columns_per_column_wise_block.clear();

}




template <typename ValType> 
void
FESystem::Numerics::MultiBlockMatrix<ValType>::setMatrixDimension(FESystemUInt row_size, FESystemUInt column_size)
{
    FESystemAssert0(!this->if_initialized, FESystem::Exception::InvalidState);
    FESystemAssert0(row_size>0, FESystem::Exception::InvalidValue);
    FESystemAssert0(column_size>0, FESystem::Exception::InvalidValue);
    
    this->row_dimension = row_size;
    this->column_dimension = column_size;

    // resize the vectors
    this->sub_matrices.resize(row_size*column_size);
    this->rows_per_row_wise_block.resize(row_size);
    this->columns_per_column_wise_block.resize(column_size);
    
    // now initialize the values
    for (FESystemUInt i=0; i<this->sub_matrices.size(); i++)
        this->sub_matrices[i] = NULL;
    for (FESystemUInt i=0; i<this->rows_per_row_wise_block.size(); i++)
        this->rows_per_row_wise_block[i] = 0;
    for (FESystemUInt i=0; i<this->columns_per_column_wise_block.size(); i++)
        this->columns_per_column_wise_block[i] = 0;
}



template <typename ValType> 
void
FESystem::Numerics::MultiBlockMatrix<ValType>::setSubMatrix(FESystemUInt row_num, FESystemUInt column_num, FESystem::Numerics::MatrixBase<ValType>& matrix)
{
    FESystemUInt mat_num = this->row_dimension*column_num+row_num;
    
    FESystemAssert0(!this->if_initialized, FESystem::Exception::InvalidState);
    FESystemAssert0(this->sub_matrices[mat_num]==NULL, FESystem::Exception::InvalidState);
    
    FESystemAssert2(row_num<this->row_dimension, FESystem::Exception::IndexOutOfBound, row_num, this->row_dimension);
    FESystemAssert2(column_num<this->column_dimension, FESystem::Exception::IndexOutOfBound, column_num, this->column_dimension);

    
    std::pair<FESystemUInt, FESystemUInt> s=matrix.getSize();
    
    // check if the number of rows and colums in the block has been set. If not, set the values, otherwise, make sure that the values match
    // rows
    if (this->rows_per_row_wise_block[row_num] > 0)
    {
        FESystemAssert2(s.first == this->rows_per_row_wise_block[row_num], FESystem::Exception::DimensionsDoNotMatch, 
                        s.first, this->rows_per_row_wise_block[row_num]);
    }
    else
        this->rows_per_row_wise_block[row_num] = s.first;
    // columns
    if (this->columns_per_column_wise_block[column_num] > 0)
    {
        FESystemAssert2(s.second == this->columns_per_column_wise_block[column_num], FESystem::Exception::DimensionsDoNotMatch, 
                        s.second, this->columns_per_column_wise_block[column_num]);
    }
    else
        this->rows_per_row_wise_block[row_num] = s.first;
    
    this->sub_matrices[mat_num] = &matrix;
    
    // increment the counter
    this->n_matrices_initialized++;
    if (this->n_matrices_initialized == this->row_dimension * this->column_dimension)
        this->if_initialized = true;
}




template <typename ValType> 
void 
FESystem::Numerics::MultiBlockMatrix<ValType>::setSubMatrixToNull(FESystemUInt row_num, FESystemUInt column_num)
{
    FESystemUInt mat_num = this->row_dimension*column_num+row_num;
    
    FESystemAssert0(!this->if_initialized, FESystem::Exception::InvalidState);
    FESystemAssert0(this->sub_matrices[mat_num]==NULL, FESystem::Exception::InvalidState);
    
    FESystemAssert2(row_num<this->row_dimension, FESystem::Exception::IndexOutOfBound, row_num, this->row_dimension);
    FESystemAssert2(column_num<this->column_dimension, FESystem::Exception::IndexOutOfBound, column_num, this->column_dimension);
    
    this->sub_matrices[mat_num] = NULL;

    // increment the counter
    this->n_matrices_initialized++;
    if (this->n_matrices_initialized == this->row_dimension * this->column_dimension)
        this->if_initialized = true;
}




template <typename ValType> 
std::pair<FESystemUInt, FESystemUInt> 
FESystem::Numerics::MultiBlockMatrix<ValType>::getBlockIDForIndices(FESystemUInt i, FESystemUInt j) const
{
    FESystemAssert0(!this->if_initialized, FESystem::Exception::InvalidState);

    std::pair<FESystemUInt, FESystemUInt> s=this->getSize();    
    FESystemAssert2(i<s.first, FESystem::Exception::IndexOutOfBound, i, s.first);
    FESystemAssert2(j<s.second, FESystem::Exception::IndexOutOfBound, j, s.second);
    
    std::pair<FESystemUInt, FESystemUInt> m(0,0);
    // identify the row
    FESystemBoolean found = false;
    FESystemUInt num=0;
    while (!found) 
    {
        num+= this->rows_per_row_wise_block[m.first];
        if (i < num)
            found = true;
        else 
            m.first++;
    }

    // identify the column
    found = false;
    num=0;
    while (!found) 
    {
        num+= this->columns_per_column_wise_block[m.second];
        if (j < num)
            found = true;
        else 
            m.second++;
    }
    
    return m;
}



template <typename ValType> 
FESystemBoolean
FESystem::Numerics::MultiBlockMatrix<ValType>::ifBlockIsNull(FESystemUInt i, FESystemUInt j) const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);

    FESystemAssert2(i<this->row_dimension, FESystem::Exception::IndexOutOfBound, i, this->row_dimension);
    FESystemAssert2(j<this->column_dimension, FESystem::Exception::IndexOutOfBound, j, this->column_dimension);
    
    if (this->sub_matrices[j*this->row_dimension+i] == NULL)
        return true;
    else
        return false;
}



template <typename ValType> 
std::pair<FESystemUInt, FESystemUInt> 
FESystem::Numerics::MultiBlockMatrix<ValType>::getRowColumnOffsetsForSubmatrix(FESystemUInt i, FESystemUInt j) const
{
    FESystemAssert0(!this->if_initialized, FESystem::Exception::InvalidState);
    
    FESystemAssert2(i<this->row_dimension, FESystem::Exception::IndexOutOfBound, i, this->row_dimension);
    FESystemAssert2(j<this->column_dimension, FESystem::Exception::IndexOutOfBound, j, this->column_dimension);
    
    std::pair<FESystemUInt, FESystemUInt> m(0,0);
    for (FESystemUInt row=0; row<i; row++) 
        m.first += this->rows_per_row_wise_block[row];
    for (FESystemUInt col=0; col<j; col++) 
        m.second += this->columns_per_column_wise_block[col];

    return m;
}



template <typename ValType> 
FESystem::Numerics::MatrixBase<ValType>& 
FESystem::Numerics::MultiBlockMatrix<ValType>::getBlockForIndices(FESystemUInt i, FESystemUInt j)
{
    FESystemAssert0(!this->if_initialized, FESystem::Exception::InvalidState);
    
    std::pair<FESystemUInt, FESystemUInt> m = this->getBlockIDForIndices(i, j);
    return this->getSubMatrix(m.first, m.second);
}



template <typename ValType> 
FESystem::Numerics::MatrixBase<ValType>& 
FESystem::Numerics::MultiBlockMatrix<ValType>::getSubMatrix(FESystemUInt row_num, FESystemUInt column_num)
{
    FESystemUInt mat_num = this->row_dimension*column_num+row_num;
    
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    FESystemAssert0(this->sub_matrices[mat_num]!=NULL, FESystem::Exception::InvalidState);
    
    if (this->ifBlockIsNull(row_num, column_num))
        FESystemAssert0(false, FESystem::Exception::InvalidValue);
    
    return *(this->sub_matrices[mat_num]);
}



template <typename ValType> 
const FESystem::Numerics::MatrixBase<ValType>& 
FESystem::Numerics::MultiBlockMatrix<ValType>::getSubMatrix(FESystemUInt row_num, FESystemUInt column_num) const
{
    FESystemUInt mat_num = this->row_dimension*column_num+row_num;
    
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    FESystemAssert0(this->sub_matrices[mat_num]!=NULL, FESystem::Exception::InvalidState);
    
    if (this->ifBlockIsNull(row_num, column_num))
        FESystemAssert0(false, FESystem::Exception::InvalidValue);
    
    return *(this->sub_matrices[mat_num]);
}


template <typename ValType> 
FESystem::Numerics::MatrixType
FESystem::Numerics::MultiBlockMatrix<ValType>::getType() const
{
    return FESystem::Numerics::MULTI_BLOCK_MATRIX; 
}


template <typename ValType> 
const std::pair<FESystemUInt, FESystemUInt>
FESystem::Numerics::MultiBlockMatrix<ValType>::getSize() const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    std::pair<FESystemUInt, FESystemUInt> s(0,0);
    for (FESystemUInt i=0; i<this->row_dimension; i++)
        s.first += this->rows_per_row_wise_block[i];
    for (FESystemUInt i=0; i<this->column_dimension; i++)
        s.second += this->columns_per_column_wise_block[i];
 
    return s;
}


template <typename ValType> 
void
FESystem::Numerics::MultiBlockMatrix<ValType>::initializeDenseMatrix(FESystem::Numerics::MatrixBase<ValType>& mat) const
{
    // make sure this is initialized
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    FESystemAssert0(mat.getType() == FESystem::Numerics::LOCAL_DENSE_MATRIX, FESystem::Exception::InvalidValue);

    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    std::pair<FESystemUInt, FESystemUInt> offsets;
    const FESystem::Numerics::MatrixBase<ValType> *mat_p=NULL;
    
    mat.resize(s.first, s.second);
    // iterate over each block and copy the values
    for (FESystemUInt row=0; row<this->row_dimension; row++)
        for (FESystemUInt col=0; col<this->column_dimension; col++)
            if (!this->ifBlockIsNull(row,col))
            {            
                mat_p = &(this->getSubMatrix(row,col));
                offsets = this->getRowColumnOffsetsForSubmatrix(row, col);
                s = mat_p->getSize();
                mat.setSubMatrixVals(offsets.first, offsets.first+s.first-1,
                                     offsets.second, offsets.second+s.second-1,
                                     0, s.first-1, 0, s.second-1,
                                     *mat_p);
            }
}



/***************************************************************************************/
// Template instantiations for some generic classes

INSTANTIATE_CLASS_FOR_ALL_DATA_TYPES(FESystem::Numerics::MultiBlockMatrix);


/***************************************************************************************/

