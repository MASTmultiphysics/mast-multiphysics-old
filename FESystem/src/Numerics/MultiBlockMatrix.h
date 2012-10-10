//
//  MultiBlockMatrix.h
//  FESystem
//
//  Created by Manav Bhatia on 4/26/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem_multi_block_matrix_h__
#define __fesystem_multi_block_matrix_h__

// C++ includes
#include <vector>

// FESystem includes
#include "Numerics/MatrixBase.h"


namespace FESystem
{
    namespace Numerics
    {
        
        /*!
         *   This class allows construction of a matrix using multiple sub-matrix 
         *   blocks. The matrix blocks should be initialized and given to this matrix.
         */
        
        template <typename ValType> 
        class MultiBlockMatrix 
        {
        public:
            
            /*!
             *   constructor
             */
            MultiBlockMatrix();
            
            virtual ~MultiBlockMatrix();
            
            /*!
             *  this method clears the data structures
             */
            virtual void clear();
            
            /*!
             *   sets the number of sub-matrices in the row and column directions respectively
             */
            void setMatrixDimension(FESystemUInt row_size, FESystemUInt column_size);
            
            
            /*!
             *   sets the matrix for the given row and colum index of the sub-matrix. For sub-matrices
             *   that will be zero, they should be explicitly set to NULL using the method setSubMatrixToNull.
             *   The row and column numbers start at zero. 
             */
            void setSubMatrix(FESystemUInt row_num, FESystemUInt column_num, FESystem::Numerics::MatrixBase<ValType>& matrix);
            
            /*!
             *  sets the sub matrix to a null matrix. It is important to call this method for all 
             *  sub-matrices which will be left as zero. The row and column numbers start at zero. 
             */
            void setSubMatrixToNull(FESystemUInt row_num, FESystemUInt column_num);
            
            /*!
             *  this returns a reference to the sub-matrix object that contains the 
             *  multi-block-sparse-matrix object
             */
            FESystem::Numerics::MatrixBase<ValType>& getSubMatrix(FESystemUInt row_num, FESystemUInt column_num);

            /*!
             *  this returns a constant reference to the sub-matrix object that contains the 
             *  multi-block-sparse-matrix object
             */
            const FESystem::Numerics::MatrixBase<ValType>& getSubMatrix(FESystemUInt row_num, FESystemUInt column_num) const;

            /*!
             *   Returns the type of matrix
             */
            virtual FESystem::Numerics::MatrixType getType() const;
            
            /*!
             *   Returns the number of rows and columns of this matrix
             */
            virtual const std::pair<FESystemUInt, FESystemUInt> getSize() const;
            
            /*!
             *   Initializes \p mat to a dense matrix containing the updated information of the block matrices.
             */
            void initializeDenseMatrix(FESystem::Numerics::MatrixBase<ValType>& mat) const;
            
    protected:
            
            /*!
             *   Returns true if the block is set to NULL, otherwise, false
             */
            FESystemBoolean ifBlockIsNull(FESystemUInt i, FESystemUInt j) const;
            
            /*!
             *   Returns the block row and column numbers for the given matrix row and column numbers
             */
            std::pair<FESystemUInt, FESystemUInt> getBlockIDForIndices(FESystemUInt i, FESystemUInt j) const;
            
            /*!
             *   Returns the first row and column number in the multiblock matrix for the block defined by i and j 
             */
            std::pair<FESystemUInt, FESystemUInt> getRowColumnOffsetsForSubmatrix(FESystemUInt i, FESystemUInt j) const;

            
            /*!
             *   Returns the sub-matrix block for the given matrix row and column numbers
             */
            FESystem::Numerics::MatrixBase<ValType>& getBlockForIndices(FESystemUInt i, FESystemUInt j);

            /*!
             *   stores the state of initialization
             */
            FESystemBoolean if_initialized;
            
            /*!
             *  number of sub-matrices that have been initialized. Used to track the status of initialization.
             */
            FESystemUInt n_matrices_initialized;
            
            /*! 
             *   number of sub-matrices in the row direction
             */
            unsigned int row_dimension;
            
            /*! 
             *   number of sub-matrices in the column direction
             */
            unsigned int column_dimension;
            
            /*!
             *   Vector of sparse matrices that form the sub-matrices of this object.
             *   The matrices are stored here in a logical 2-D array in a column major format, where all values from 
             *   the first column are followed by the values from the second column and so on.
             */
            std::vector<FESystem::Numerics::MatrixBase<ValType>*> sub_matrices;
            
            /*!
             *   stores the number of rows in each row-wise block
             */
            std::vector<FESystemUInt> rows_per_row_wise_block;
            
            /*!
             *   stores the number of columns in each column-wise block
             */
            std::vector<FESystemUInt> columns_per_column_wise_block;            
        };
    }
}



#endif // __fesystem_multi_block_matrix_h__

