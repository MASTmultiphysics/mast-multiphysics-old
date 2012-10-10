//
//  SparseMatrix.h
//  FESystem
//
//  Created by Manav Bhatia on 4/26/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem_sparse_matrix_h__
#define __fesystem_sparse_matrix_h__


// FESystem includes
#include "Numerics/LocalMatrixBase.h"


namespace FESystem
{
	namespace Numerics
	{
        // Forward declerations
        class SparsityPattern;
        
        /*!
         *    This specializes the MatrixBase class for sparse matrices. The sparsit pattern is used to initialize the matrix, and it controls the 
         *    number of nonzero entries per row.
         */
		template <typename ValType>
		class SparseMatrix: public FESystem::Numerics::LocalMatrixBase<ValType>
		{
		public:
            /*!
             *    Constructor
             */
            SparseMatrix();
            
            
            virtual ~SparseMatrix();

            /*!
             *   Returns the type of matrix
             */
            virtual FESystem::Numerics::MatrixType getType() const;
                        
            /*!
             *   Returns a constant reference to the sparsity pattern data structure
             */
            const FESystem::Numerics::SparsityPattern& getSparsityPattern() const;

            /*!
             *   resizes the matrix to the specified given number of rows, \p i, and columns, \p j, and zeros the contents
             */
            virtual void resize(FESystemUInt i, FESystemUInt j); 
            
            
            /*!
             *    Initializes itself for the specified sparsity pattern
             */ 
            virtual void resize(const FESystem::Numerics::SparsityPattern& s_pattern); 
            
			
            /*!
             *   Only copies the matrix values, and does not change the sparsity patter or size
             */
            virtual void copyMatrixValues(const FESystem::Numerics::MatrixBase<ValType>& m);

            
            /*!
             *   Resizes and copies the contents of matrix \m to this matrix. 
             */
            virtual void copyMatrix(const FESystem::Numerics::MatrixBase<ValType>& m);
            
            /*!
             *   Resizes and copies the contents of matrix \m to this matrix.
             */
            virtual void copyRealMatrix(const FESystem::Numerics::MatrixBase<typename RealOperationType(ValType)>& m);

            /*!
             *   Resizes and copies the contents of the transpose of matrix \m to this matrix. 
             */
            virtual void copyMatrixTranspose(const FESystem::Numerics::MatrixBase<ValType>& m);
            
            /*!
             *   returns the value of A(i,j)
             */ 
			virtual ValType getVal(const FESystemUInt i, const FESystemUInt j) const; 
			
            /*!
             *   A(i,i) = 1.0 for all i
             */ 
            virtual void setToIdentity();
            
            /*!
             *   A(i,j) =  val
             */ 
			virtual void setVal(const FESystemUInt i, 
                                const FESystemUInt j, 
                                const ValType& val); 
            
            /*!
             *   A(i,j) = A(i,j) + val
             */ 
			virtual void addVal(const FESystemUInt i, 
                                const FESystemUInt j, 
                                const ValType& val); 
            
            /*!
             *   Adds the values of the matrix to the specified indices 
             */
            virtual void addVal(const std::vector<FESystemUInt>& i, 
                                const std::vector<FESystemUInt>& j, 
                                const FESystem::Numerics::MatrixBase<ValType>& m); 

            /*!
             *   copies the diagonal of the matrix to \p vec
             */
            virtual void getDiagonal(FESystem::Numerics::VectorBase<ValType>& vec) const;
            
            
            /*!
             *   copies \p vec to the diagonal of the matrix
             */
            virtual void setDiagonal(const FESystem::Numerics::VectorBase<ValType>& vec);
            
            
            /*!
             *   Copies the contents of the vector \p vec into the columns col1-col2 for row number \p row_num. 
             *   Note that all numbers are zero offset. 
             */
            virtual void setRowVals(const FESystemUInt row_num, 
                                    const FESystemUInt col1,
                                    const FESystemUInt col2,
                                    const FESystem::Numerics::VectorBase<ValType>& vec); 
            
            /*!
             *   Copies the contents of the columns col1-col2 for row number \p row_num into vector \p vec
             *   Note that all numbers are zero offset. 
             */
            virtual void getRowVals(const FESystemUInt row_num, 
                                    const FESystemUInt col1,
                                    const FESystemUInt col2,
                                    FESystem::Numerics::VectorBase<ValType>& vec) const; 
            
            
            /*!
             *   Copies the contents of the rows row1-row2 for column number \p col_num into vector \p vec
             *   Note that all numbers are zero offset. 
             */
            virtual void getColumnVals(const FESystemUInt col_num, 
                                       const FESystemUInt row1,
                                       const FESystemUInt row2,
                                       FESystem::Numerics::VectorBase<ValType>& vec) const; 
            
            /*!
             *   Copies the contents of the vector \p vec into rows row1-row2 for column number \p col_num. 
             *   Note that all numbers are zero offset. 
             */
            virtual void setColumnVals(const FESystemUInt col_num, 
                                       const FESystemUInt row1,
                                       const FESystemUInt row2,
                                       const FESystem::Numerics::VectorBase<ValType>& vec); 
            
            /*!
             *   Copies the contents of the submatrix (row1:row2, col1:col2) of this matrix into the submatrix 
             *   (m_row1:m_row2, m_col1:m_col2) of the given matrix \p m
             *   Note that all numbers are zero offset. 
             */
            virtual void getSubMatrixVals(const FESystemUInt row1,
                                          const FESystemUInt row2,
                                          const FESystemUInt col1,
                                          const FESystemUInt col2,
                                          const FESystemUInt m_row1,
                                          const FESystemUInt m_row2,
                                          const FESystemUInt m_col1,
                                          const FESystemUInt m_col2,
                                          FESystem::Numerics::MatrixBase<ValType>& m) const ;
            
            /*!
             *   Copies the contents of the submatrix defined by the row and column indices given in \p rows and 
             *   \p cols. This is useful for extracting submatrices after application of Dirichlet boundary condition. 
             *   It is assumed that the rows and columns ids are in ascending order. 
             */
            virtual void getSubMatrixValsFromRowAndColumnIndices(const std::vector<FESystemUInt>& rows,
                                                                 const std::vector<FESystemUInt>& cols,
                                                                 const std::map<FESystemUInt, FESystemUInt>& old_to_new_id_map,
                                                                 FESystem::Numerics::MatrixBase<ValType>& mat) const;
            
            /*!
             *   Zeros the contents of the submatrix (m_row1:m_row2, m_col1:m_col2) 
             */
            virtual void zeroSubMatrixVals(const FESystemUInt row1, const FESystemUInt row2,
                                           const FESystemUInt col1, const FESystemUInt col2);

            /*!
             *   Copies the contents of the submatrix (m_row1:m_row2, m_col1:m_col2) of the given matrix \p m into 
             *   (row1:row2, col1:col2) of this matrix.  
             *   Note that all numbers are zero offset. 
             */
            virtual void setSubMatrixVals(const FESystemUInt row1,
                                          const FESystemUInt row2,
                                          const FESystemUInt col1,
                                          const FESystemUInt col2,
                                          const FESystemUInt m_row1,
                                          const FESystemUInt m_row2,
                                          const FESystemUInt m_col1,
                                          const FESystemUInt m_col2,
                                          const FESystem::Numerics::MatrixBase<ValType>& m);

            /*!
             *   Adds the contents of the submatrix (m_row1:m_row2, m_col1:m_col2) of the given matrix \p m into 
             *   (row1:row2, col1:col2) of this matrix.  
             *   Note that all numbers are zero offset. 
             */
            virtual void addSubMatrixVals(const FESystemUInt row1,
                                          const FESystemUInt row2,
                                          const FESystemUInt col1,
                                          const FESystemUInt col2,
                                          const FESystemUInt m_row1,
                                          const FESystemUInt m_row2,
                                          const FESystemUInt m_col1,
                                          const FESystemUInt m_col2,
                                          const ValType f, const FESystem::Numerics::MatrixBase<ValType>& m);

            /*!
             *   Performs the dot product of the \p col_num^th column with the given vector \p vec. 
             */
            virtual ValType dotProductWithColumn(const FESystemUInt col_num, 
                                                 FESystem::Numerics::VectorBase<ValType>& vec) const; 
            
			/*!
             *    Scales the \p i^th row with the factor \p t.
             */
			virtual  void scaleRow(const FESystemInt i, const ValType& t);
            
			/*!
             *    Scales the \p i^th columns with the factor \p t.
             */
            virtual  void scaleColumn(const FESystemInt i, const ValType& t);
            
            /*!
             *   this->(row1, :) = f1*this->(row1, :) + f2*(row2,:)
             */
            virtual void addScaledSubRow(FESystemUInt row1, ValType f1, FESystemUInt row2, ValType f2, FESystemUInt col1, FESystemUInt col2);
            
			/*!
             *    Creates the orthogonal projector onto the vector \p v and stores it in this matrix. This is given as 
             *    \[$ P = v v^T \]$
             */
            virtual void createOrthogonalProjectorOntoVector(const FESystem::Numerics::VectorBase<ValType>& v);
            
			/*!
             *    Creates the orthogonal projector along the vector \p v and stores it in this matrix. This is given as 
             *    \[$ P = I - v v^T \]$
             */
            virtual void createOrthogonalProjectorAlongVector(const FESystem::Numerics::VectorBase<ValType>& v);
            
			/*!
             *    Creates the orthogonal reflector such that multiplication of the vector \p v with the reflector rotates it to
             *    a vector of same length along the e_1 axis. Since the rotated vector could be chosen to lie on either the 
             *    positive or the negative e_1, the larger movement of the following two is chosen as the 'along' vector (|v|e_1 - v) or (-|v|e_1-v). 
             *    \[$ P = I - 2 v1 v1^T \]$ where \[$ v1 =  +/- |v|e_1 - v \]$
             */
            virtual void createOrthogonalReflectorForVector(const FESystem::Numerics::VectorBase<ValType>& v);
            
			/*!
             *    shifts the diagonal of this matrix by the vector v.
             */
            virtual void shiftDiagonal(const ValType& v);
                                    
            /*!
             *   \brief \f$ \{x\}_{res} = [A] \{x\}_t \f$
             */
            virtual  void rightVectorMultiply(const FESystem::Numerics::VectorBase<ValType>& t,
                                              FESystem::Numerics::VectorBase<ValType>& res) const;
            
            /*!
             *   \brief \f$ \{x\}_{res} = [A] \{x\}_t \f$, multiply with the submatrix (row1:row2, col1:col2)
             */
            virtual  void rightSubVectorMultiply(FESystemUInt row1, FESystemUInt row2, FESystemUInt col1, FESystemUInt col2, 
                                                 FESystemUInt v_row1, FESystemUInt v_row2, FESystemUInt res_row1, FESystemUInt res_row2,
                                                 const FESystem::Numerics::VectorBase<ValType>& t, FESystem::Numerics::VectorBase<ValType>& res) const;

            /*!
             *   \brief \f$ \{x\}_{res}^{T} = \{x\}_t^T  [A] \f$
             */
            virtual  void leftVectorMultiply(const FESystem::Numerics::VectorBase<ValType>& t,
                                             FESystem::Numerics::VectorBase<ValType>& res) const;
            
            
            virtual ValType multiplySubVectorWithSubRow(const FESystemUInt row, const FESystemUInt col1, const FESystemUInt col2,
                                                        const FESystemUInt vec_el1, const FESystemUInt vec_el2, const FESystem::Numerics::VectorBase<ValType>& vec) const;
            

            virtual ValType multiplySubVectorWithSubColumn(const FESystemUInt col, const FESystemUInt row1, const FESystemUInt row2,
                                                           const FESystemUInt vec_el1, const FESystemUInt vec_el2, const FESystem::Numerics::VectorBase<ValType>& vec) const;

            /*!
             *   \brief \f$ [X]_{res} = [A] [X]_t \f$
             */
            virtual void matrixRightMultiply (ValType f, const MatrixBase<ValType>& t,
                                              MatrixBase<ValType>& res) const;
            
            
            /*!
             *   \brief \f$ [X]_{res} = [A] [X](row1:row2,col1:col2)_t \f$
             */
            virtual void matrixRightMultiplySubMatrix(FESystemUInt row1, FESystemUInt row2, FESystemUInt col1, FESystemUInt col2,
                                                      ValType f, const MatrixBase<ValType>& t, MatrixBase<ValType>& res) const;
            
            /*!
             *   \brief \f$ [X]_{res} = [A]^T [X]_t \f$
             */
            virtual void matrixTransposeRightMultiply (ValType f, const MatrixBase<ValType>& t,
                                                       MatrixBase<ValType>& res) const;
            
            /*!
             *   \brief \f$ [X]_{res} = [A] [X]_t^T \f$
             */            
            virtual void matrixRightMultiplyTranspose (ValType f, const MatrixBase<ValType>& t,
                                                       MatrixBase<ValType>& res) const;
            
            /*!
             *   \f$ A = f * A_t \f$
             */
			virtual void add (ValType f, const MatrixBase<ValType>& t);
            
            /*!
             *    Initializes the LU factored matrices based on the current matrix
             */  
            virtual void initializeLUFactoredMatrices(FESystem::Numerics::MatrixBase<ValType>& l_mat, FESystem::Numerics::MatrixBase<ValType>& u_mat) const;

		protected:
			            
            /*!
             *    Sparsity pattern for this matrix
             */
            const FESystem::Numerics::SparsityPattern* sparsity_pattern;
            
		};
	}
}




#endif // __fesystem_sparse_matrix_h__
