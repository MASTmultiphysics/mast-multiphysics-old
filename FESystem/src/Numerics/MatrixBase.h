
/*
 *  MatrixBase.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/27/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

#ifndef __fesystem_matrix_base_h__
#define __fesystem_matrix_base_h__

// C++ includes
#include <memory>
#include <vector>
#include <map>

// FESystem includes
#include "Base/FESystemTypes.h"
#include "Base/FESystemExceptions.h"


namespace FESystem
{
	namespace Numerics
	{
        // Forward declerations
        template <typename ValType> class VectorBase;
        
        enum MatrixType
        {
            LOCAL_DENSE_MATRIX,
            LOCAL_SPARSE_MATRIX,
            MULTI_BLOCK_MATRIX
        };
        
        /*!
         *    This class provides the required interface for matrices. Different matrix classes derived from this provide the 
         *    specialized implementations of the methods. Note that all row and column numbers are zero-offset.
         */
		template <typename ValType>
		class MatrixBase
		{
		public:
            
            virtual ~MatrixBase(){}
            
            virtual void clear() = 0;
            
            /*!
             *   Returns the type of matrix
             */
            virtual FESystem::Numerics::MatrixType getType() const = 0;
            
            /*!
             *   Returns the number of rows and columns of this matrix
             */
            virtual const std::pair<FESystemUInt, FESystemUInt> getSize() const = 0;
            
            /*!
             *   resizes the matrix to the specified given number of rows, \p i, and columns, \p j, and zeros the contents
             */
            virtual void resize(FESystemUInt i, FESystemUInt j) = 0; 
			
            /*!
             *   Only copies the matrix values, without changing the sparsity pattern. This is meant only for sparse matrix data structures
             *   where for efficiency the values from one matrix of different sparsity pattern to another can be copied, as opposed to elementwise
             *   copy function calls. 
             */
            virtual void copyMatrixValues(const FESystem::Numerics::MatrixBase<ValType>& m) = 0;

            /*!
             *   Resizes and copies the contents of matrix \m to this matrix. 
             */
            virtual void copyMatrix(const FESystem::Numerics::MatrixBase<ValType>& m) = 0;

            /*!
             *   Resizes and copies the contents of the transpose of matrix \m to this matrix. 
             */
            virtual void copyMatrixTranspose(const FESystem::Numerics::MatrixBase<ValType>& m) = 0;

            /*!
             *   zeros this matrix
             */ 
            virtual void zero() = 0;
            
            /*!
             *   returns the value of A(i,j)
             */ 
			virtual ValType getVal(const FESystemUInt i, const FESystemUInt j) const = 0; 
			
            /*!
             *   A(i,i) = 1.0 for all i
             */ 
            virtual void setToIdentity() = 0;
            
            /*!
             *   A(i,j) =  val
             */ 
			virtual void setVal(const FESystemUInt i, 
                                const FESystemUInt j, 
                                const ValType& val) = 0; 
            
            /*!
             *   A(i,j) = A(i,j) + val
             */ 
			virtual void addVal(const FESystemUInt i, 
                                const FESystemUInt j, 
                                const ValType& val) = 0; 

            /*!
             *   Adds the values of the matrix to the specified indices 
             */
            virtual void addVal(const std::vector<FESystemUInt>& i, 
                                const std::vector<FESystemUInt>& j, 
                                const FESystem::Numerics::MatrixBase<ValType>& m) = 0; 

            /*!
             *   Returns the minimum value
             */
            virtual ValType getMinVal() const=0;
            
            /*!
             *   Returns the maximum value
             */
            virtual ValType getMaxVal() const=0;

            /*!
             *   copies the diagonal of the matrix to \p vec
             */
            virtual void getDiagonal(FESystem::Numerics::VectorBase<ValType>& vec) const = 0;

            
            /*!
             *   copies \p vec to the diagonal of the matrix
             */
            virtual void setDiagonal(const FESystem::Numerics::VectorBase<ValType>& vec) = 0;

            
            /*!
             *   Copies the contents of the vector \p vec into the columns col1-col2 for row number \p row_num. 
             *   Note that all numbers are zero offset. 
             */
            virtual void setRowVals(const FESystemUInt row_num, 
                                    const FESystemUInt col1,
                                    const FESystemUInt col2,
                                    const FESystem::Numerics::VectorBase<ValType>& vec) = 0; 
            
            /*!
             *   Copies the contents of the columns col1-col2 for row number \p row_num into vector \p vec
             *   Note that all numbers are zero offset. 
             */
            virtual void getRowVals(const FESystemUInt row_num, 
                                    const FESystemUInt col1,
                                    const FESystemUInt col2,
                                    FESystem::Numerics::VectorBase<ValType>& vec) const = 0; 

            
            /*!
             *   Copies the contents of the rows row1-row2 for column number \p col_num into vector \p vec
             *   Note that all numbers are zero offset. 
             */
            virtual void getColumnVals(const FESystemUInt col_num, 
                                       const FESystemUInt row1,
                                       const FESystemUInt row2,
                                       FESystem::Numerics::VectorBase<ValType>& vec) const = 0; 

            /*!
             *   Copies the contents of the vector \p vec into rows row1-row2 for column number \p col_num. 
             *   Note that all numbers are zero offset. 
             */
            virtual void setColumnVals(const FESystemUInt col_num, 
                                       const FESystemUInt row1,
                                       const FESystemUInt row2,
                                       const FESystem::Numerics::VectorBase<ValType>& vec) = 0; 
            
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
                                          FESystem::Numerics::MatrixBase<ValType>& m) const =0;

            /*!
             *   Copies the contents of the submatrix defined by the row and column indices given in \p rows and 
             *   \p cols. This is useful for extracting submatrices after application of Dirichlet boundary condition. 
             *   It is assumed that the rows and columns ids are in ascending order. 
             */
            virtual void getSubMatrixValsFromRowAndColumnIndices(const std::vector<FESystemUInt>& rows,
                                                                 const std::vector<FESystemUInt>& cols,
                                                                 const std::map<FESystemUInt, FESystemUInt>& old_to_new_id_map,
                                                                 FESystem::Numerics::MatrixBase<ValType>& mat) const= 0;
            
            /*!
             *   Zeros the contents of the submatrix (m_row1:m_row2, m_col1:m_col2) 
             */
            virtual void zeroSubMatrixVals(const FESystemUInt row1, const FESystemUInt row2,
                                           const FESystemUInt col1, const FESystemUInt col2)=0;

            
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
                                          const FESystem::Numerics::MatrixBase<ValType>& m)=0;

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
                                          const ValType f, const FESystem::Numerics::MatrixBase<ValType>& m)=0;

            /*!
             *   Performs the dot product of the \p col_num^th column with the given vector \p vec. 
             */
            virtual ValType dotProductWithColumn(const FESystemUInt col_num, 
                                                 FESystem::Numerics::VectorBase<ValType>& vec) const=0; 

            /*!
             *    Scales this matrix with factor \p t.
             */
			virtual  void scale(const ValType& t)=0;

			/*!
             *    Scales the \p i^th row with the factor \p t.
             */
			virtual  void scaleRow(const FESystemInt i, const ValType& t)=0;
            
			/*!
             *    Scales the \p i^th columns with the factor \p t.
             */
            virtual  void scaleColumn(const FESystemInt i, const ValType& t)=0;
            
            /*!
             *   this->(row1, col1:col2) = f1*this->(row1, col1:col2) + f2*(row2,col1:col2)
             */
            virtual void addScaledSubRow(FESystemUInt row1, ValType f1, FESystemUInt row2, ValType f2, FESystemUInt col1, FESystemUInt col2)=0;
            
			/*!
             *    Creates the orthogonal projector onto the vector \p v and stores it in this matrix. This is given as 
             *    \[$ P = v v^T \]$
             */
            virtual void createOrthogonalProjectorOntoVector(const FESystem::Numerics::VectorBase<ValType>& v)=0;
            
			/*!
             *    Creates the orthogonal projector along the vector \p v and stores it in this matrix. This is given as 
             *    \[$ P = I - v v^T \]$
             */
            virtual void createOrthogonalProjectorAlongVector(const FESystem::Numerics::VectorBase<ValType>& v)=0;
            
			/*!
             *    Creates the orthogonal reflector such that multiplication of the vector \p v with the reflector rotates it to
             *    a vector of same length along the e_1 axis. Since the rotated vector could be chosen to lie on either the 
             *    positive or the negative e_1, the larger movement of the following two is chosen as the 'along' vector (|v|e_1 - v) or (-|v|e_1-v). 
             *    \[$ P = I - 2 v1 v1^T \]$ where \[$ v1 =  +/- |v|e_1 - v \]$
             */
            virtual void createOrthogonalReflectorForVector(const FESystem::Numerics::VectorBase<ValType>& v)=0;

			/*!
             *    shifts the diagonal of this matrix by the vector v.
             */
            virtual void shiftDiagonal(const ValType& v)=0;
            
			/*!
             *    Returns \f$ v^T A v \f$
             */
			virtual ValType getRayleighQuotient(const FESystem::Numerics::VectorBase<ValType>& v) =0;

            /*! 
             *   This method scales each row to a unit vector
             */ 
            virtual void normalizeRowBasis() = 0;
            
            /*! 
             *   This method scales each column to a unit vector
             */ 
            virtual void normalizeColumnBasis() = 0;
            
            
            /*!
             *   Calculates the matrix 1-norm. 
             */            
            virtual typename RealOperationType(ValType) getL1Norm() const = 0;
            
            /*!
             *   Calculates the matrix Frobenius
             */            
            virtual typename RealOperationType(ValType) getFrobeniusNorm() const = 0;

            /*!
             *   Calculates the matrix determinant
             */            
            virtual ValType getDeterminant() const = 0;
            
            /*!
             *   Calculates the matrix inverse and returns it in \p mat
             */            
            virtual void getInverse(FESystem::Numerics::MatrixBase<ValType>& mat) const = 0;
            
            /*!
             *   Calculates the eigenvalue of the matrix diagonal block \f$ \begin{array}{cc} m_{i,i} & m_{i,i+1}\\ m_{i+1,i} & m_{i+1,i+1} \end{array} \f$. 
             *   This method returns real values, and should be used only if the uses is certain about the symmetry of the block. The value are returned in the \p val1 and \p val2. 
             */
            template <typename ValTypeEig> void getDiagonalSymmetricBlockEigenvalue(FESystemUInt i, ValTypeEig& val1, ValTypeEig& val2) const;
            
            /*!
             *   Calculates the eigenvalue of the matrix diagonal block \f$ \begin{array}{cc} m_{i,i} & m_{i,i+1}\\ m_{i+1,i} & m_{i+1,i+1} \end{array} \f$. 
             *   This method returns complex values, and should ideally be used when the block is likely to be asymmetric. This will also work for 
             */
            template <typename ValTypeEig> void getDiagonalBlockEigenvalue(FESystemUInt i, ValTypeEig& val1, ValTypeEig& val2) const;

            /*!
             *   \brief \f$ \{x\}_{res} = [A] \{x\}_t \f$
             */
            virtual  void rightVectorMultiply(const FESystem::Numerics::VectorBase<ValType>& t,
                                              FESystem::Numerics::VectorBase<ValType>& res) const=0;
            /*!
             *   \brief \f$ \{x\}_{res} = [A] \{x\}_t \f$, multiply with the submatrix (row1:row2, col1:col2)
             */
            virtual  void rightSubVectorMultiply(FESystemUInt row1, FESystemUInt row2, FESystemUInt col1, FESystemUInt col2, 
                                                 FESystemUInt v_row1, FESystemUInt v_row2, FESystemUInt res_row1, FESystemUInt res_row2,
                                                 const FESystem::Numerics::VectorBase<ValType>& t, FESystem::Numerics::VectorBase<ValType>& res) const=0;
            
            /*!
             *   \brief \f$ \{x\}_{res}^{T} = \{x\}_t^T  [A] \f$
             */
            virtual  void leftVectorMultiply(const FESystem::Numerics::VectorBase<ValType>& t,
                                             FESystem::Numerics::VectorBase<ValType>& res) const=0;
            

            virtual ValType multiplySubVectorWithSubRow(const FESystemUInt row, const FESystemUInt col1, const FESystemUInt col2,
                                                        const FESystemUInt vec_el1, const FESystemUInt vec_el2, const FESystem::Numerics::VectorBase<ValType>& vec) const=0;

            virtual ValType multiplySubVectorWithSubColumn(const FESystemUInt col, const FESystemUInt row1, const FESystemUInt row2,
                                                           const FESystemUInt vec_el1, const FESystemUInt vec_el2, const FESystem::Numerics::VectorBase<ValType>& vec) const=0;

            /*!
             *   \brief \f$ [X]_{res} = [A] [X]_t \f$
             */
            virtual void matrixRightMultiply (ValType f, const MatrixBase<ValType>& t,
                                              MatrixBase<ValType>& res) const=0;
                        
            /*!
             *   \brief \f$ [X]_{res} = [A] [X](row1:row2,col1:col2)_t \f$
             */
            virtual void matrixRightMultiplySubMatrix(FESystemUInt row1, FESystemUInt row2, FESystemUInt col1, FESystemUInt col2,
                                                      ValType f, const MatrixBase<ValType>& t, MatrixBase<ValType>& res) const=0;
            
            /*!
             *   \brief \f$ [X]_{res} = [A]^T [X]_t \f$
             */
            virtual void matrixTransposeRightMultiply (ValType f, const MatrixBase<ValType>& t,
                                                       MatrixBase<ValType>& res) const=0;
            
            /*!
             *   \brief \f$ [X]_{res} = [A] [X]_t^T \f$
             */            
            virtual void matrixRightMultiplyTranspose (ValType f, const MatrixBase<ValType>& t,
                                                      MatrixBase<ValType>& res) const=0;
            
            /*!
             *    Initializes the LU factored matrices based on the current matrix. The matrix L^{-1} is returned in l_mat and the matrix U is returned in u_mat.
             */  
            virtual void initializeLUFactoredMatrices(FESystem::Numerics::MatrixBase<ValType>& l_mat, FESystem::Numerics::MatrixBase<ValType>& u_mat) const = 0;
            
            /*!
             *   \f$ A = f * A_t \f$
             */
			virtual void add (ValType f, const MatrixBase<ValType>& t)=0;
            
            /*!
             *   Returns a constant pointer to the vector of values of this matrix. Should be used with care. 
             */
            virtual const ValType* getMatrixValues() const = 0;
            
            /*!
             *   Returns a pointer to the vector of values of this matrix. Should be used with care. 
             */
            virtual ValType* getMatrixValues() = 0;
            
            /*!
             *   Write a formatted output to the stream \p out
             */
            virtual void write(std::ostream& out) const = 0;
            
            /*!
             *   Write a formatted output to the stream \p out
             */
            virtual void writeDetailed(std::ostream& out) const = 0;
            			
		protected:
			
            
		};
        
        DeclareException0(MatrixNotInitialized, 
                          << "Matrix Not Initialized Before Usage\n");

        
        DeclareException0(MatrixMustBeSquareForOperation, 
                          << "Operation is valid only for square matrix. This is not a square matrix.\n");

        DeclareException2(RowNumberSizeMismatch, 
                          FESystemUInt, FESystemUInt, 
                          << "Row Number Exceedes Matrix Size\n"
                          << "Matrix Row Dimension: " << Arg1 << "\n"
                          << "Row Number: " << Arg2); 
        
        DeclareException2(ColumnNumberSizeMismatch, 
                          FESystemUInt, FESystemUInt, 
                          << "Column Number Exceedes Matrix Size\n"
                          << "Matrix Column Dimension: " << Arg1 << "\n"
                          << "Column Number: " << Arg2); 
        
        DeclareException3(MatrixVectorSizeMismatch, 
                          FESystemUInt, FESystemUInt, FESystemUInt, 
                          << "Matrix and Vector Size Mismatch\n"
                          << "Matrix Dimension: " << Arg1 << ", " << Arg2 <<"\n"
                          << "Vector Dimension: " << Arg3); 
        
        DeclareException4(MatrixMultiplicationSizeMismatch, 
                          FESystemUInt, FESystemUInt, FESystemUInt, FESystemUInt, 
                          << "Matrix and Multiplication Size Mismatch\n"
                          << "Matrix-1 Dimension: " << Arg1 << ", " << Arg2 <<"\n"
                          << "Matrix-2 Dimension: " << Arg3 << ", " << Arg4); 
        
        DeclareException4(ElementLocationExceedesMatrixSize, 
                          FESystemUInt, FESystemUInt, FESystemUInt, FESystemUInt, 
                          << "Element Location Exceedes Matrix Size\n"
                          << "Matrix Dimension: " << Arg1 << ", " << Arg2 <<"\n"
                          << "Element Location: " << Arg3 << ", " << Arg4);

        DeclareException4(MatrixSizeMismatch, 
                          FESystemUInt, FESystemUInt, FESystemUInt, FESystemUInt, 
                          << "Matrix Size Does Not Match\n"
                          << "Should Be: " << Arg1 << ", " << Arg2 <<"\n"
                          << "Found: " << Arg3 << ", " << Arg4);

        
        template <typename ValType>
        std::auto_ptr<FESystem::Numerics::MatrixBase<ValType> > 
        MatrixCreate(FESystem::Numerics::MatrixType mat_type);
	}
}




template <typename ValType>
template <typename ValTypeEig>
void inline FESystem::Numerics::MatrixBase<ValType>::getDiagonalSymmetricBlockEigenvalue(FESystemUInt i, ValTypeEig& val1, ValTypeEig& val2) const
{   
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    
    FESystemAssert4( ((i < s.first) && (i+1 < s.second)),
                    FESystem::Numerics::ElementLocationExceedesMatrixSize,
                    s.first, s.second, i, i+1);
    
    // get the four entries of the block
    ValTypeEig a = FESystem::Base::real(this->getVal(i,i)),
    d = FESystem::Base::real(this->getVal(i+1,i+1)), val3, val4; // for symmetric matrices, diagonal values are always, and so are the eigenvalues
    ValType b = this->getVal(i,i+1), c=this->getVal(i+1,i);
    // the polynomial is created from [a b; c d] and the A, B, C terms of the quadratic polynomial A l^2 + B l + C = 0 are
    // A = 1, B = -(a+d), C = (ad-bc)
    val3 = sqrt(pow(a+d,2)-4.0*(a*d-b*c)); // val3 = sqrt((a+d)^2-4(ad-bc))
    
    val1 = (a+d+val3)*0.5;
    val2 = (a+d-val3)*0.5;
}



template <typename ValType>
template <typename ValTypeEig>
void inline FESystem::Numerics::MatrixBase<ValType>::getDiagonalBlockEigenvalue(FESystemUInt i, ValTypeEig& val1, ValTypeEig& val2) const
{   
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    
    FESystemAssert4( ((i < s.first) && (i+1 < s.second)),
                    FESystem::Numerics::ElementLocationExceedesMatrixSize,
                    s.first, s.second, i, i+1);
    
    // get the four entries of the block
    ValTypeEig a = this->getVal(i,i), b = this->getVal(i,i+1), c=this->getVal(i+1,i), d = this->getVal(i+1,i+1);
    ValTypeEig val3 = pow(a+d,2), val4 = (a*d-b*c);
    val4 *= -4.0; 
    val3 += val4; val3 = sqrt(val3);
    // the polynomial is created from [a b; c d] and the A, B, C terms of the quadratic polynomial A l^2 + B l + C = 0 are
    // A = 1, B = -(a+d), C = (ad-bc)
    val1 = val3; val1 += (a+d); val1*= 0.5;
    val2 = -val3; val2 += (a+d); val2*= 0.5;
}




#endif // __fesystem_Matrix_base_h__
