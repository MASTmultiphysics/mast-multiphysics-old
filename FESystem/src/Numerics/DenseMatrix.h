/*
 *  DenseMatrix.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/5/10.
 *  Copyright 2010 . All rights reserved.
 *
 */

#ifndef __fesystem_dense_matrix_h__
#define __fesystem_dense_matrix_h__

// C++ includes
#include <vector>
#include <set>

// FESystem includes
#include "Numerics/LocalMatrixBase.h"

namespace FESystem
{
	namespace Numerics
	{
		template <typename ValType>
		class DenseMatrix: public FESystem::Numerics::LocalMatrixBase<ValType>
		{
		public:
			DenseMatrix();
			
			virtual ~DenseMatrix();
			
            virtual FESystem::Numerics::MatrixType getType() const;
            
            virtual const std::pair<FESystemUInt, FESystemUInt> getSize() const;
            
            void resize(FESystemUInt i, FESystemUInt j); 
            
            void resize(ValType* v, FESystemUInt i, FESystemUInt j); 
			
            virtual void copyMatrixValues(const FESystem::Numerics::MatrixBase<ValType>& m);

            virtual void copyMatrix(const FESystem::Numerics::MatrixBase<ValType>& m);
            
            /*!
             *   Resizes and copies the contents of matrix \m to this matrix.
             */
            virtual void copyRealMatrix(const FESystem::Numerics::MatrixBase<typename RealOperationType(ValType)>& m);

            virtual void copyMatrixTranspose(const FESystem::Numerics::MatrixBase<ValType>& m);

            void zero();
            
            virtual void setToIdentity();

			virtual ValType getVal(const FESystemUInt i, const FESystemUInt j) const; 
			
			virtual void setVal(const FESystemUInt i,
                                const FESystemUInt j, 
                                const ValType& val) ; 
            
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

            virtual void getDiagonal(FESystem::Numerics::VectorBase<ValType>& vec) const;
                        
            virtual void setDiagonal(const FESystem::Numerics::VectorBase<ValType>& vec);


            virtual void setRowVals(const FESystemUInt row_num, 
                                    const FESystemUInt col1,
                                    const FESystemUInt col2,
                                    const FESystem::Numerics::VectorBase<ValType>& vec); 
            
            virtual void getRowVals(const FESystemUInt row_num, 
                                    const FESystemUInt col1,
                                    const FESystemUInt col2,
                                    FESystem::Numerics::VectorBase<ValType>& vec) const; 
            
            virtual void getColumnVals(const FESystemUInt col_num, 
                                       const FESystemUInt row1,
                                       const FESystemUInt row2,
                                       FESystem::Numerics::VectorBase<ValType>& vec) const; 

            virtual void setColumnVals(const FESystemUInt col_num, 
                                       const FESystemUInt row1,
                                       const FESystemUInt row2,
                                       const FESystem::Numerics::VectorBase<ValType>& vec); 

            virtual void getSubMatrixVals(const FESystemUInt row1,
                                          const FESystemUInt row2,
                                          const FESystemUInt col1,
                                          const FESystemUInt col2,
                                          const FESystemUInt m_row1,
                                          const FESystemUInt m_row2,
                                          const FESystemUInt m_col1,
                                          const FESystemUInt m_col2,
                                          FESystem::Numerics::MatrixBase<ValType>& m) const;

            /*!
             *   Zeros the contents of the submatrix (m_row1:m_row2, m_col1:m_col2) 
             */
            virtual void zeroSubMatrixVals(const FESystemUInt row1, const FESystemUInt row2,
                                           const FESystemUInt col1, const FESystemUInt col2);

            virtual void setSubMatrixVals(const FESystemUInt row1,
                                          const FESystemUInt row2,
                                          const FESystemUInt col1,
                                          const FESystemUInt col2,
                                          const FESystemUInt m_row1,
                                          const FESystemUInt m_row2,
                                          const FESystemUInt m_col1,
                                          const FESystemUInt m_col2,
                                          const FESystem::Numerics::MatrixBase<ValType>& m);

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
             *   Copies the contents of the submatrix defined by the row and column indices given in \p rows and 
             *   \p cols. This is useful for extracting submatrices after application of Dirichlet boundary condition. 
             *   It is assumed that the row and column ids are in ascending order. 
             */
            virtual void getSubMatrixValsFromRowAndColumnIndices(const std::vector<FESystemUInt>& rows,
                                                                 const std::vector<FESystemUInt>& cols,
                                                                 const std::map<FESystemUInt, FESystemUInt>& old_to_new_id_map,
                                                                 FESystem::Numerics::MatrixBase<ValType>& mat) const;

			virtual  void scaleRow(const FESystemInt row_num, const ValType& t);
            
            virtual  void scaleColumn(const FESystemInt col_num, const ValType& t);
            
            /*!
             *   this->(row1, :) = f1*this->(row1, :) + f2*(row2,:)
             */
            virtual void addScaledSubRow(FESystemUInt row1, ValType f1, FESystemUInt row2, ValType f2, FESystemUInt col1, FESystemUInt col2);

            
            virtual void shiftDiagonal(const ValType& v);

            virtual void createOrthogonalProjectorOntoVector(const FESystem::Numerics::VectorBase<ValType>& v);
            
            virtual void createOrthogonalProjectorAlongVector(const FESystem::Numerics::VectorBase<ValType>& v);
            
            virtual void createOrthogonalReflectorForVector(const FESystem::Numerics::VectorBase<ValType>& v);

            virtual ValType dotProductWithColumn(const FESystemUInt col_num, 
                                                 FESystem::Numerics::VectorBase<ValType>& vec) const; 

			virtual void rightVectorMultiply(const FESystem::Numerics::VectorBase<ValType>& t,
                                             FESystem::Numerics::VectorBase<ValType>& res) const;
            
            /*!
             *   \brief \f$ \{x\}_{res} = [A] \{x\}_t \f$, multiply with the submatrix (row1:row2, col1:col2)
             */
            virtual  void rightSubVectorMultiply(FESystemUInt row1, FESystemUInt row2, FESystemUInt col1, FESystemUInt col2, 
                                                 FESystemUInt v_row1, FESystemUInt v_row2, FESystemUInt res_row1, FESystemUInt res_row2,
                                                 const FESystem::Numerics::VectorBase<ValType>& t, FESystem::Numerics::VectorBase<ValType>& res) const;

            virtual void leftVectorMultiply(const FESystem::Numerics::VectorBase<ValType>& t,
                                            FESystem::Numerics::VectorBase<ValType>& res) const;
            
            virtual ValType multiplySubVectorWithSubRow(const FESystemUInt row, const FESystemUInt col1, const FESystemUInt col2,
                                                        const FESystemUInt vec_el1, const FESystemUInt vec_el2, const FESystem::Numerics::VectorBase<ValType>& vec) const;
            
            
            virtual ValType multiplySubVectorWithSubColumn(const FESystemUInt col, const FESystemUInt row1, const FESystemUInt row2,
                                                           const FESystemUInt vec_el1, const FESystemUInt vec_el2, const FESystem::Numerics::VectorBase<ValType>& vec) const;

            virtual void matrixRightMultiply (ValType f, 
                                              const MatrixBase<ValType>& t,
                                              MatrixBase<ValType>& res) const;
            
            /*!
             *   \brief \f$ [X]_{res} = [A] [X]_t \f$
             */
            virtual void matrixTransposeRightMultiply (ValType f, 
                                                       const MatrixBase<ValType>& t,
                                                       MatrixBase<ValType>& res) const;
            
            /*!
             *   \brief \f$ [X]_{res} = [A] [X](row1:row2,col1:col2)_t \f$
             */
            virtual void matrixRightMultiplySubMatrix(FESystemUInt row1, FESystemUInt row2, FESystemUInt col1, FESystemUInt col2,
                                                      ValType f, const MatrixBase<ValType>& t, MatrixBase<ValType>& res) const;
            
            virtual void matrixRightMultiplyTranspose (ValType f, 
                                                       const MatrixBase<ValType>& t,
                                                       MatrixBase<ValType>& res) const;
            
			virtual void add (ValType f, const MatrixBase<ValType>& t);
            
            /*!
             *    Initializes the LU factored matrices based on the current matrix. 
             */  
            virtual void initializeLUFactoredMatrices(FESystem::Numerics::MatrixBase<ValType>& l_mat, FESystem::Numerics::MatrixBase<ValType>& u_mat) const;
            
		protected:
			

            /*!
             *   Returns the value of the cofactor for a submatrix formed by the first column to the last, and the 
             *   set of rows included in \p rows . The number of columns and number of retained rows for this submatrix
             *   should form a submatrix. 
             */
            ValType getCofactor(FESystemUInt first_column, const std::set<FESystemUInt>& rows) const;

                        
		};
                
	}
}

#endif

