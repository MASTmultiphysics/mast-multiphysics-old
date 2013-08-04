//
//  structural_basis_matrix.h
//  MAST
//
//  Created by Manav Bhatia on 7/28/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef MAST_structural_basis_matrix_h
#define MAST_structural_basis_matrix_h

// libmesh includes
#include "libmesh/shell_matrix.h"

template <typename T>
class BasisMatrix: public ShellMatrix<T>
{
public:
    BasisMatrix():
    ShellMatrix<T>()
    { }
    
    virtual ~BasisMatrix()
    { }
    
    /**
     * @returns \p m, the row-dimension of the matrix where the marix is
     * \f$ M \times N \f$.
     */
    virtual numeric_index_type m () const;
    
    /**
     * @returns \p n, the column-dimension of the matrix where the marix
     * is \f$ M \times N \f$.
     */
    virtual numeric_index_type n () const;
    
    /**
     * Multiplies the matrix with \p arg and stores the result in \p
     * dest.
     */
    virtual void vector_mult (NumericVector<T>& dest,
                              const NumericVector<T>& arg) const;
    
    /**
     * Multiplies the transpose of matrix with \p arg and stores the
     * result in \p dest.
     */
    virtual void vector_mult_transpose (NumericVector<T>& dest,
                                        const NumericVector<T>& arg) const;

    /**
     * Multiplies the matrix with \p arg and adds the result to \p dest.
     */
    virtual void vector_mult_add (NumericVector<T>& dest,
                                  const NumericVector<T>& arg) const;
    
    /**
     * Copies the diagonal part of the matrix into \p dest.
     */
    virtual void get_diagonal (NumericVector<T>& dest) const;
    
    /*!
     *  Returns the vector that defines the \p i^th basis vector
     */
    
    virtual NumericVector<T>& basis(unsigned int i);
};


#endif
