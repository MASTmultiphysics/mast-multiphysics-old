//
//  structural_basis_matrix.h
//  MAST
//
//  Created by Manav Bhatia on 7/28/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef MAST_structural_basis_matrix_h
#define MAST_structural_basis_matrix_h

// C++ includes
#include <vector>

// libmesh includes
#include "libmesh/shell_matrix.h"
#include "libmesh/numeric_vector.h"


template <typename T>
class BasisMatrix: public ShellMatrix<T>
{
public:
    BasisMatrix(const Parallel::Communicator &comm_in):
    ShellMatrix<T>(comm_in)
    { }
    
    virtual ~BasisMatrix()
    { }
    
    /**
     * @returns \p m, the row-dimension of the matrix where the marix is
     * \f$ M \times N \f$.
     */
    virtual numeric_index_type m () const
    {
        libmesh_assert(modes.size() > 0);
        return modes[0]->size();
    }
    
    /**
     * @returns \p n, the column-dimension of the matrix where the marix
     * is \f$ M \times N \f$.
     */
    virtual numeric_index_type n () const
    {
        libmesh_assert(modes.size() > 0);
        return modes.size();
    }

    
    /**
     * Multiplies the matrix with \p arg and stores the result in \p
     * dest.
     */
    template <typename VecType>
    void vector_mult (NumericVector<T>& dest,
                      const VecType& arg) const
    {
        libmesh_assert(modes.size() > 0);
        libmesh_assert_equal_to(m(), dest.size());
        libmesh_assert_equal_to(n(), arg.size());
                                                        
        dest.zero();
        for (unsigned int i=0; i<n(); i++)
            dest.add(arg(i), *(modes[i]));
    }
    
    /**
     * Multiplies the matrix with \p arg and stores the result in \p
     * dest.
     */
    virtual void vector_mult (NumericVector<T>& dest,
                      const NumericVector<T>& arg) const
    {
        // not defined for multiplcation with NumericVector
        libmesh_assert(false);
    }
    

    /**
     * Multiplies the matrix with \p arg and stores the result in \p
     * dest.
     */
    template <typename VecType>
    void vector_mult_transpose (VecType& dest,
                                const NumericVector<T>& arg) const
    {
        libmesh_assert(modes.size() > 0);
        libmesh_assert_equal_to(m(), arg.size());
        libmesh_assert_equal_to(n(), dest.size());
        
        dest.setZero(n());
        for (unsigned int i=0; i<n(); i++)
            dest(i) = arg.dot(*(modes[i]));
    }

    
    /**
     * Multiplies the transpose of matrix with \p arg and stores the
     * result in \p dest.
     */
    virtual void vector_mult_transpose (NumericVector<T>& dest,
                                        const NumericVector<T>& arg) const
    {
        // not defined for multiplcation with NumericVector
        libmesh_assert(false);
    }

    /**
     * Multiplies the matrix with \p arg and adds the result to \p dest.
     */
    virtual void vector_mult_add (NumericVector<T>& dest,
                                  const NumericVector<T>& arg) const
    {
        // not defined for multiplcation with NumericVector
        libmesh_assert(false);
    }

    
    /**
     * Copies the diagonal part of the matrix into \p dest.
     */
    virtual void get_diagonal (NumericVector<T>& dest) const
    {
        // not defined for multiplcation with NumericVector
        libmesh_assert(false);
    }

    
    /*!
     *  Returns the vector that defines the \p i^th basis vector
     */
    
    virtual NumericVector<T>& basis(unsigned int i)
    {
        // not defined for multiplcation with NumericVector
        libmesh_assert(modes.size() > 0);
        libmesh_assert_less(i, modes.size());
        
        return *(modes[i]);
    }
   
    /*!
     *   vector of modes
     */
    std::vector<NumericVector<T>*> modes;
};


#endif
