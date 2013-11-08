//
//  fem_operator_matrix.h
//  MAST
//
//  Created by Manav Bhatia on 6/17/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST__fem_operator_matrix_h__
#define __MAST__fem_operator_matrix_h__

// C++ includes
#include <vector>


// libMesh includes
#include "libmesh/dense_vector.h"
#include "libmesh/dense_matrix.h"


using namespace libMesh;


class FEMOperatorMatrix
{
public:
    FEMOperatorMatrix();
    
    virtual ~FEMOperatorMatrix();
    
    /*!
     *   clears the data structures
     */
    void clear();
    
    
    unsigned int m() const {return _n_interpolated_vars;}
    
    unsigned int n() const {return _n_discrete_vars*_n_dofs_per_var;}
    
    /*!
     *   this initializes the operator for number of rows and variables, assuming 
     *   that all variables has the same number of dofs. This is typically the case
     *   for structural strain operator matrices. Note that when this method is used
     *   the user must set the matrix entries by calling set_shape_functions
     */
    void reinit(unsigned int n_interpolated_vars,
                unsigned int n_discrete_vars,
                unsigned int n_discrete_dofs_per_var);

    /*!
     *   sets the shape function values for the block corresponding to 
     *   \par interpolated_var and \par discrete_var. This means that the row
     *   \par interpolated_var, the value in columns 
     *   \par discrete_vars*n_discrete_dofs_per_var - (discrete_vars+1)*n_discrete_dofs_per_var-1)
     *    will be set equal to \par shape_func .
     */
    void set_shape_function(unsigned int interpolated_var,
                            unsigned int discrete_var,
                            const DenseVector<Real>& shape_func);
    
    /*!
     *   this initializes all variables to use the same interpolation function. 
     *   It is assumed that the number of discrete vars is same as the number of
     *   interpolated vars. This is typically the case for fluid elements and 
     *   for structural element inertial matrix calculations
     */
    void reinit(unsigned int n_interpolated_vars, const DenseVector<Real>& shape_func);
    
    /*!
     *   res = [this] * v
     */
    template <typename T>
    void vector_mult(DenseVector<T>& res, const DenseVector<T>& v) const;
    

    /*!
     *   res = v^T * [this]
     */
    template <typename T>
    void vector_mult_transpose(DenseVector<T>& res, const DenseVector<T>& v) const;

    
    /*!
     *   [R] = [this] * [M]
     */
    template <typename T>
    void right_multiply(DenseMatrix<T>& r, const DenseMatrix<T>& m) const;
    
    
    /*!
     *   [R] = [this]^T * [M]
     */
    template <typename T>
    void right_multiply_transpose(DenseMatrix<T>& r, const DenseMatrix<T>& m) const;


    /*!
     *   [R] = [this]^T * [M]
     */
    template <typename T>
    void right_multiply_transpose(DenseMatrix<T>& r, const FEMOperatorMatrix& m) const;

    
    /*!
     *   [R] = [M] * [this]
     */
    template <typename T>
    void left_multiply(DenseMatrix<T>& r, const DenseMatrix<T>& m) const;
    
    
    /*!
     *   [R] = [M] * [this]^T
     */
    template <typename T>
    void left_multiply_transpose(DenseMatrix<T>& r, const DenseMatrix<T>& m) const;

    
protected:
    
    /*!
     *    number of rows of the operator
     */
    unsigned int _n_interpolated_vars;

    /*!
     *    number of discrete variables in the system
     */
    unsigned int _n_discrete_vars;

    /*!
     *    number of dofs for each variable
     */
    unsigned int _n_dofs_per_var;

    /*!
     *    stores the shape function values that defines the coupling 
     *    of i_th interpolated var and j_th discrete var. Stored in 
     *    column major format. NULL, if values are zero, otherwise the 
     *    value is set in the vector.
     */
    std::vector<DenseVector<Real>*>  _var_shape_functions;
};




inline
void
FEMOperatorMatrix::clear()
{
    _n_interpolated_vars = 0;
    _n_discrete_vars     = 0;
    _n_dofs_per_var      = 0;

    // iterate over the shape function entries and delete the non-NULL values
    std::vector<DenseVector<Real>*>::iterator it = _var_shape_functions.begin(),
    end = _var_shape_functions.end();
    
    for ( ; it!=end; it++)
        if ( *it != NULL)
            delete *it;
    
    _var_shape_functions.clear();
}




inline
void
FEMOperatorMatrix::reinit(unsigned int n_interpolated_vars,
                          unsigned int n_discrete_vars,
                          unsigned int n_discrete_dofs_per_var)
{
    this->clear();
    _n_interpolated_vars = n_interpolated_vars;
    _n_discrete_vars = n_discrete_vars;
    _n_dofs_per_var = n_discrete_dofs_per_var;
    _var_shape_functions.resize(_n_interpolated_vars*_n_discrete_vars);
    for (unsigned int i=0; i<_var_shape_functions.size(); i++)
        _var_shape_functions[i] = NULL;
}



inline
void
FEMOperatorMatrix::set_shape_function(unsigned int interpolated_var,
                                      unsigned int discrete_var,
                                      const DenseVector<Real>& shape_func)
{
    // make sure that reinit has been called.
    libmesh_assert(_var_shape_functions.size());
    
    // also make sure that the specified indices are within bounds
    libmesh_assert(interpolated_var < _n_interpolated_vars);
    libmesh_assert(discrete_var < _n_discrete_vars);
    
    DenseVector<Real>* vec = new DenseVector<Real>;
    *vec = shape_func;
    _var_shape_functions[discrete_var*_n_interpolated_vars+interpolated_var] = vec;
}



inline
void FEMOperatorMatrix::reinit(unsigned int n_vars, const DenseVector<Real>& shape_func)
{
    this->clear();
    
    _n_interpolated_vars = n_vars;
    _n_discrete_vars = n_vars;
    _n_dofs_per_var = shape_func.size();
    _var_shape_functions.resize(n_vars*n_vars);
    
    for (unsigned int i=0; i<n_vars; i++)
    {
        DenseVector<Real>*  vec = new DenseVector<Real>;
        *vec = shape_func;
        _var_shape_functions[i*n_vars+i] = vec;
    }
}



template <typename T>
inline
void FEMOperatorMatrix::vector_mult(DenseVector<T>& res, const DenseVector<T>& v) const
{
    libmesh_assert_equal_to(res.size(), _n_interpolated_vars);
    libmesh_assert_equal_to(v.size(), n());

    res.zero();
    unsigned int index = 0;

    for (unsigned int i=0; i<_n_interpolated_vars; i++) // row
        for (unsigned int j=0; j<_n_discrete_vars; j++) { // column
            index = j*_n_interpolated_vars+i;
            if (_var_shape_functions[index]) // check if this is non-NULL
                for (unsigned int k=0; k<_n_dofs_per_var; k++)
                    res(i) +=
                    (*_var_shape_functions[index])(k) * v(j*_n_dofs_per_var+k);
        }
}


template <typename T>
inline
void FEMOperatorMatrix::vector_mult_transpose(DenseVector<T>& res, const DenseVector<T>& v) const
{
    libmesh_assert_equal_to(res.size(), n());
    libmesh_assert_equal_to(v.size(), _n_interpolated_vars);
    
    res.zero();
    unsigned int index = 0;

    for (unsigned int i=0; i<_n_interpolated_vars; i++) // row
        for (unsigned int j=0; j<_n_discrete_vars; j++) { // column
            index = j*_n_interpolated_vars+i;
            if (_var_shape_functions[index]) // check if this is non-NULL
                for (unsigned int k=0; k<_n_dofs_per_var; k++)
                    res(j*_n_dofs_per_var+k) +=
                    (*_var_shape_functions[index])(k) * v(i);
        }
}


template <typename T>
inline
void FEMOperatorMatrix::right_multiply(DenseMatrix<T>& r, const DenseMatrix<T>& m) const
{
    libmesh_assert_equal_to(r.m(), _n_interpolated_vars);
    libmesh_assert_equal_to(r.n(), m.n());
    libmesh_assert_equal_to(m.m(), n());
    
    r.zero();
    unsigned int index = 0;

    for (unsigned int i=0; i<_n_interpolated_vars; i++) // row
        for (unsigned int j=0; j<_n_discrete_vars; j++) { // column of operator
            index = j*_n_interpolated_vars+i;
            if (_var_shape_functions[index]) { // check if this is non-NULL
                for (unsigned int l=0; l<m.n(); l++) // column of matrix
                    for (unsigned int k=0; k<_n_dofs_per_var; k++)
                        r(i,l) +=
                        (*_var_shape_functions[index])(k) * m(j*_n_dofs_per_var+k,l);
            }
        }
}


template <typename T>
inline
void FEMOperatorMatrix::right_multiply_transpose(DenseMatrix<T>& r, const DenseMatrix<T>& m) const
{
    libmesh_assert_equal_to(r.m(), n());
    libmesh_assert_equal_to(r.n(), m.n());
    libmesh_assert_equal_to(m.m(), _n_interpolated_vars);
    
    r.zero();
    unsigned int index = 0;

    for (unsigned int i=0; i<_n_interpolated_vars; i++) // row
        for (unsigned int j=0; j<_n_discrete_vars; j++) { // column of operator
            index = j*_n_interpolated_vars+i;
            if (_var_shape_functions[index]) { // check if this is non-NULL
                for (unsigned int l=0; l<m.n(); l++) // column of matrix
                    for (unsigned int k=0; k<_n_dofs_per_var; k++)
                        r(j*_n_dofs_per_var+k,l) +=
                        (*_var_shape_functions[index])(k) * m(i,l);
            }
        }
}



template <typename T>
inline
void FEMOperatorMatrix::right_multiply_transpose(DenseMatrix<T>& r, const FEMOperatorMatrix& m) const
{
    libmesh_assert_equal_to(r.m(), n());
    libmesh_assert_equal_to(r.n(), m.n());
    libmesh_assert_equal_to(_n_interpolated_vars, m._n_interpolated_vars);
    
    r.zero();
    unsigned int index_i, index_j = 0;

    for (unsigned int i=0; i<_n_discrete_vars; i++) // row of result
        for (unsigned int j=0; j<m._n_discrete_vars; j++) // column of result
            for (unsigned int k=0; k<_n_interpolated_vars; k++) {
                index_i = i*_n_interpolated_vars+k;
                index_j = j*m._n_interpolated_vars+k;
                if (_var_shape_functions[index_i] &&
                    m._var_shape_functions[index_j]) { // if shape function exists for both
                    const DenseVector<Real> &n1 = *_var_shape_functions[index_i],
                    &n2 = *m._var_shape_functions[index_j];
                    for (unsigned int i_n1=0; i_n1<n1.size(); i_n1++)
                        for (unsigned int i_n2=0; i_n2<n2.size(); i_n2++)
                            r (i*_n_dofs_per_var+i_n1,
                              j*m._n_dofs_per_var+i_n2) += n1(i_n1) * n2(i_n2);
                }
            }
}




template <typename T>
inline
void FEMOperatorMatrix::left_multiply(DenseMatrix<T>& r, const DenseMatrix<T>& m) const
{
    libmesh_assert_equal_to(r.m(), m.m());
    libmesh_assert_equal_to(r.n(), n());
    libmesh_assert_equal_to(m.n(), _n_interpolated_vars);
    
    r.zero();
    unsigned int index = 0;
    
    for (unsigned int i=0; i<_n_interpolated_vars; i++) // row
        for (unsigned int j=0; j<_n_discrete_vars; j++) { // column of operator
            index = j*_n_interpolated_vars+i;
            if (_var_shape_functions[index]) { // check if this is non-NULL
                for (unsigned int l=0; l<m.m(); l++) // rows of matrix
                    for (unsigned int k=0; k<_n_dofs_per_var; k++)
                        r(l,j*_n_dofs_per_var+k) +=
                        (*_var_shape_functions[index])(k) * m(l,i);
            }
        }
}


template <typename T>
inline
void FEMOperatorMatrix::left_multiply_transpose(DenseMatrix<T>& r, const DenseMatrix<T>& m) const
{
    libmesh_assert_equal_to(r.m(), m.m());
    libmesh_assert_equal_to(r.n(), _n_interpolated_vars);
    libmesh_assert_equal_to(m.n(), n());
    
    r.zero();
    unsigned int index = 0;
    
    for (unsigned int i=0; i<_n_interpolated_vars; i++) // row
        for (unsigned int j=0; j<_n_discrete_vars; j++) { // column of operator
            index = j*_n_interpolated_vars+i;
            if (_var_shape_functions[index]) { // check if this is non-NULL
                for (unsigned int l=0; l<m.m(); l++) // column of matrix
                    for (unsigned int k=0; k<_n_dofs_per_var; k++)
                        r(l,i) +=
                        (*_var_shape_functions[index])(k) * m(l,j*_n_dofs_per_var+k);
            }
        }
}




#endif // __MAST__fem_operator_matrix_h__
