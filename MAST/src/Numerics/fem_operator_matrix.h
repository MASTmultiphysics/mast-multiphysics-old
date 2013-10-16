//
//  fem_operator_matrix.h
//  RealSolver
//
//  Created by Manav Bhatia on 6/17/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __RealSolver__fem_operator_matrix__
#define __RealSolver__fem_operator_matrix__

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
    
    ~FEMOperatorMatrix();
    
    
    unsigned int m() const {return _n_vars;}
    
    unsigned int n() const {return _total_dofs;}
    
    /*!
     *   this initializes all variables to use the same interpolation function
     */
    void reinit(unsigned int n_vars, const DenseVector<Real>& shape_func);
    
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
    
    unsigned int _n_vars;
    
    unsigned int _total_dofs;
    
    std::vector<unsigned int> _n_dofs_per_var;
    
    std::vector<DenseVector<Real> > _var_shape_functions;
    
};



inline
FEMOperatorMatrix::FEMOperatorMatrix():
_n_vars(0),
_total_dofs(0)
{
    
}


inline
FEMOperatorMatrix::~FEMOperatorMatrix()
{
    
}



inline
void FEMOperatorMatrix::reinit(unsigned int n_vars, const DenseVector<Real>& shape_func)
{
    _n_vars = n_vars;
    _n_dofs_per_var.resize(n_vars);
    _var_shape_functions.resize(n_vars);
    _total_dofs = 0;
    
    for (unsigned int i=0; i<n_vars; i++)
    {
        _var_shape_functions[i] = shape_func;
        _n_dofs_per_var[i] = shape_func.size();
        _total_dofs += shape_func.size();
    }
}



template <typename T>
inline
void FEMOperatorMatrix::vector_mult(DenseVector<T>& res, const DenseVector<T>& v) const
{
    libmesh_assert_equal_to(res.size(), _n_vars);
    libmesh_assert_equal_to(v.size(), _total_dofs);

    res.zero();

    unsigned int dof = 0;
    
    for (unsigned int i=0; i<_n_vars; i++) 
        for (unsigned int j=0; j<_n_dofs_per_var[i]; j++)
        {
            res(i) += _var_shape_functions[i](j) * v(dof);
            dof++;
        }
}


template <typename T>
inline
void FEMOperatorMatrix::vector_mult_transpose(DenseVector<T>& res, const DenseVector<T>& v) const
{
    libmesh_assert_equal_to(res.size(), _total_dofs);
    libmesh_assert_equal_to(v.size(), _n_vars);
    
    res.zero();
    
    unsigned int dof = 0;
    
    for (unsigned int i=0; i<_n_vars; i++)
        for (unsigned int j=0; j<_n_dofs_per_var[i]; j++)
        {
            res(dof) = _var_shape_functions[i](j) * v(i);
            dof++;
        }
}


template <typename T>
inline
void FEMOperatorMatrix::right_multiply(DenseMatrix<T>& r, const DenseMatrix<T>& m) const
{
    libmesh_assert_equal_to(r.m(), _n_vars);
    libmesh_assert_equal_to(r.n(), m.n());
    
    libmesh_assert_equal_to(m.m(), _total_dofs);
    
    r.zero();
    
    unsigned int off = 0;
    
    for (unsigned int i=0; i<_n_vars; i++)
    {
        for (unsigned int j=0; j<m.n(); j++)
            for (unsigned int k=0; k<_n_dofs_per_var[i]; k++)
                r(i, j) += _var_shape_functions[i](k) * m(off+k, j);
        off+= _n_dofs_per_var[i];
    }
}


template <typename T>
inline
void FEMOperatorMatrix::right_multiply_transpose(DenseMatrix<T>& r, const DenseMatrix<T>& m) const
{
    libmesh_assert_equal_to(r.m(), _total_dofs);
    libmesh_assert_equal_to(r.n(), m.n());
    
    libmesh_assert_equal_to(m.m(), _n_vars);
    
    r.zero();
    
    unsigned int dof = 0;
    
    for (unsigned int j=0; j<m.n(); j++)
    {
        dof = 0;
        for (unsigned int i=0; i<_n_vars; i++)
            for (unsigned int k=0; k<_n_dofs_per_var[i]; k++)
            {
                r(dof, j) = _var_shape_functions[i](k) * m(i, j);
                dof ++;
            }
    }
}



template <typename T>
inline
void FEMOperatorMatrix::right_multiply_transpose(DenseMatrix<T>& r, const FEMOperatorMatrix& m) const
{
    libmesh_assert_equal_to(r.m(), _total_dofs);
    libmesh_assert_equal_to(r.n(), m._total_dofs);
    
    libmesh_assert_equal_to(_n_vars, m._n_vars);
    
    r.zero();
    
    unsigned int row_off = 0, col_off=0;

    for (unsigned int i_var=0; i_var<_n_vars; i_var++)
    {
        for (unsigned int i=0; i<_n_dofs_per_var[i_var]; i++)
            for (unsigned int j=0; j<m._n_dofs_per_var[i_var]; j++)
                r(row_off+i, col_off+j) =
                (_var_shape_functions[i_var](i) * m._var_shape_functions[i_var](j));
                
        row_off += _n_dofs_per_var[i_var];
        col_off += m._n_dofs_per_var[i_var];
    }
}




template <typename T>
inline
void FEMOperatorMatrix::left_multiply(DenseMatrix<T>& r, const DenseMatrix<T>& m) const
{
    libmesh_assert_equal_to(r.m(), m.m());
    libmesh_assert_equal_to(r.n(), _total_dofs);
    
    libmesh_assert_equal_to(m.n(), _n_vars);
    
    r.zero();
    
    unsigned int dof = 0;
    
    for (unsigned int i=0; i<_n_vars; i++)
        for (unsigned int k=0; k<_n_dofs_per_var[i]; k++)
        {
            for (unsigned int j=0; j<m.m(); j++)
                r(j, dof) = _var_shape_functions[i](k) * m(j, i);
            dof++;
        }
}


template <typename T>
inline
void FEMOperatorMatrix::left_multiply_transpose(DenseMatrix<T>& r, const DenseMatrix<T>& m) const
{
    libmesh_assert_equal_to(r.m(), m.m());
    libmesh_assert_equal_to(r.n(), _n_vars);
    
    libmesh_assert_equal_to(m.n(), _total_dofs);
    
    r.zero();
    
    unsigned int off = 0;
    
    for (unsigned int i=0; i<_n_vars; i++) // column of the B (and r) matrix
    {
        for (unsigned int j=0; j<m.m(); j++) // row of the m (and r) matrix
            for (unsigned int k=0; k<_n_dofs_per_var[i]; k++)
                r(j, i) += _var_shape_functions[i](k) * m(j, off+k);
        off += _n_dofs_per_var[i];
    }
}




#endif /* defined(__RealSolver__fem_operator_matrix__) */
