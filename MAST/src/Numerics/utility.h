//
//  utility.h
//  MAST
//
//  Created by Manav Bhatia on 6/10/14.
//  Copyright (c) 2014 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_utility_h__
#define __MAST_utility_h__


// MAST includes
#include "Base/MAST_data_types.h"

// libMesh includes
#include "libmesh/dense_vector.h"


namespace MAST {
    
    template <typename ValType>
    void transform_to_elem_vector(libMesh::DenseVector<ValType>& v,
                                  const DenseRealVector& v_real);
    
    
    template <>
    inline void transform_to_elem_vector(DenseRealVector& v,
                                         const DenseRealVector& v_real) {
        // make sure that the real vector is twice the size of the dense vector
        const unsigned int n = v.size();
        libmesh_assert_equal_to(v_real.size(), n);
        v = v_real;
    }
    
    
    template <>
    inline void transform_to_elem_vector(DenseComplexVector& v,
                                  const DenseRealVector& v_real) {
        // make sure that the real vector is twice the size of the dense vector
        const unsigned int n = v.size();
        libmesh_assert_equal_to(v_real.size(), 2*n);
        
        for (unsigned int i=0; i<n; i++)
            v(i) = libMesh::Complex(v_real(i), v_real(i+n));
    }

    
    /*!
     *    All calculations in MAST are done using Real numbers. The complex
     *    variables are divided into two unknowns, one each for the real and
     *    imaginary variables. This provides a template method to add a real
     *    or complex vector to the assembled vector.
     */
    template <typename ValType>
    void add_to_assembled_vector(DenseRealVector& assembled_vec,
                                 const libMesh::DenseVector<ValType>& elem_vec);
    
    
    /*!
     *    All calculations in MAST are done using Real numbers. The complex
     *    variables are divided into two unknowns, one each for the real and
     *    imaginary variables. This provides a template method to add a real
     *    or complex matrix to the assembled matrix.
     */
    template <typename ValType>
    inline void add_to_assembled_matrix(DenseRealMatrix& assembled_mat,
                                 const libMesh::DenseMatrix<ValType>& elem_mat);
    
    
    template <>
    inline void
    add_to_assembled_matrix(DenseRealMatrix& assembled_mat,
                            const DenseRealMatrix& elem_mat) {
        assembled_mat += elem_mat;
    }
    
    
    template <>
    inline void
    add_to_assembled_vector(DenseRealVector& assembled_vec,
                            const DenseRealVector& elem_vec) {
        assembled_vec += elem_vec;
    }
    
    
    template <>
    inline void
    add_to_assembled_matrix(DenseRealMatrix& assembled_mat,
                            const DenseComplexMatrix& elem_mat) {
        
        // make sure the the assembled mat is twice the size of the elem mat
        const unsigned int m = elem_mat.m(), n = elem_mat.n();
        libmesh_assert_equal_to(assembled_mat.m(), m*2);
        libmesh_assert_equal_to(assembled_mat.n(), n*2);
        for (unsigned int i=0; i<m; i++)
            for (unsigned int j=0; j<n; j++) {
                assembled_mat(i,j)     +=  std::real(elem_mat(i,j));
                assembled_mat(i+m,j+n) +=  std::real(elem_mat(i,j));
                assembled_mat(i,j+n)   += -std::imag(elem_mat(i,j));
                assembled_mat(i+m,j)   +=  std::imag(elem_mat(i,j));
            }
    }
    
    
    template <>
    inline void
    add_to_assembled_vector(DenseRealVector& assembled_vec,
                            const DenseComplexVector& elem_vec) {
        // make sure the the assembled mat is twice the size of the elem mat
        const unsigned int n = elem_vec.size();
        libmesh_assert_equal_to(assembled_vec.size(), n*2);
        
        for (unsigned int i=0; i<n; i++) {
            assembled_vec(i)     +=   std::real(elem_vec(i));
            assembled_vec(i+n)   +=  std::imag(elem_vec(i));
        }
    }
    
}




#endif // __MAST_utility_h__
