//
//  MAST_data_types.h
//  MAST
//
//  Created by Manav Bhatia on 7/26/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef MAST_data_types_h
#define MAST_data_types_h

// libMesh includes
#include "libmesh/libmesh_common.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"

// Eigen includes
#include "Eigen/Dense"
using namespace Eigen;

typedef libMesh::Real Real;
typedef libMesh::Complex Complex;

typedef Matrix<Real, Dynamic, 1> RealVectorX;
typedef Matrix<Real, 3, 1> RealVector3;
typedef Matrix<Complex, Dynamic, 1> ComplexVectorX;
typedef Matrix<Complex, 3, 1> ComplexVector3;

typedef Matrix<Real, Dynamic, Dynamic> RealMatrixX;
typedef Matrix<Real, 3, 3> RealMatrix3;
typedef Matrix<Complex, Dynamic, Dynamic> ComplexMatrixX;
typedef Matrix<Complex, 3, 3> ComplexMatrix3;

typedef libMesh::DenseMatrix<Real> DenseRealMatrix;
typedef libMesh::DenseMatrix<Complex> DenseComplexMatrix;

typedef libMesh::DenseVector<Real> DenseRealVector;
typedef libMesh::DenseVector<Complex> DenseComplexVector;

#endif
