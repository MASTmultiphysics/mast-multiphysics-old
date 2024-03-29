/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2014  Manav Bhatia
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

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
