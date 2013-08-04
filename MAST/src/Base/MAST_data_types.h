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

// Eigen includes
#include "Eigen/Dense"
using namespace Eigen;
using namespace libMesh;

typedef Matrix<Real, Dynamic, 1> RealVectorX;
typedef Matrix<Real, 3, 1> RealVector3;
typedef Matrix<Complex, Dynamic, 1> ComplexVectorX;
typedef Matrix<Complex, 3, 1> ComplexVector3;

typedef Matrix<Real, Dynamic, Dynamic> RealMatrixX;
typedef Matrix<Real, 3, 3> RealMatrix3;
typedef Matrix<Complex, Dynamic, Dynamic> ComplexMatrixX;
typedef Matrix<Complex, 3, 3> ComplexMatrix3;

#endif
