//
//  structural_basis_matrix.cpp
//  MAST
//
//  Created by Manav Bhatia on 11/08/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//


// MAST includes
#include "Base/MAST_data_types.h"
#include "Numerics/basis_matrix.h"

template <typename T>
BasisMatrix<T>::BasisMatrix(const libMesh::Parallel::Communicator &comm_in):
libMesh::ShellMatrix<T>(comm_in)
{ }

template <typename T>
BasisMatrix<T>::~BasisMatrix()
{ }


// explicit instantiations
template BasisMatrix<Real>::BasisMatrix(const libMesh::Parallel::Communicator &comm_in);
template BasisMatrix<Complex>::BasisMatrix(const libMesh::Parallel::Communicator &comm_in);
