//
//  structural_basis_matrix.cpp
//  MAST
//
//  Created by Manav Bhatia on 11/08/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//


// MAST includes
#include "Numerics/basis_matrix.h"

template <typename T>
BasisMatrix<T>::BasisMatrix(const libMesh::Parallel::Communicator &comm_in):
libMesh::ShellMatrix<T>(comm_in)
{ }

template <typename T>
BasisMatrix<T>::~BasisMatrix()
{ }


// explicit instantiations
template BasisMatrix<libMesh::Real>::BasisMatrix(const libMesh::Parallel::Communicator &comm_in);
template BasisMatrix<libMesh::Complex>::BasisMatrix(const libMesh::Parallel::Communicator &comm_in);
