//
//  fem_operator_matrix.cpp
//  MAST
//
//  Created by Manav Bhatia on 6/17/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

// MAST includes
#include "Numerics/fem_operator_matrix.h"


FEMOperatorMatrix::FEMOperatorMatrix():
_n_interpolated_vars(0),
_n_discrete_vars(0),
_n_dofs_per_var(0)
{
    
}


FEMOperatorMatrix::~FEMOperatorMatrix()
{
    this->clear();
}


