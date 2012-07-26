//
//  StructuralDiscretization.cpp
//  FESystem
//
//  Created by Manav Bhatia on 4/27/12.
//  Copyright (c) 2012. All rights reserved.
//

// FESystem includes
#include "Disciplines/StructuralDiscretization.h"


FESystem::Discretization::StructuralDiscretization::StructuralDiscretization(FESystem::Mesh::MeshBase& m, FESystem::Base::DegreeOfFreedomMap& dofs):
FESystem::Discretization::DiscretizationBase(m, dofs)
{
    
}

FESystem::Discretization::StructuralDiscretization::~StructuralDiscretization()
{
    
}


void
FESystem::Discretization::StructuralDiscretization::evaluateJacobian(FESystemUInt iter_num, 
                                                                     const FESystem::Numerics::VectorBase<FESystemDouble>& vec, 
                                                                     FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    
}



void
FESystem::Discretization::StructuralDiscretization::evaluateTransientSolverData(FESystemUInt transient_iter_num, FESystemUInt sub_iter_num, FESystemDouble time,
                                                                                const FESystem::Numerics::VectorBase<FESystemDouble>& x_vec,
                                                                                FESystem::Numerics::VectorBase<FESystemDouble>& x_dot_vec,
                                                                                FESystem::Numerics::MatrixBase<FESystemDouble>& x_dot_jac_mat)
{
    
}

