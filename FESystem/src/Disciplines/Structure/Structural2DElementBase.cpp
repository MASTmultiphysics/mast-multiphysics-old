//
//  Structural2DElementBase.cpp
//  FESystem
//
//  Created by Manav Bhatia on 8/15/12.
//
//

// FESystem includes
#include "Disciplines/Structure/Structural2DElementBase.h"


FESystem::Structures::Structural2DElementBase::Structural2DElementBase():
FESystem::Structures::StructuralElementBase(),
th_val(0.0)
{
    
}

FESystem::Structures::Structural2DElementBase::~Structural2DElementBase()
{
    
}


void
FESystem::Structures::Structural2DElementBase::clear()
{
    FESystem::Structures::StructuralElementBase::clear();
    this->th_val = 0.0;
}


void
FESystem::Structures::Structural2DElementBase::initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe,
                                                          const FESystem::Quadrature::QuadratureBase& q_rule, FESystemDouble E, FESystemDouble nu, FESystemDouble rho, FESystemDouble th)
{
    FESystem::Structures::StructuralElementBase::initialize(elem, fe, q_rule, E, nu, rho);
    this->th_val = th;
}

