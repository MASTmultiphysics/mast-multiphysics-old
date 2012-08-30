//
//  Structural1DElementBase.cpp
//  FESystem
//
//  Created by Manav Bhatia on 8/15/12.
//
//


// FESystem includes
#include "Disciplines/Structure/Structural1DElementBase.h"
#include "Base/FESystemExceptions.h"


FESystem::Structures::Structural1DElementBase::Structural1DElementBase():
FESystem::Structures::StructuralElementBase(),
area_val(0.0)
{
    
}

FESystem::Structures::Structural1DElementBase::~Structural1DElementBase()
{
    
}


FESystemDouble
FESystem::Structures::Structural1DElementBase::getArea() const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    return this->area_val;
}


void
FESystem::Structures::Structural1DElementBase::clear()
{
    FESystem::Structures::StructuralElementBase::clear();
    this->area_val = 0.0;
}




void
FESystem::Structures::Structural1DElementBase::initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, const FESystem::Quadrature::QuadratureBase& q_rule,
                                                          FESystemDouble E, FESystemDouble nu, FESystemDouble rho, FESystemDouble area)
{
    FESystem::Structures::StructuralElementBase::initialize(elem, fe, q_rule, E, nu, rho);
    this->area_val = area;
}


