//
//  LinearBeamElementBase.cpp
//  FESystem
//
//  Created by Manav Bhatia on 8/23/12.
//
//


// FESystem includes
#include "Disciplines/Structure/LinearBeamElementBase.h"
#include "Mesh/ElemBase.h"


FESystem::Structures::LinearBeamElementBase::LinearBeamElementBase():
FESystem::Structures::Structural1DElementBase(),
if_include_chordwise_stiffness(false),
I_tr_val(0.0),
I_ch_val(0.0)
{
    
}
            


FESystem::Structures::LinearBeamElementBase::~LinearBeamElementBase()
{
    
}


FESystemUInt
FESystem::Structures::LinearBeamElementBase::getNElemDofs() const
{
    if (this->if_include_chordwise_stiffness)
        return 4*this->geometric_elem->getNNodes();
    else
        return 2*this->geometric_elem->getNNodes();
}



FESystemBoolean
FESystem::Structures::LinearBeamElementBase::ifIncludeChordwiseStiffness()
{
    return this->if_include_chordwise_stiffness;
}
            



void
FESystem::Structures::LinearBeamElementBase::getActiveElementMatrixIndices(std::vector<FESystemUInt>& vec)
{
    FESystemUInt n = this->geometric_elem->getNNodes();
    
    if (this->if_include_chordwise_stiffness)
    {
        vec.resize(4*n);
        for (FESystemUInt i=0; i<4*n; i++) vec[i] = n + i; // v-, w-displacement, theta-x and theta-y
    }
    else
    {
        vec.resize(2*n);
        for (FESystemUInt i=0; i<n; i++) vec[i] = 2*n + i; // w-displacement
        for (FESystemUInt i=0; i<n; i++) vec[i+n] = 4*n + i; // theta-y
    }
    
}
            

void
FESystem::Structures::LinearBeamElementBase::clear()
{
    this->if_include_chordwise_stiffness = false;
    this->I_tr_val = 0.0;
    this->I_ch_val = 0.0;
}


void
FESystem::Structures::LinearBeamElementBase::initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, const FESystem::Quadrature::QuadratureBase& q_rule,
                                                        FESystemDouble E, FESystemDouble nu, FESystemDouble rho, FESystemDouble I_tr, FESystemDouble I_ch, FESystemBoolean if_chordwise)
{
    FESystem::Structures::StructuralElementBase::initialize(elem,fe,q_rule,E,nu,rho);
    this->if_include_chordwise_stiffness = if_include_chordwise_stiffness;
    this->I_tr_val = I_tr;
    this->I_ch_val = I_ch;
}
            
