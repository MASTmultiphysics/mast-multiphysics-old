//
//  StructuralElementBase.cpp
//  FESystem
//
//  Created by Manav Bhatia on 8/12/12.
//
//


// FESystem includes
#include "Disciplines/Structure/StructuralElementBase.h"
#include "Base/FESystemExceptions.h"


FESystem::Structures::StructuralElementBase::StructuralElementBase():
if_initialized(false),
geometric_elem(NULL),
quadrature(NULL),
finite_element(NULL),
E_val(0.0),
nu_val(0.0),
G_val(0.0),
rho_val(0.0)
{
    
}


FESystem::Structures::StructuralElementBase::~StructuralElementBase()
{
    
}


std::string
FESystem::Structures::StructuralElementBase::getVariableName(FESystem::Structures::StructuralVariable var)
{
    std::string name;
    
    switch (var)
    {
        case FESystem::Structures::U_DISP:
            name = "U_DISP";
            break;
            
        case FESystem::Structures::V_DISP:
            name = "V_DISP";
            break;

        case FESystem::Structures::W_DISP:
            name = "W_DISP";
            break;

        case FESystem::Structures::THETA_X:
            name = "THETA_X";
            break;

        case FESystem::Structures::THETA_Y:
            name = "THETA_Y";
            break;

        case FESystem::Structures::THETA_Z:
            name = "THETA_Z";
            break;
        
        default:
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
    }
    
    return name;
}


FESystem::Structures::StructuralVariable
FESystem::Structures::StructuralElementBase::getVariableEnum(const std::string &var)
{
    FESystem::Structures::StructuralVariable name;
    
    if (var == "U_DISP")
        name = FESystem::Structures::U_DISP;
    else if (var == "V_DISP")
        name = FESystem::Structures::V_DISP;
    else if (var == "W_DISP")
        name = FESystem::Structures::W_DISP;
    else if (var == "THETA_X")
        name = FESystem::Structures::THETA_X;
    else if (var == "THETA_Y")
        name = FESystem::Structures::THETA_Y;
    else if (var == "THETA_Z")
        name = FESystem::Structures::THETA_Z;
    else
        FESystemAssert0(false, FESystem::Exception::InvalidValue);
    
    return name;
}


void
FESystem::Structures::StructuralElementBase::clear()
{
    this->geometric_elem = NULL;
    this->quadrature = NULL;
    this->finite_element = NULL;
    this->E_val  = 0.0;
    this->nu_val = 0.0;
    this->G_val = 0.0;
    this->rho_val = 0.0;

    this->if_initialized = false;
}



void
FESystem::Structures::StructuralElementBase::initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe,
                                                        const FESystem::Quadrature::QuadratureBase& q_rule, FESystemDouble E, FESystemDouble nu, FESystemDouble rho)
{
    FESystemAssert0(!this->if_initialized, FESystem::Exception::InvalidFunctionCall);
    
    this->if_initialized = true;
    this->geometric_elem = &elem;
    this->quadrature = &q_rule;
    this->finite_element = &fe;

    this->E_val  = E;
    this->nu_val = nu;
    this->G_val = E/2.0/(1.0+nu);
    this->rho_val = rho;
}


