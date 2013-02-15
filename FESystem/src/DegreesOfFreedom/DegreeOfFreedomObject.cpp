//
//  DegreeOfFreedomObject.cpp
//  FESystem
//
//  Created by Manav Bhatia on 4/9/12.
//  Copyright (c) 2012. All rights reserved.
//

// FESystem includes
#include "DegreesOfFreedom/DegreeOfFreedomObject.h"
#include "DegreesOfFreedom/DegreeOfFreedomUnit.h"
#include "Utils/AutoPtrTable.h"
#include "Base/FESystemExceptions.h"


FESystem::DegreeOfFreedom::DegreeOfFreedomObject::DegreeOfFreedomObject():
variable_table(NULL)
{
    
}
            

FESystem::DegreeOfFreedom::DegreeOfFreedomObject::~DegreeOfFreedomObject()
{
    if (this->variable_table != NULL)
        delete this->variable_table;
}


void
FESystem::DegreeOfFreedom::DegreeOfFreedomObject::init(FESystemUInt n_vars)
{
    FESystemAssert0(this->variable_table == NULL, FESystem::Exception::InvalidState);

    std::vector<FESystemUInt> el(1);
    el[0] = n_vars;
    this->variable_table = new FESystem::Utility::AutoPtrTable<FESystem::DegreeOfFreedom::DegreeOfFreedomUnit>;
    this->variable_table->reinit(el);
}



FESystemUInt
FESystem::DegreeOfFreedom::DegreeOfFreedomObject::getNDegreeOfFreedomUnits() const
{
    FESystemAssert0(this->variable_table != NULL, FESystem::Exception::InvalidState); // make sure that the table has been initialized

    return this->variable_table->getSizeVector()[0]; // this is a 1-D vector and thus the first dimension is the required data
}


FESystem::DegreeOfFreedom::DegreeOfFreedomUnit&
FESystem::DegreeOfFreedom::DegreeOfFreedomObject::addDegreeOfFreedomUnit(FESystemUInt var_id)
{
    std::vector<FESystemUInt> el(1);
    el[0] = var_id;

    FESystemAssert0(this->variable_table != NULL, FESystem::Exception::InvalidState); // make sure that the table has been initialized
    FESystemAssert0(!this->variable_table->ifValExists(el), FESystem::Exception::InvalidState); // make sure that the variable does not alredy exist

    FESystem::DegreeOfFreedom::DegreeOfFreedomUnit* dof = new FESystem::DegreeOfFreedom::DegreeOfFreedomUnit;
    
    this->variable_table->resetVal(el, dof); 
    return *dof;
}



FESystem::DegreeOfFreedom::DegreeOfFreedomUnit& 
FESystem::DegreeOfFreedom::DegreeOfFreedomObject::getDegreeOfFreedomUnit(FESystemUInt var_id)
{
    std::vector<FESystemUInt> el(1);
    el[0] = var_id;
    
    FESystemAssert0(this->variable_table != NULL, FESystem::Exception::InvalidState); // make sure that the table has been initialized
    FESystemAssert0(this->variable_table->ifValExists(el), FESystem::Exception::InvalidState); // make sure that the variable does not alredy exist
    
    return *(this->variable_table->getVal(el)); 
}


const FESystem::DegreeOfFreedom::DegreeOfFreedomUnit& 
FESystem::DegreeOfFreedom::DegreeOfFreedomObject::getDegreeOfFreedomUnit(FESystemUInt var_id) const
{
    std::vector<FESystemUInt> el(1);
    el[0] = var_id;
    
    FESystemAssert0(this->variable_table != NULL, FESystem::Exception::InvalidState); // make sure that the table has been initialized
    FESystemAssert0(this->variable_table->ifValExists(el), FESystem::Exception::InvalidState); // make sure that the variable does not alredy exist
    
    return *(this->variable_table->getVal(el)); 
}
            