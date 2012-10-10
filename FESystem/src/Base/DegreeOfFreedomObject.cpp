//
//  DegreeOfFreedomObject.cpp
//  FESystem
//
//  Created by Manav Bhatia on 4/9/12.
//  Copyright (c) 2012. All rights reserved.
//

// FESystem includes
#include "Base/DegreeOfFreedomObject.h"
#include "Base/DegreeOfFreedomUnit.h"
#include "Utils/AutoPtrTable.h"
#include "Base/FESystemExceptions.h"


FESystem::Base::DegreeOfFreedomObject::DegreeOfFreedomObject():
variable_table(NULL)
{
    
}
            

FESystem::Base::DegreeOfFreedomObject::~DegreeOfFreedomObject()
{

}


void
FESystem::Base::DegreeOfFreedomObject::init(FESystemUInt n_vars)
{
    FESystemAssert0(this->variable_table.get() == NULL, FESystem::Exception::InvalidState);

    static std::vector<FESystemUInt> el(1);
    el[0] = n_vars;
    this->variable_table.reset(new FESystem::Utility::AutoPtrTable<FESystem::Base::DegreeOfFreedomUnit>);
    this->variable_table->reinit(el);
}



FESystemUInt
FESystem::Base::DegreeOfFreedomObject::getNDegreeOfFreedomUnits() const
{
    FESystemAssert0(this->variable_table.get() != NULL, FESystem::Exception::InvalidState); // make sure that the table has been initialized

    return this->variable_table->getSizeVector()[0]; // this is a 1-D vector and thus the first dimension is the required data
}


FESystem::Base::DegreeOfFreedomUnit&
FESystem::Base::DegreeOfFreedomObject::addDegreeOfFreedomUnit(FESystemUInt var_id)
{
    static std::vector<FESystemUInt> el(1);
    el[0] = var_id;

    FESystemAssert0(this->variable_table.get() != NULL, FESystem::Exception::InvalidState); // make sure that the table has been initialized
    FESystemAssert0(!this->variable_table->ifValExists(el), FESystem::Exception::InvalidState); // make sure that the variable does not alredy exist

    FESystem::Base::DegreeOfFreedomUnit* dof = new FESystem::Base::DegreeOfFreedomUnit;
    
    this->variable_table->resetVal(el, dof); 
    return *dof;
}



FESystem::Base::DegreeOfFreedomUnit& 
FESystem::Base::DegreeOfFreedomObject::getDegreeOfFreedomUnit(FESystemUInt var_id)
{
    static std::vector<FESystemUInt> el(1);
    el[0] = var_id;
    
    FESystemAssert0(this->variable_table.get() != NULL, FESystem::Exception::InvalidState); // make sure that the table has been initialized
    FESystemAssert0(this->variable_table->ifValExists(el), FESystem::Exception::InvalidState); // make sure that the variable does not alredy exist
    
    return *(this->variable_table->getVal(el)); 
}


const FESystem::Base::DegreeOfFreedomUnit& 
FESystem::Base::DegreeOfFreedomObject::getDegreeOfFreedomUnit(FESystemUInt var_id) const
{
    static std::vector<FESystemUInt> el(1);
    el[0] = var_id;
    
    FESystemAssert0(this->variable_table.get() != NULL, FESystem::Exception::InvalidState); // make sure that the table has been initialized
    FESystemAssert0(this->variable_table->ifValExists(el), FESystem::Exception::InvalidState); // make sure that the variable does not alredy exist
    
    return *(this->variable_table->getVal(el)); 
}
            
