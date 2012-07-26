/*
 *  PolynomialCoefficientListIterator.cpp
 *  FESystem
 *
 *  Created by Manav  Bhatia on 1/10/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "Surrogates/PolynomialCoefficientListIterator.h"
#include "Surrogates/PolynomialCoefficientList.h"


FESystem::Surrogates::PolynomialCoefficientListIterator::PolynomialCoefficientListIterator
(const FESystem::Surrogates::PolynomialCoefficientList& c_list):
coeff_list(c_list),
active_term_counter(0),
counter(0)
{
  
}


FESystem::Surrogates::PolynomialCoefficientListIterator::~PolynomialCoefficientListIterator()
{
  
}


void 
FESystem::Surrogates::PolynomialCoefficientListIterator::increment()
{
  this->counter++;
  if (this->ifEnd()) return;
  while (!(this->coeff_list.ifKeepPolynomialTerm(this->counter)))
  {
    if (this->ifEnd()) 
      return;
    else 
      this->counter++;
  }

  this->active_term_counter++;
}


void
FESystem::Surrogates::PolynomialCoefficientListIterator::decrement()
{
  this->counter--;
  if (this->ifReverseEnd()) return;
  while (!(this->coeff_list.ifKeepPolynomialTerm(this->counter)))
  {
    if (this->ifReverseEnd()) 
      return;
    else 
      this->counter--;
  }
  
  this->active_term_counter--;
}



void 
FESystem::Surrogates::PolynomialCoefficientListIterator::begin()
{
  this->counter = -1;
  this->active_term_counter = -1;
  this->increment();
}



void
FESystem::Surrogates::PolynomialCoefficientListIterator::reverseBegin()
{
  this->counter = this->coeff_list.getTotalNumberOfCoefficients();
  this->active_term_counter = this->coeff_list.getTotalNumberOfCoefficients();
  this->decrement();
}



FESystemBoolean
FESystem::Surrogates::PolynomialCoefficientListIterator::ifBegin() const
{
  return (this->counter == 0);
}


FESystemBoolean
FESystem::Surrogates::PolynomialCoefficientListIterator::ifReverseEnd() const
{
  return (this->counter < 0);
}


FESystemBoolean
FESystem::Surrogates::PolynomialCoefficientListIterator::ifEnd() const
{
  return (this->counter > (this->coeff_list.getTotalNumberOfCoefficients()-1));
}



FESystemInt 
FESystem::Surrogates::PolynomialCoefficientListIterator::getID() const
{
  return this->active_term_counter;
}



template <typename ValType>
ValType 
FESystem::Surrogates::PolynomialCoefficientListIterator::getCurrentValue
(const FESystem::Numerics::VectorBase<ValType>& vec)
{
  const PolynomialTerm& p = this->coeff_list.getPolynomialTerm(this->getID());
  return p.calculateValue(vec);
}


/***************************************************************************************/
// Template instantiations for some generic classes

template FESystemDouble FESystem::Surrogates::PolynomialCoefficientListIterator::getCurrentValue
(const FESystem::Numerics::VectorBase<FESystemDouble>& vec);



/***************************************************************************************/
