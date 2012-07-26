/*
 *  PolynomialCoefficientList.cpp
 *  FESystem
 *
 *  Created by Manav  Bhatia on 1/10/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

// FESystem includes
#include "Surrogates/PolynomialCoefficientList.h"
#include "Surrogates/PolynomialCoefficientListIterator.h"
#include "Field/FieldBase.h"


FESystem::Surrogates::PolynomialCoefficientList::PolynomialCoefficientList():
field(NULL),
initialized(false)
{
  this->polynomial_iterator.reset(new FESystem::Surrogates::PolynomialCoefficientListIterator(*this));
}


FESystem::Surrogates::PolynomialCoefficientList::~PolynomialCoefficientList()
{
  this->polynomial_terms.clear();
}


void
FESystem::Surrogates::PolynomialCoefficientList::clear()
{
  this->polynomial_terms.clear();
  this->initialized = false;
  this->number_of_parameters = 0;
}


void
FESystem::Surrogates::PolynomialCoefficientList::prepareCoefficientList
(FESystem::Field::FieldBase<FESystemDouble>& f, FESystem::Surrogates::PolynomialType t)
{
  FESystemAssert0(!this->initialized,
                  PolynomialCoefficientListNotClearedBeforeReInitialization);
  
  this->field = &f;
  this->polynomial_type = t;  
  
  switch (t) 
  {
    case FESystem::Surrogates::ZEROTH_ORDER:
      this->addTerms(0);
      break;
      
    case FESystem::Surrogates::FIRST_ORDER:
    {
      this->addTerms(0);
      this->addTerms(1);
    }
      break;
      
    case FESystem::Surrogates::SECOND_ORDER_COMPLETE:
    {
      this->addTerms(0);
      this->addTerms(2);
    }
      break;
      
    case FESystem::Surrogates::SECOND_ORDER_NO_COUPLING:
    {
      this->addTerms(0);
      this->addTerms(2);
      this->removeCouplingTerms();
    }
      break;
      
    case FESystem::Surrogates::THIRD_ORDER_COMPLETE:
    {
      this->addTerms(0);
      this->addTerms(3);
    }
      break;
      
    case FESystem::Surrogates::THIRD_ORDER_NO_COUPLING:
    {
      this->addTerms(0);
      this->addTerms(3);
      this->removeCouplingTerms();
    }
      break;
      
    default:
      break;
  }
  
  this->initialized = true;
}


void
FESystem::Surrogates::PolynomialCoefficientList::addTerms(FESystemUInt order)
{
  FESystemUInt n_params = this->field->getNumberOfParameters();
  
  FESystem::Surrogates::PolynomialTerm p, p2;
  p.setNumberOfParameters(n_params);
  p.resetTermExponents();
  p2.setNumberOfParameters(n_params);
  p2.resetTermExponents();
  
  switch (order) 
  {
    case 0:
      this->addPolynomialTerm(p);
      break;
      
    default:
    {
      for (FESystemUInt k=0; k<order; k++)
      {
        std::vector<FESystem::Surrogates::PolynomialTerm> polynomial_terms_copy(this->polynomial_terms);
        
        for (FESystemUInt i=0; i<n_params; i++)
        {          
          p.resetTermExponents();
          p.polynomial_terms[i].second = 1;
          
          for (FESystemUInt j=0; j<polynomial_terms_copy.size(); j++)
          {
            polynomial_terms_copy[j].multiplyPolynomials(p, p2);
            this->addPolynomialTerm(p2);
          }
        }
      }
      break;
    }
  }
}




void
FESystem::Surrogates::PolynomialCoefficientList::removeCouplingTerms()
{
  // iterate over all the terms and if any term has exponent for more than one parameter, 
  // make it to not be considered
  
  std::vector<FESystem::Surrogates::PolynomialTerm>::iterator it, end;
  it = this->polynomial_terms.begin();
  end = this->polynomial_terms.end();
  
  for ( ; it != end; it++)
    if (it->ifCouplingTerm()) 
      it->if_keep = false;
}


void
FESystem::Surrogates::PolynomialCoefficientList::addPolynomialTerm
(const FESystem::Surrogates::PolynomialTerm& p)
{
  FESystemUInt n_params = this->field->getNumberOfParameters();

  std::pair<FESystemBoolean, FESystem::Surrogates::PolynomialTerm*> 
  term_pair = this->findPolynomialTerm(p);
  
  if (term_pair.first)
    term_pair.second->coefficient += p.coefficient;
  else 
    this->polynomial_terms.push_back(p);
}



std::pair<FESystemBoolean, FESystem::Surrogates::PolynomialTerm*>
FESystem::Surrogates::PolynomialCoefficientList::findPolynomialTerm
(const FESystem::Surrogates::PolynomialTerm& p)
{
  std::pair<FESystemBoolean, FESystem::Surrogates::PolynomialTerm*> r_val;
  std::vector<FESystem::Surrogates::PolynomialTerm>::iterator it, end;
  it = this->polynomial_terms.begin();
  end = this->polynomial_terms.end();
  
  for ( ; it != end; it++)
    if (it->ifSame(p)) break;
  
  if (it == end)
    r_val = std::pair<FESystemBoolean, FESystem::Surrogates::PolynomialTerm*>(false, NULL);
  else 
    r_val = std::pair<FESystemBoolean, FESystem::Surrogates::PolynomialTerm*>(true, &(*it));
  
  return r_val;
}



FESystem::Surrogates::PolynomialCoefficientListIterator& 
FESystem::Surrogates::PolynomialCoefficientList::getIterator()
{
  this->polynomial_iterator->begin();
  return *this->polynomial_iterator;
}



FESystemUInt 
FESystem::Surrogates::PolynomialCoefficientList::getTotalNumberOfCoefficients() const
{
  FESystemAssert0(this->initialized, 
                  FESystem::Surrogates::PolynomialCoefficientListNotInitialized);
  
  // only the number of active coefficients is returned
  return this->polynomial_terms.size();
}


FESystemUInt 
FESystem::Surrogates::PolynomialCoefficientList::getNumberOfActiveCoefficients() const
{
  FESystemAssert0(this->initialized, 
                  FESystem::Surrogates::PolynomialCoefficientListNotInitialized);

  // only the number of active coefficients is returned
  FESystemUInt n = 0;
  for (FESystemUInt i=0; i<this->polynomial_terms.size(); i++)
    if (this->polynomial_terms[i].if_keep)
      n += 1;

  return n;
}



const FESystem::Surrogates::PolynomialTerm&
FESystem::Surrogates::PolynomialCoefficientList::getPolynomialTerm(FESystemUInt n) const
{
  FESystemAssert2(n < this->getTotalNumberOfCoefficients(),
                  PolynomialTermNumberExceedsCurrentCoefficientList, 
                  this->getTotalNumberOfCoefficients(), n);
  
  return this->polynomial_terms[n];
}



FESystemBoolean
FESystem::Surrogates::PolynomialCoefficientList::ifKeepPolynomialTerm(FESystemUInt n) const
{
//  FESystemAssert0(this->initialized, 
//                  FESystem::Surrogates::PolynomialCoefficientListNotInitialized);
  FESystemAssert2(n < this->getTotalNumberOfCoefficients(), 
                  FESystem::Surrogates::PolynomialTermNumberExceedsCurrentCoefficientList,
                  this->getTotalNumberOfCoefficients(), n);
  
  return this->polynomial_terms[n].if_keep;
}


void
FESystem::Surrogates::PolynomialCoefficientList::write(std::ostream& out) const
{
  FESystemUInt total = this->polynomial_terms.size(),
  active = this->getNumberOfActiveCoefficients();
  
  out << "***********  Polynomial *********** " << std::endl
  << "Total number of terms : " << total << std::endl
  << "# Active terms: " << active << std::endl; 
  for (FESystemUInt i=0; i < total; i++)
  {
    if (this->ifKeepPolynomialTerm(i)) this->getPolynomialTerm(i).writeShort(out);
    if (i < total-1) out << "  +  " << std::endl; 
  }
  out << "# Inactive terms: " << (total-active) << std::endl; 
  for (FESystemUInt i=0; i < total; i++)
  {
    if (!(this->ifKeepPolynomialTerm(i))) 
    {
      this->getPolynomialTerm(i).writeShort(out);
      if (i < total-1) out << "  +  " << std::endl; 
    }
  }
  out << "*********************************** " << std::endl;
}


