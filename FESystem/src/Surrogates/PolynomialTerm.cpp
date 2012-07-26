/*
 *  PolynomialTerm.cpp
 *  FESystem
 *
 *  Created by Manav Bhatia on 2/17/11.
 *  Copyright 2011 . All rights reserved.
 *
 */

// FESystem includes
#include "Surrogates/PolynomialTerm.h"
#include "Numerics/VectorBase.h"


FESystem::Surrogates::PolynomialTerm::PolynomialTerm():
n_params(0),
coefficient(0.0),
if_keep(true)
{
  
}

FESystem::Surrogates::PolynomialTerm::PolynomialTerm(const PolynomialTerm& p):
n_params(p.n_params),
coefficient(p.coefficient),
if_keep(p.if_keep)
{
  this->polynomial_terms.resize(this->n_params);
  
  for( FESystemUInt i=0; i < this->n_params; i++)
    this->polynomial_terms[i] = 
    std::pair<FESystemUInt, FESystemUInt>(p.polynomial_terms[i].first,
                                          p.polynomial_terms[i].second);
}


void 
FESystem::Surrogates::PolynomialTerm::multiplyPolynomials(FESystem::Surrogates::PolynomialTerm& p1, 
                                                          FESystem::Surrogates::PolynomialTerm& p_res)
{
  FESystemAssert2((this->n_params == p1.n_params), 
                  FESystem::Surrogates::PolynomialMismatch,
                  this->n_params, p1.n_params); 
  for (FESystemUInt i=0; i < this->n_params; i++)
  {
    p_res.polynomial_terms[i].first = this->polynomial_terms[i].first;
    p_res.polynomial_terms[i].second = this->polynomial_terms[i].second + p1.polynomial_terms[i].second;
  }
  p_res.coefficient = this->coefficient * p_res.coefficient;
}

void 
FESystem::Surrogates::PolynomialTerm::setNumberOfParameters(FESystemUInt n)
{
  this->n_params = n;
  this->polynomial_terms.resize(n);
  this->resetTermExponents();
}


FESystemUInt
FESystem::Surrogates::PolynomialTerm::getNumberOfParameters() const
{
  return this->n_params;
}


void
FESystem::Surrogates::PolynomialTerm::resetTermExponents()
{
  this->coefficient = 1.0;
  for( FESystemUInt i=0; i < this->n_params; i++)
    this->polynomial_terms[i] = std::pair<FESystemUInt, FESystemUInt>(i,0);
}



FESystemBoolean
FESystem::Surrogates::PolynomialTerm::ifCouplingTerm() const
{
  std::vector<std::pair<FESystemUInt, FESystemUInt> >::const_iterator it, end;
  it = this->polynomial_terms.begin();
  end = this->polynomial_terms.end();
  
  FESystemUInt n_terms = 0;
  for ( ; it != end; it++)
    if (it->second > 0)
      n_terms += 1;
  
  if (n_terms > 1) return true;
  
  return false;
}



FESystemDouble 
FESystem::Surrogates::PolynomialTerm::calculateValue
(const FESystem::Numerics::VectorBase<FESystemDouble>& vec) const
{
  FESystemAssert2(vec.getSize() == this->getNumberOfParameters(),
                  FESystem::Surrogates::ParameterSizeMismatch,
                  this->getNumberOfParameters(), vec.getSize());
  
  FESystemDouble val = 1.0;
  
  for (FESystemUInt i=0; i < vec.getSize(); i++)
    if (this->polynomial_terms[i].second != 0)
      val *= pow(vec.getVal(i), this->polynomial_terms[i].second);
      
  val *= this->coefficient;

  return val;
}


FESystemBoolean
FESystem::Surrogates::PolynomialTerm::ifSame(const PolynomialTerm& p) const
{
  FESystemBoolean same = true;
  if (p.n_params != this->n_params)
    same = false;
  else 
  {
    for (FESystemUInt i=0; i<this->n_params; i++)
      if ((p.polynomial_terms[i].first != this->polynomial_terms[i].first) || 
          (p.polynomial_terms[i].second != this->polynomial_terms[i].second))
      {
        same = false;
        break;
      }
  }
  
  return same;
}


void
FESystem::Surrogates::PolynomialTerm::writeLong(std::ostream& out) const
{
  out << "Polynomial Term: " << std::endl
  << "if keep: " << this->if_keep << std::endl 
  << "# coefficient: " << this->coefficient << std::endl
  << "# parameters: " << this->n_params << std::endl;
  
  FESystemUInt n=0;
  for( FESystemUInt i=0; i < this->n_params; i++)
  {
    out << "x" << this->polynomial_terms[i].first 
    << "^" <<  this->polynomial_terms[i].second;
    if (i < this->n_params-1) out << "   ";
    n++;
    if (n == 5)
    {
      n = 0;
      out << std::endl;
    }
  }
  if (n > 0) out << std::endl << std::endl;
}


void
FESystem::Surrogates::PolynomialTerm::writeShort(std::ostream& out) const
{
  FESystemUInt n=0;
  out << this->coefficient << " * ";
  for( FESystemUInt i=0; i < this->n_params; i++)
  {
    out << "x" << this->polynomial_terms[i].first 
    << "^" <<  this->polynomial_terms[i].second;
    if (i < this->n_params-1) out << "  *  ";
    n++;
    if (n == 5)
    {
      n = 0;
      out << std::endl;
    }
  }
  if (n > 0) out << std::endl;
}


