/*
 *  PolynomialTerm.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 2/17/11.
 *  Copyright 2011 . All rights reserved.
 *
 */

#ifndef __fesystem_polynomial_term_h__
#define __fesystem_polynomial_term_h__

// C++ includes
#include <vector>


// FESystem includes
#include "Base/FESystemTypes.h"
#include "Base/FESystemExceptions.h"
#include "Surrogates/PolynomialType.h"


namespace FESystem
{
  // Forward declerations
  namespace Numerics {template <typename ValType> class VectorBase;}
  
  namespace Surrogates
  {
    
    class PolynomialTerm 
    {
    public:
      PolynomialTerm();
      
      PolynomialTerm(const PolynomialTerm& p);
      
      
      void multiplyPolynomials(FESystem::Surrogates::PolynomialTerm& p1, 
                               FESystem::Surrogates::PolynomialTerm& p_res);
      
      void setNumberOfParameters(FESystemUInt n);
      
      FESystemUInt getNumberOfParameters() const;
      
      void resetTermExponents();
      
      FESystemBoolean ifSame(const FESystem::Surrogates::PolynomialTerm& p) const;
      
      FESystemBoolean ifCouplingTerm() const;
      
      FESystemDouble calculateValue(const FESystem::Numerics::VectorBase<FESystemDouble>& vec) const;
      
      FESystemDouble getCoefficient() const;      
      
      void writeLong(std::ostream& out) const;

      void writeShort(std::ostream& out) const;
      
      std::vector<std::pair<FESystemUInt, FESystemUInt> > polynomial_terms; 
      
      FESystemDouble coefficient;
      
      FESystemUInt n_params; 
      
      FESystemBoolean if_keep; 
    };

    
    DeclareException2(PolynomialMismatch, FESystemUInt, FESystemUInt,
                      << "Polynomial Terms Have Different Number Of Parameters\n"
                      << "Polynomial 1: " << Arg1 << "\n"
                      << "Polynomial 2: " << Arg2 << "\n");

    DeclareException2(ParameterSizeMismatch, FESystemUInt, FESystemUInt,
                      << "Number of parameters different current term.\n"
                      << "Current Number of Parameters: " << Arg1 << "\n"
                      << "Given Number: " << Arg2 << "\n");
    
  }
}




#endif // __fesystem_polynomial_term_h__

