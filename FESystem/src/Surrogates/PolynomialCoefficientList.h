/*
 *  PolynomialCoefficientList.h
 *  FESystem
 *
 *  Created by Manav  Bhatia on 1/10/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __fesystem_polynomial_coefficient_list_h__
#define __fesystem_polynomial_coefficient_list_h__

// C++ includes
#include <vector>
#include <memory>

// FESystem includes
#include "Base/FESystemTypes.h"
#include "Base/FESystemExceptions.h"
#include "Surrogates/PolynomialTerm.h"
#include "Surrogates/PolynomialType.h"

// Forward declerations
namespace FESystem
{
  namespace Field
  {
    template <typename ValType> class FieldBase;
  }
}

namespace FESystem
{
  namespace Surrogates
  {
    // Forward Declerations
    class PolynomialCoefficientListIterator;
    
    class PolynomialCoefficientList
    {
    public:
      PolynomialCoefficientList();
      
      ~PolynomialCoefficientList();
      
      void clear();
      
      void prepareCoefficientList(FESystem::Field::FieldBase<FESystemDouble>&,
                                  FESystem::Surrogates::PolynomialType t);
      
      const FESystem::Surrogates::PolynomialTerm& getPolynomialTerm(FESystemUInt n) const;
      
      FESystem::Surrogates::PolynomialCoefficientListIterator& getIterator();
      
      FESystemUInt getTotalNumberOfCoefficients() const;

      FESystemUInt getNumberOfActiveCoefficients() const;

      FESystemDouble getCurrentParameter() const;
      
      FESystemBoolean ifKeepPolynomialTerm(FESystemUInt n) const;
            
      void write(std::ostream& out) const;
      
    protected:
      
      void addTerms(FESystemUInt order);

      void addPolynomialTerm(const FESystem::Surrogates::PolynomialTerm& p);
      
      void removeCouplingTerms();

      std::pair<FESystemBoolean, FESystem::Surrogates::PolynomialTerm*>
      findPolynomialTerm(const FESystem::Surrogates::PolynomialTerm& p);
      
      std::vector<FESystem::Surrogates::PolynomialTerm> polynomial_terms;
      
      FESystemUInt number_of_parameters;
      
      FESystem::Surrogates::PolynomialType polynomial_type;
      
      FESystem::Field::FieldBase<FESystemDouble>* field;
      
      FESystemBoolean initialized;
      
      std::auto_ptr<FESystem::Surrogates::PolynomialCoefficientListIterator> polynomial_iterator;
    };
    
    DeclareException0(PolynomialCoefficientListNotInitialized, 
                      << "Polynomial Coefficient List Not Initialized Before Usage\n");

    DeclareException0(PolynomialCoefficientListNotClearedBeforeReInitialization, 
                      << "Polynomial Coefficient List Not Cleared Before Re-initialization\n");

    DeclareException2(PolynomialTermNumberExceedsCurrentCoefficientList, FESystemUInt, FESystemUInt, 
                      << "Given number exceeds the number of terms in current polynomial\n"
                      << "Number of Terms: " << Arg1 << "\n"
                      << "Given number: " << Arg2 << "\n");
    
  }
}


#endif // __fesystem_polynomial_coefficient_list_h__

