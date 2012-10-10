/*
 *  PolynomialCoefficientListIterator.h
 *  FESystem
 *
 *  Created by Manav  Bhatia on 1/10/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __fesystem_polynomial_coefficient_list_iterator_h__
#define __fesystem_polynomial_coefficient_list_iterator_h__

// FESystem includes
#include "Base/FESystemTypes.h"
#include "Base/IteratorBase.h"

namespace FESystem
{
  // Forward declerations
  namespace Numerics { template <typename ValType> class VectorBase;}
  
  namespace Surrogates
  {
    // Forward declerations
    class PolynomialCoefficientList;

    
    class PolynomialCoefficientListIterator
    {
    public:
      PolynomialCoefficientListIterator(const FESystem::Surrogates::PolynomialCoefficientList& c_list);
      
      virtual ~PolynomialCoefficientListIterator();
      
      virtual void increment();
      
      virtual void decrement();
      
      virtual void begin();
      
      virtual void reverseBegin();
      
      virtual FESystemBoolean ifBegin() const;
      
      virtual FESystemBoolean ifReverseEnd() const;
      
      virtual FESystemBoolean ifEnd() const;
      
      FESystemInt getID() const;
      
      template <typename ValType>
      ValType getCurrentValue(const FESystem::Numerics::VectorBase<ValType>& vec);
      
    protected:
      
      const FESystem::Surrogates::PolynomialCoefficientList& coeff_list;

      FESystemInt active_term_counter;
      
      FESystemInt counter;
    };
    
  }
}


#endif // __fesystem_polynomial_coefficient_list_iterator_h__
