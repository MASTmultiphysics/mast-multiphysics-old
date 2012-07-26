/*
 *  FieldPointIterator.h
 *  FESystem
 *
 *  Created by Manav  Bhatia on 1/10/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __fesystem_field_parameter_iterator_h__
#define __fesystem_field_parameter_iterator_h__

// C++ includes
#include <memory>

// FESystem includes
#include "Base/FESystemTypes.h"
#include "Base/IteratorBase.h"

namespace FESystem
{
  // Forward declerations
  namespace Numerics { template <typename ValType> class VectorBase;}

  namespace Field
  {
    // Forward declerations
    template <typename ValType> class FieldBase;
    
    template <typename ValType>
    class FieldPointIterator: FESystem::Base::IteratorBase
    {
    public:
      FieldPointIterator(const FESystem::Field::FieldBase<ValType>& f);
      
      virtual ~FieldPointIterator();
      
      void initialize();

      virtual void increment();
      
      virtual void decrement();
      
      virtual void begin();
      
      virtual void reverseBegin();
      
      virtual FESystemBoolean ifBegin() const;
      
      virtual FESystemBoolean ifReverseEnd() const;
      
      virtual FESystemBoolean ifEnd() const;
      
      FESystemUInt getCurrentPointID() const;
      
      const FESystem::Numerics::VectorBase<ValType>& getCurrentParameter();
      
      const FESystem::Numerics::VectorBase<ValType>& getResponsesAtPoint();

      FESystemInt getID() const;
      
    protected:
      
      
      const FESystem::Field::FieldBase<ValType>& field;
      
      std::auto_ptr<FESystem::Numerics::VectorBase<ValType> > param;
      
      std::auto_ptr<FESystem::Numerics::VectorBase<ValType> > resp;
      
      FESystemInt counter;
    };
  }
}

#endif // __fesystem_field_parameter_iterator_h__
