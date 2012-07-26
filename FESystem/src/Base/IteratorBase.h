/*
 *  IteratorBase.h
 *  FESystem
 *
 *  Created by Manav  Bhatia on 2/16/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __fesystem_iterator_base_h__
#define __fesystem_iterator_base_h__

// FESystem includes
#include "Base/FESystemTypes.h"

namespace FESystem
{
  namespace Base
  {
    class IteratorBase
    {
    public:

      virtual void increment() = 0;
      
      virtual void decrement() = 0;
      
      virtual void begin() = 0;
      
      virtual void reverseBegin() = 0;
      
      virtual FESystemBoolean ifBegin() const = 0;
      
      virtual FESystemBoolean ifReverseEnd() const = 0;
      
      virtual FESystemBoolean ifEnd() const = 0;
            
    protected:
      
    };
  }
}

#endif // __fesystem_field_parameter_iterator_h__
