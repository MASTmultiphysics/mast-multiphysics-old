/*
 *  SurrogateBase.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 1/9/11.
 *  Copyright 2011 . All rights reserved.
 *
 */

#ifndef __fesystem_surrogate_base_h__
#define __fesystem_surrogate_base_h__

// FESystem includes
#include "Base/FESystemTypes.h"
#include "Base/FESystemExceptions.h"

namespace FESystem
{
  // Forward declerations
  namespace Field {template <typename ValType> class FieldBase; }
  namespace Numerics { template <typename ValType> class VectorBase; }
  
  namespace Surrogates
  {
    template <typename ValType>
    class SurrogateBase
    {
    public:
      SurrogateBase();
      
      virtual ~SurrogateBase();
      
      void setField(FESystem::Field::FieldBase<ValType>& f);
      
      FESystem::Field::FieldBase<ValType>& getField();
      
      virtual void initialize() = 0;
      
      virtual ValType getEstimate(FESystemUInt resp_num, 
                                  const FESystem::Numerics::VectorBase<ValType>& vec) = 0;
      
    protected:
      
      
      FESystem::Field::FieldBase<ValType>* field;
      
    };
    
    DeclareException0(SurrogateFieldAlreadyPresent,
                      << "Field pointer already present. Cannot be reset unless cleared.");
    
    DeclareException0(SurrogateFieldIsNull,
                      << "Field pointer is NULL. Cannot be accessed.");
    
    DeclareException0(SurrogateNotInitialized,
                      << "Surrogate is uninitialized. Cannot proceed with calculations.");
    
    DeclareException2(ResponseNumberExceedsCurrentField, FESystemUInt, FESystemUInt,
                      << "Requested response number exceeds current limit.\n" 
                      << "Number of responses: " << Arg1 << "\n"
                      << "Given number: " << Arg2 << "\n");
  }
}



#endif // __fesystem_surrogate_base_h__

