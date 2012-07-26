/*
 *  SurrogateBase.cpp
 *  FESystem
 *
 *  Created by Manav Bhatia on 1/9/11.
 *  Copyright 2011 . All rights reserved.
 *
 */

#include "Surrogates/SurrogateBase.h"

template <typename ValType>
FESystem::Surrogates::SurrogateBase<ValType>::SurrogateBase():
field(NULL)
{
  
}
      
template <typename ValType>
FESystem::Surrogates::SurrogateBase<ValType>::~SurrogateBase()
{
  
}

      
template <typename ValType>
void 
FESystem::Surrogates::SurrogateBase<ValType>::setField(FESystem::Field::FieldBase<ValType>& f)
{
  FESystemAssert0(this->field == NULL,
                  FESystem::Surrogates::SurrogateFieldAlreadyPresent);
  
  this->field = &f;
}
      

template <typename ValType>
FESystem::Field::FieldBase<ValType>&
FESystem::Surrogates::SurrogateBase<ValType>::getField()
{
  FESystemAssert0(this->field != NULL,
                  FESystem::Surrogates::SurrogateFieldIsNull);
  
  return (*(this->field));
}
      

/***************************************************************************************/
// Template instantiations for some generic classes

template class FESystem::Surrogates::SurrogateBase<FESystemDouble>;


/***************************************************************************************/
