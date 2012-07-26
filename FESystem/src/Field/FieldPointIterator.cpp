/*
 *  FieldPointIterator.cpp
 *  FESystem
 *
 *  Created by Manav  Bhatia on 1/10/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */


// FESystem includes
#include "Field/FieldPointIterator.h"
#include "Field/FieldBase.h"
#include "Numerics/VectorBase.h"


template <typename ValType>
FESystem::Field::FieldPointIterator<ValType>::FieldPointIterator
(const FESystem::Field::FieldBase<ValType>& f):
field(f),
counter(0)
{
  this->param.reset(FESystem::Numerics::VectorCreate<ValType>
                    (FESystem::Numerics::LOCAL_VECTOR).release());
  this->resp.reset(FESystem::Numerics::VectorCreate<ValType>
                   (FESystem::Numerics::LOCAL_VECTOR).release());
}


template <typename ValType>
FESystem::Field::FieldPointIterator<ValType>::~FieldPointIterator()
{
  
}
      

template <typename ValType>
void 
FESystem::Field::FieldPointIterator<ValType>::initialize()
{
  this->param->resize(this->field.getNumberOfParameters());
  this->resp->resize(this->field.getNumberOfResponses());
}


template <typename ValType>
void 
FESystem::Field::FieldPointIterator<ValType>::increment()
{
  this->counter++;
}
  

template <typename ValType>
void 
FESystem::Field::FieldPointIterator<ValType>::decrement()
{
  this->counter--;
}
      
template <typename ValType>
void 
FESystem::Field::FieldPointIterator<ValType>::begin()
{
  this->counter = 0;
}
      

template <typename ValType>
void 
FESystem::Field::FieldPointIterator<ValType>::reverseBegin()
{
  this->counter = this->field.getNumberOfPoints()-1;
}
      
template <typename ValType>
FESystemBoolean 
FESystem::Field::FieldPointIterator<ValType>::ifBegin() const
{
  return (this->counter == 0);
}


template <typename ValType>
FESystemBoolean 
FESystem::Field::FieldPointIterator<ValType>::ifReverseEnd() const
{
  return (this->counter < 0);
}


template <typename ValType>
FESystemBoolean 
FESystem::Field::FieldPointIterator<ValType>::ifEnd() const
{
  return (this->counter > this->field.getNumberOfPoints()-1);
}
      

template <typename ValType>
FESystemUInt 
FESystem::Field::FieldPointIterator<ValType>::getCurrentPointID() const
{
  return this->counter;
}
      

template <typename ValType>
const FESystem::Numerics::VectorBase<ValType>& 
FESystem::Field::FieldPointIterator<ValType>::getCurrentParameter() 
{
  this->field.getParametersAtPoint(this->getCurrentPointID(), *(this->param));
  return *(this->param);
}
      

template <typename ValType>
const FESystem::Numerics::VectorBase<ValType>& 
FESystem::Field::FieldPointIterator<ValType>::getResponsesAtPoint() 
{
  this->field.getResponsesAtPoint(this->getCurrentPointID(), *(this->resp));
  return *(this->resp);
}
      
      

/***************************************************************************************/
// Template instantiations for some generic classes

template class FESystem::Field::FieldPointIterator<FESystemDouble>;


/***************************************************************************************/
