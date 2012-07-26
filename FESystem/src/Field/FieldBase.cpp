/*
 *  FieldBase.cpp
 *  FESystem
 *
 *  Created by Manav  Bhatia on 1/10/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "Field/FieldBase.h"



template <typename ValType>
FESystem::Field::FieldBase<ValType>::FieldBase():
n_parameters(0),
n_responses(0),
n_points(0)
{
  
}
      

template <typename ValType>
FESystem::Field::FieldBase<ValType>::~FieldBase()
{
  
}
      

template <typename ValType>
FESystemUInt 
FESystem::Field::FieldBase<ValType>::getNumberOfParameters() const
{
  return this->n_parameters;
}
      

template <typename ValType>
FESystemUInt 
FESystem::Field::FieldBase<ValType>::getNumberOfResponses() const
{
  return this->n_responses;
}


template <typename ValType>
FESystemUInt 
FESystem::Field::FieldBase<ValType>::getNumberOfPoints() const
{
  return this->n_points;
}


template <typename ValType>
void
FESystem::Field::FieldBase<ValType>::setDimensions(FESystemUInt points,
                                                   FESystemUInt params, 
                                                   FESystemUInt resps)
{
  this->n_points = points;
  this->n_parameters = params;
  this->n_responses = resps;
}
      


/***************************************************************************************/
// Template instantiations for some generic classes

template class FESystem::Field::FieldBase<FESystemDouble>;


/***************************************************************************************/

