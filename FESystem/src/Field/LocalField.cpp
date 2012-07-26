/*
 *  LocalField.cpp
 *  FESystem
 *
 *  Created by Manav  Bhatia on 1/10/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

// FESystem includes
#include "Field/LocalField.h"
#include "Field/FieldPointIterator.h"
#include "Numerics/LocalVector.h"
#include "Numerics/DenseMatrix.h"


template <typename ValType>
FESystem::Field::LocalField<ValType>::LocalField():
FESystem::Field::FieldBase<ValType>()
{
  this->field_parameter_iterator.reset(new FESystem::Field::FieldPointIterator<ValType>(*this));
  
  this->parameters.reset(FESystem::Numerics::MatrixCreate<ValType>
                         (FESystem::Numerics::LOCAL_DENSE_MATRIX).release());
  this->responses.reset(FESystem::Numerics::MatrixCreate<ValType>
                        (FESystem::Numerics::LOCAL_DENSE_MATRIX).release());
}

template <typename ValType>
FESystem::Field::LocalField<ValType>::~LocalField()
{
  
}




template <typename ValType>
void
FESystem::Field::LocalField<ValType>::setDimensions(FESystemUInt n_points,
                                                    FESystemUInt n_params,
                                                    FESystemUInt n_resps)
{
  FESystem::Field::FieldBase<ValType>::setDimensions(n_points, n_params, n_resps);
  
  this->parameters->resize(n_points, n_params);
  this->responses->resize(n_points, n_resps);
  
  this->field_parameter_iterator->initialize();
}




template <typename ValType>
FESystem::Field::FieldPointIterator<ValType>&
FESystem::Field::LocalField<ValType>::getFieldPointIterator()
{
  return *(this->field_parameter_iterator);
}



template <typename ValType>
void 
FESystem::Field::LocalField<ValType>::addParameterAndResponses
(FESystemUInt point_num,
 const FESystem::Numerics::VectorBase<ValType>& param,
 const FESystem::Numerics::VectorBase<ValType>& resps)
{
  FESystemAssert2(point_num < this->getNumberOfPoints(),
                  FESystem::Field::PointNumberExceedsLimit,
                  this->getNumberOfPoints(), point_num);
  FESystemAssert2(param.getSize() == this->getNumberOfParameters(),
                  FESystem::Field::NumberOfParameterDifferentFromFieldSize,
                  this->getNumberOfParameters(), param.getSize());
  FESystemAssert2(resps.getSize() == this->getNumberOfResponses(),
                  FESystem::Field::NumberOfRespnsesDifferentFromFieldSize,
                  this->getNumberOfResponses(), resps.getSize());
  
  // add the point 
  this->parameters->setRowVals(point_num, 0, param.getSize()-1, param);
  this->responses->setRowVals(point_num, 0, resps.getSize()-1, resps);
}




template <typename ValType>
void 
FESystem::Field::LocalField<ValType>::getParametersAtPoint(FESystemUInt n, 
                                                           FESystem::Numerics::VectorBase<ValType>& vec) const
{
  FESystemAssert2(n < this->getNumberOfPoints(), 
                  FESystem::Field::PointNumberExceedsLimit,
                  this->getNumberOfPoints(), n);
  FESystemAssert2(vec.getSize() == this->getNumberOfParameters(),
                  FESystem::Field::VectorSizeMismatch,
                  this->getNumberOfParameters(), vec.getSize());
  
  this->parameters->getRowVals(n, 0, this->getNumberOfParameters()-1, vec);
}



template <typename ValType>
void
FESystem::Field::LocalField<ValType>::getResponsesAtPoint(FESystemUInt n, 
                                                          FESystem::Numerics::VectorBase<ValType>& vec) const
{
  FESystemAssert2(n < this->getNumberOfPoints(), 
                  FESystem::Field::PointNumberExceedsLimit,
                  this->getNumberOfPoints(), n);
  FESystemAssert2(vec.getSize() == this->getNumberOfResponses(), 
                  FESystem::Field::VectorSizeMismatch,
                  this->getNumberOfResponses(), vec.getSize());
  
  this->responses->getRowVals(n, 0, this->getNumberOfResponses()-1, vec);
}


template <typename ValType>
void 
FESystem::Field::LocalField<ValType>::getParameterAtAllPoints(FESystemUInt n, 
                                                              FESystem::Numerics::VectorBase<ValType>& vec) const
{
  FESystemAssert2(n < this->getNumberOfParameters(), 
                  FESystem::Field::NumberOfParameterDifferentFromFieldSize,
                  this->getNumberOfParameters(), n);
  FESystemAssert2(vec.getSize() == this->getNumberOfPoints(),
                  FESystem::Field::VectorSizeMismatch,
                  this->getNumberOfPoints(), vec.getSize());
  
  this->parameters->getColumnVals(n, 0, this->getNumberOfPoints()-1, vec);
}



template <typename ValType>
void
FESystem::Field::LocalField<ValType>::getResponseAtAllPoints(FESystemUInt n, 
                                                             FESystem::Numerics::VectorBase<ValType>& vec) const
{
  FESystemAssert2(n < this->getNumberOfResponses(), 
                  FESystem::Field::NumberOfRespnsesDifferentFromFieldSize,
                  this->getNumberOfResponses(), n);
  FESystemAssert2(vec.getSize() == this->getNumberOfPoints(), 
                  FESystem::Field::VectorSizeMismatch,
                  this->getNumberOfPoints(), vec.getSize());
  
  this->responses->getColumnVals(n, 0, this->getNumberOfPoints()-1, vec);
}


/***************************************************************************************/
// Template instantiations for some generic classes

template class FESystem::Field::LocalField<FESystemDouble>;


/***************************************************************************************/
