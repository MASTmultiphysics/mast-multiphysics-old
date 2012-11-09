/*
 *  LeastSquareSurrogate.cpp
 *  FESystem
 *
 *  Created by Manav  Bhatia on 1/10/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

// FESystem includes
#include "Surrogates/LeastSquareSurrogate.h"
#include "Field/FieldBase.h"
#include "Field/FieldPointIterator.h"
#include "Solvers/LinearSolvers/LapackLinearSolver.h"
#include "Solvers/LinearSolvers/LinearLeastSquareSolver.h"
#include "Numerics/DenseMatrix.h"
#include "Numerics/LocalVector.h"
#include "Surrogates/PolynomialCoefficientList.h"
#include "Surrogates/PolynomialCoefficientListIterator.h"


template <typename ValType>
FESystem::Surrogates::LeastSquareSurrogate<ValType>::LeastSquareSurrogate():
FESystem::Surrogates::SurrogateBase<ValType>(),
initialized(false)
{
  
}



template <typename ValType>
FESystem::Surrogates::LeastSquareSurrogate<ValType>::~LeastSquareSurrogate()
{
  
}



template <typename ValType>
void 
FESystem::Surrogates::LeastSquareSurrogate<ValType>::initialize()
{
  FESystem::Field::FieldBase<ValType>& f = this->getField();
  FESystemUInt n_resps = f.getNumberOfResponses(), n_points = f.getNumberOfPoints(); 
    
  // create the least square matrix
  this->ls_matrix.reset(FESystem::Numerics::MatrixCreate<ValType>
                        (FESystem::Numerics::LOCAL_DENSE_MATRIX).release());

  // create the coefficients
  this->coeffs.reset(new FESystem::Surrogates::PolynomialCoefficientList());
  this->coeffs->prepareCoefficientList(f, this->getPolynomialType());
  FESystem::Surrogates::PolynomialCoefficientListIterator& coeff_it = this->coeffs->getIterator();
  this->coeffs->write(std::cout);

  // now iterate over the parameters and create the matrix
  FESystem::Field::FieldPointIterator<ValType>& point_it = f.getFieldPointIterator();
  
  this->ls_matrix->resize(n_points, this->coeffs->getNumberOfActiveCoefficients());
  
  
  for (point_it.begin(); !point_it.ifEnd(); point_it.increment())
    for (coeff_it.begin(); !coeff_it.ifEnd(); coeff_it.increment())
      this->ls_matrix->setVal(point_it.getCurrentPointID(), coeff_it.getID(), 
                               coeff_it.getCurrentValue(point_it.getCurrentParameter()));
  
  // create the solver
  this->least_square_solver.reset(new FESystem::LinearSolvers::LinearLeastSquareSolver<ValType>());
  this->linear_solver.reset(new FESystem::LinearSolvers::LapackLinearSolver<ValType>());
  this->ls_solved_polynomial_coefficient.resize(n_resps);
  for (FESystemUInt i=0; i<n_resps; i++)
  {
    this->ls_solved_polynomial_coefficient[i] = 
    FESystem::Numerics::VectorCreate<ValType>(FESystem::Numerics::LOCAL_VECTOR).release();
    this->ls_solved_polynomial_coefficient[i]->resize(this->coeffs->getNumberOfActiveCoefficients()); 
  }

  // attach the matrix and solver
  this->least_square_solver->setSystemMatrix(*(this->ls_matrix), false);
  this->least_square_solver->setLinearSolver(*(this->linear_solver));

  // calculate rhs
  std::auto_ptr<FESystem::Numerics::VectorBase<ValType> > 
  rhs(FESystem::Numerics::VectorCreate<ValType>(FESystem::Numerics::LOCAL_VECTOR).release());
  rhs->resize(n_points);

  for (FESystemUInt i=0; i<n_resps; i++)
  {
    // zero the vectors
    rhs->zero();
    this->ls_solved_polynomial_coefficient[i]->zero();
    
    // get the RHS for the least square system and solve. The multiplication is performed
    // by the LS solver. 
    f.getResponseAtAllPoints(i, *rhs); 
    this->least_square_solver->solve(*rhs, *(this->ls_solved_polynomial_coefficient[i])); 
  }
  
  // set the flag to true
  this->initialized = true;
}



template <typename ValType>
ValType 
FESystem::Surrogates::LeastSquareSurrogate<ValType>::getEstimate
(FESystemUInt resp_num, const FESystem::Numerics::VectorBase<ValType>& vec)
{
  FESystemAssert0( this->initialized, SurrogateNotInitialized);
  
  FESystemAssert2(resp_num < this->getField().getNumberOfResponses(),
                  FESystem::Surrogates::ResponseNumberExceedsCurrentField,
                  this->getField().getNumberOfResponses(), resp_num);
  
  // calculate the value using the polynomial coefficients calculated through the solver
  FESystem::Surrogates::PolynomialCoefficientListIterator& coeff_it = this->coeffs->getIterator();
  ValType val=0.0;
  for (coeff_it.begin(); !coeff_it.ifEnd(); coeff_it.increment())
    val += coeff_it.getCurrentValue(vec) * this->ls_solved_polynomial_coefficient[resp_num]->getVal(coeff_it.getID());
  
  return val;
}


template <typename ValType>
void
FESystem::Surrogates::LeastSquareSurrogate<ValType>::setPolynomialType
(FESystem::Surrogates::PolynomialType p)
{
  this->polynomial_type = p;
}


template <typename ValType>
FESystem::Surrogates::PolynomialType
FESystem::Surrogates::LeastSquareSurrogate<ValType>::getPolynomialType() const
{
  return this->polynomial_type;
}



/***************************************************************************************/
// Template instantiations for some generic classes

template class FESystem::Surrogates::LeastSquareSurrogate<double>;


/***************************************************************************************/

