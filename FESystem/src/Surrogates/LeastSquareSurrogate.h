/*
 *  LeastSquareSurrogate.h
 *  FESystem
 *
 *  Created by Manav  Bhatia on 1/10/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __fesystem_least_square_surrogate_h__
#define __fesystem_least_square_surrogate_h__

// C++ includes
#include <memory>
#include <vector>


// FESystem includes
#include "Surrogates/SurrogateBase.h"
#include "Surrogates/PolynomialType.h"

namespace FESystem
{
  namespace LinearSolvers { template <typename ValType>  class LinearLeastSquareSolver; }
  namespace LinearSolvers { template <typename ValType>  class LinearSolverBase; }
  namespace Numerics { template <typename ValType>  class MatrixBase; }
  
  namespace Surrogates 
  {
    // Forward declerations
    class PolynomialCoefficientList;
    
    template <typename ValType>
    class LeastSquareSurrogate: public FESystem::Surrogates::SurrogateBase<ValType>
    {
    public: 
      LeastSquareSurrogate();
      
      virtual ~LeastSquareSurrogate();
      
      void initialize();
      
      ValType getEstimate(FESystemUInt resp_num, 
                          const FESystem::Numerics::VectorBase<ValType>& vec);  
      
      void setPolynomialType(FESystem::Surrogates::PolynomialType p); 

      FESystem::Surrogates::PolynomialType getPolynomialType() const; 
      
    protected:
      
      std::auto_ptr<FESystem::LinearSolvers::LinearLeastSquareSolver<ValType> > least_square_solver;

      std::auto_ptr<FESystem::LinearSolvers::LinearSolverBase<ValType> > linear_solver;

      std::auto_ptr<FESystem::Numerics::MatrixBase<ValType> > ls_matrix;

      std::auto_ptr<FESystem::Surrogates::PolynomialCoefficientList> coeffs;

      std::vector<FESystem::Numerics::VectorBase<ValType>* > ls_solved_polynomial_coefficient;
      
      FESystem::Surrogates::PolynomialType polynomial_type;
      
      FESystemBoolean initialized;
    };
  }
}


#endif // __fesystem_least_square_surrogate_h__
