/*
 *  FieldBase.h
 *  FESystem
 *
 *  Created by Manav  Bhatia on 1/10/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __fesystem_field_base_h__
#define __fesystem_field_base_h__

// FESystem includes 
#include "Base/FESystemTypes.h"
#include "Base/FESystemExceptions.h"


namespace FESystem
{
  namespace Numerics { template <typename ValType> class VectorBase; }
  
  namespace Field
  {
    
    // Forward declerations 
    template <typename ValType> class FieldPointIterator;

    template <typename ValType>
    class FieldBase
    {
    public:
      FieldBase();
      
      virtual ~FieldBase();
      
      FESystemUInt getNumberOfParameters() const;
      
      FESystemUInt getNumberOfResponses() const;
      
      FESystemUInt getNumberOfPoints() const;
      
      virtual void setDimensions(FESystemUInt n_points, FESystemUInt n_params, FESystemUInt n_resps);
    
      virtual FESystem::Field::FieldPointIterator<ValType>& getFieldPointIterator()=0;

      virtual void getParametersAtPoint(FESystemUInt n, FESystem::Numerics::VectorBase<ValType>& vec) const = 0;

      virtual void getResponsesAtPoint(FESystemUInt n, FESystem::Numerics::VectorBase<ValType>& vec) const = 0;

      virtual void getParameterAtAllPoints(FESystemUInt n, FESystem::Numerics::VectorBase<ValType>& vec) const = 0;
      
      virtual void getResponseAtAllPoints(FESystemUInt n, FESystem::Numerics::VectorBase<ValType>& vec) const = 0;

      virtual void addParameterAndResponses(FESystemUInt point,
                                            const FESystem::Numerics::VectorBase<ValType>& param,
                                            const FESystem::Numerics::VectorBase<ValType>& resps)=0;
      
    protected:
      
      FESystemUInt n_parameters;
      
      FESystemUInt n_responses;

      FESystemUInt n_points;
    };
    
    DeclareException2(PointNumberExceedsLimit, FESystemUInt, FESystemUInt,
                      << "Given point number is greater than number of points in field.\n"
                      << "Number of points: " << Arg1 << "\n"
                      << "Given point number: " << Arg2 << "\n");
    
    DeclareException2(VectorSizeMismatch, FESystemUInt, FESystemUInt,
                      << "Given vector does not match needed size.\n"
                      << "Required size: " << Arg1 << "\n"
                      << "Given size: " << Arg2 << "\n");
        
    DeclareException2(NumberOfParameterDifferentFromFieldSize, FESystemUInt, FESystemUInt,
                      << "Number of parameters not equal to size of this field.\n"
                      << "Number of parameters: " << Arg1 << "\n"
                      << "Given number: " << Arg2 << "\n");
    
    DeclareException2(NumberOfRespnsesDifferentFromFieldSize, FESystemUInt, FESystemUInt,
                      << "Number of responses not equal to size of this field.\n"
                      << "Number of responses: " << Arg1 << "\n"
                      << "Given number: " << Arg2 << "\n");
  }
}

#endif // __fesystem_field_base_h__
