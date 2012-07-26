/*
 *  LocalField.h
 *  FESystem
 *
 *  Created by Manav  Bhatia on 1/10/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __fesystem_local_field_h__
#define __fesystem_local_field_h__

// C++ includes
#include <memory>


// FESystem includes
#include "Base/FESystemExceptions.h"
#include "Field/FieldBase.h"

namespace FESystem
{
  // Forward declerations
  namespace Numerics {template <typename ValType> class MatrixBase;}
  namespace Numerics {template <typename ValType> class VectorBase;}
  
  namespace Field
  {
    // Forward declerations
    template <typename ValType> class FieldPointIterator;
    
    template <typename ValType>
    class LocalField: public FESystem::Field::FieldBase<ValType>
    {
    public:
      LocalField();
      
      virtual ~LocalField();
      
      virtual FESystem::Field::FieldPointIterator<ValType>& getFieldPointIterator();
      
      virtual void setDimensions(FESystemUInt n_points, FESystemUInt n_params, FESystemUInt n_resps);
      
      virtual void getParametersAtPoint(FESystemUInt n, FESystem::Numerics::VectorBase<ValType>& vec) const;
      
      virtual void getResponsesAtPoint(FESystemUInt n, FESystem::Numerics::VectorBase<ValType>& vec) const;
      
      virtual void getParameterAtAllPoints(FESystemUInt n, FESystem::Numerics::VectorBase<ValType>& vec) const;
      
      virtual void getResponseAtAllPoints(FESystemUInt n, FESystem::Numerics::VectorBase<ValType>& vec) const;

      virtual void addParameterAndResponses(FESystemUInt point,
                                            const FESystem::Numerics::VectorBase<ValType>& param,
                                            const FESystem::Numerics::VectorBase<ValType>& resps);
      
    protected:
      
      std::auto_ptr<FESystem::Numerics::MatrixBase<ValType> > parameters;
      
      std::auto_ptr<FESystem::Numerics::MatrixBase<ValType> > responses;
      
      std::auto_ptr<FESystem::Field::FieldPointIterator<ValType> > field_parameter_iterator;
    };
  }
}


#endif

