//
//  LinearFunctionMapping.h
//  FESystem
//
//  Created by Manav Bhatia on 4/13/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem_linear_function_mapping_h__
#define __fesystem_linear_function_mapping_h__

// C++ includes
#include <memory>

// FESystem includes
#include "Functions/FunctionMappingBase.h"

namespace FESystem
{
    namespace Functions
    {
        /*!
         *   This specializes the FunctionMappingBase class for linear functions for a specified basis of 
         *   transformation (Jacobian). 
         */
        template <typename ValType>
        class LinearFunctionMapping: public FESystem::Functions::FunctionMappingBase<ValType>
        {
        public:
            /*!
             *   Constructor
             */
            LinearFunctionMapping();
            
            virtual ~LinearFunctionMapping();
            
            /*!
             *   clears the data structures use with another mapping
             */
            void clear();

            /*!
             *   initializes the data structures for a specified Jacobian
             */
            void reinit(const FESystem::Numerics::MatrixBase<ValType>& mat);
            
            /*!
             *   Returns the dimensions of the original and the mapped space
             */
            virtual std::pair<FESystemUInt, FESystemUInt> getDimensions() const; 

            /*!
             *   Maps the input vector \p vin to the mapped space vector \p vout. The user must ensure that the vector 
             *   dimensions are consistent
             */
            virtual void map(const FESystem::Numerics::VectorBase<ValType>& vin, FESystem::Numerics::VectorBase<ValType>& vout) const;
                        
            /*!
             *   Inverse maps the input vector \p vin to the original space vector \p vout. The user must ensure that the vector 
             *   dimensions are consistent
             */
            virtual void inverseMap(const FESystem::Numerics::VectorBase<ValType>& vin, FESystem::Numerics::VectorBase<ValType>& vout) const;

            /*!
             *   Calculates the Jacobian of the mapping at the specified point \p vin and returns it in the matrix \p mat. These are essentially
             *   the first order derivatives of the mapped functions with respect to the original variables. 
             */
            virtual void getMappingJacobian(const FESystem::Numerics::VectorBase<ValType>& vin, FESystem::Numerics::MatrixBase<ValType>& mat) const;
            
            /*!
             *   Returns the value of derivatives of the mapped function at \p val in X1, with respect to the derivative information provided 
             *   in \p derivative_orders, whose i^th element is the order of the derivative with respect to the i^th coordinate. The derivative 
             *   values are returned in \p deriv.
             */
            virtual void getFunctionDerivative(const std::vector<FESystemUInt>& derivative_orders, const FESystem::Numerics::VectorBase<ValType>& val, 
                                               FESystem::Numerics::VectorBase<ValType>& deriv) const;
            
        protected:
            
            /*!
             *   Stores if the mapping is initialized
             */
            FESystemBoolean if_initialized;
            
            /*!
             *   Transformation Jacobian
             */
            std::auto_ptr<FESystem::Numerics::MatrixBase<ValType> > jac;
            
        };
    }
}

#endif // __fesystem_linear_function_mapping_h__

