//
//  CylindricalFunctionMapping.h
//  FESystem
//
//  Created by Manav Bhatia on 4/13/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem_polar_function_mapping_h__
#define __fesystem_polar_function_mapping_h__

// FESystem includes
#include "Functions/FunctionMappingBase.h"

namespace FESystem
{
    namespace Functions
    {
        /*!
         *   This specializes the FunctionMappingBase class for cylindrical functions for a specified basis of 
         *   transformation (Jacobian). 
         */
        template <typename ValType>
        class CylindricalFunctionMapping: public FESystem::Functions::FunctionMappingBase<ValType>
        {
        public:
            /*!
             *   Constructor
             */
            CylindricalFunctionMapping();
            
            virtual ~CylindricalFunctionMapping();
            
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
             *   Calculates the Jacobian of the mapping at the specified point \p vin and returns it in the matrix \p mat
             */
            virtual void getMappingJacobian(const FESystem::Numerics::VectorBase<ValType>& vin, FESystem::Numerics::MatrixBase<ValType>& mat) const;
            
        protected:
            
        };
    }
}


#endif // __fesystem_cylindrical_function_mapping_h__
