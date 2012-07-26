//
//  FunctionMappingBase.h
//  FESystem
//
//  Created by Manav Bhatia on 4/13/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem_function_mapping_base_h__
#define __fesystem_function_mapping_base_h__

// C++ includes
#include <vector>

// FESystem includes
#include "Base/FESystemTypes.h"
#include "Base/FESystemExceptions.h"

namespace FESystem
{
    // Forward declerations
    namespace Numerics {template <typename ValType> class VectorBase;}
    namespace Numerics {template <typename ValType> class MatrixBase;}
    
    namespace Functions
    {
        /*!
         *   This is a base class for creation of function mappings. The basic idea of mapping is that 
         *   if X is an R^(m1 x n1) coordinate, then the function F(X) is a mapped R^(m2 x n2) coordinate such that given a 
         *   value of X, F(X) and its derivatives can be calculated. This can be used to create coordinate system 
         *   basis, finite element shape functions, etc. The function is specialized for n1=n2=1 with vector data types, 
         *   and for general matrix types. 
         */
        template <typename ValType>
        class FunctionMappingBase
        {
        public:
            /*!
             *   Constructor
             */
            FunctionMappingBase();
            
            virtual ~FunctionMappingBase();

            /*!
             *   Returns the dimensions of the original and the mapped space
             */
            virtual std::pair<FESystemUInt, FESystemUInt> getDimensions() const = 0; 
            
            /*!
             *   Maps the input vector \p vin to the mapped space vector \p vout. The user must ensure that the vector 
             *   dimensions are consistent
             */
            virtual void map(const FESystem::Numerics::VectorBase<ValType>& vin, FESystem::Numerics::VectorBase<ValType>& vout) const=0;
            
            /*!
             *   Inverse maps the input vector \p vin to the original space vector \p vout. The user must ensure that the vector 
             *   dimensions are consistent
             */
            virtual void inverseMap(const FESystem::Numerics::VectorBase<ValType>& vin, FESystem::Numerics::VectorBase<ValType>& vout) const=0;
            
            /*!
             *   Calculates the Jacobian of the mapping at the specified point \p vin and returns it in the matrix \p mat. These are essentially
             *   the first order derivatives of the mapped functions with respect to the original variables. 
             */
            virtual void getMappingJacobian(const FESystem::Numerics::VectorBase<ValType>& vin, FESystem::Numerics::MatrixBase<ValType>& mat) const=0;
            
            /*!
             *   Returns the value of derivatives of the mapped function at \p val in X1, with respect to the derivative information provided 
             *   in \p derivative_orders, whose i^th element is the order of the derivative with respect to the i^th coordinate. The derivative 
             *   values are returned in \p deriv.
             */
            virtual void getFunctionDerivative(const std::vector<FESystemUInt>& derivative_orders, const FESystem::Numerics::VectorBase<ValType>& val, FESystem::Numerics::VectorBase<ValType>& deriv) const=0;

            
        protected:
            
        };
        
        
        /*!
         *   This defines a specialized class for discrete mapping functions where a separate interpolation functions exist for degrees of freedoms,
         *   such as the Lagrange interpolation function. The mapping function provides an additional method to return these interpolation functions
         *   in a vector form (as opposed to proving the interpolated quantity). 
         */
        template <typename ValType>
        class DiscreteFunctionMappingBase: public FESystem::Functions::FunctionMappingBase<ValType>
        {
        public:
            /*!
             *   Constructor
             */
            DiscreteFunctionMappingBase();
            
            ~DiscreteFunctionMappingBase();
            
            /*!
             *   Returns the number of discrete functions, which is the number of shape functions that this method will return. 
             */
            virtual FESystemUInt getNDiscreteFunctions() const=0;
            
            /*!
             *   This is an invalid function call for this method, and an error will be thrown if this method is called
             */
            virtual void inverseMap(const FESystem::Numerics::VectorBase<ValType>& vin, FESystem::Numerics::VectorBase<ValType>& vout) const;

            /*!
             *   Sets the values of the discrete function in the mapped domain so that the discrete mapping function can perform the interpolation. 
             *   The individual discrete function mapping derived from this method will specify the ordering of the dofs that should be ensured in 
             *   the input value. THe value of i_dim should be less than or equal to the value of the dimension of the mapped domain. 
             */
            virtual void setDiscreteFunctionValues(const FESystem::Numerics::MatrixBase<ValType>& vals);

            /*!
             *   Returns the value of the discrete interpolation functions.
             */
            virtual void getDiscreteFunction(const FESystem::Numerics::VectorBase<ValType>& val, FESystem::Numerics::VectorBase<ValType>& deriv) const=0;

            
            /*!
             *   Returns the derivative of the discrete interpolation functions. The derivative order with respect to the i^th coordinate is specified in the 
             *   i^th element of \p derivative_order. The location at which this is to be evaluated is \p val. The dimension of \p deriv is equal the number 
             *   of interpolation function values. 
             */
            virtual void getDiscreteFunctionDerivative(const std::vector<FESystemUInt>& derivative_orders, const FESystem::Numerics::VectorBase<ValType>& val, 
                                                       FESystem::Numerics::VectorBase<ValType>& deriv) const=0;

            
        protected:
            
            /*!
             *   clears the data structures
             */
            virtual void clear();
            
            /*!
             *   initializes the data structures
             */
            virtual void reinit(FESystemUInt dim, FESystemUInt n_pts);
            
            /*!
             *   stores the state of initialization of discrete function mapping
             */
            FESystemBoolean if_discrete_function_initialized;
            
            /*!
             *   This stores the values of the mapped domain discrete function values
             */
            FESystem::Numerics::MatrixBase<ValType>* discrete_function_values; 
            
        };

        /*!
         *   makes sure that the dimensions are consistent
         */
        DeclareException2(InconsistentDimensions, 
                          FESystemUInt,
                          FESystemUInt,
                          << "Found : " << Arg1 
                          << "Should be: " << Arg2 );
        
    }
}


#endif //__fesystem_function_mapping_base_h__

