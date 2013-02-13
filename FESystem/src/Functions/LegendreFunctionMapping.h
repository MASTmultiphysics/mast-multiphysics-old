//
//  LegendreFunctionMapping.h
//  FESystem
//
//  Created by Manav Bhatia on 2/12/13.
//
//

#ifndef __FESystem__LegendreFunctionMapping__
#define __FESystem__LegendreFunctionMapping__


// FESystem includes
#include "Functions/FunctionMappingBase.h"


namespace FESystem
{
    // Forward declerations
    namespace Functions {template <typename ValType> class LegendreFunction;}
    
    namespace Functions
    {
        /*!
         *   This class derives from the base class DiscreteFunctionMappingBase, and implements the Legendre function
         *   with abcissa and ordinate types as FESystemDoubles.
         */
        template <typename ValType>
        class LegendreFunctionMapping: public FESystem::Functions::DiscreteFunctionMappingBase<ValType>
        {
        public:
            /*!
             *    default constructor
             */
            LegendreFunctionMapping();
            
            /*!
             *    destructor
             */
            virtual ~LegendreFunctionMapping();
            
            /*!
             *   clears the data structures
             */
            virtual void clear();
            
            /*!
             *   This initializes the mapping using the given dimension, \p dim, for the mapped space, and the
             *   Legendre functions for the domain. One Legendre function per domain is needed.
             */
            void reinit(const FESystemUInt dim, const FESystemUInt o, const FESystem::Functions::LegendreFunction<ValType>& func);
            
            /*!
             *   Returns the dimensions of the original and the mapped space
             */
            virtual std::pair<FESystemUInt, FESystemUInt> getDimensions() const;
            
            /*!
             *   Returns the number of discrete functions, which is the number of shape functions that this method will return.
             */
            virtual FESystemUInt getNDiscreteFunctions() const;
            
            /*!
             *   Maps the input vector \p vin to the mapped space vector \p vout. The user must ensure that the vector
             *   dimensions are consistent
             */
            virtual void map(const FESystem::Numerics::VectorBase<ValType>& vin, FESystem::Numerics::VectorBase<ValType>& vout) const;
            
            /*!
             *   Calculates the Jacobian of the mapping at the specified point \p vin and returns it in the matrix \p mat
             */
            virtual void getMappingJacobian(const FESystem::Numerics::VectorBase<ValType>& vin, FESystem::Numerics::MatrixBase<ValType>& mat) const;
            
            /*!
             *   Returns the value of derivatives of the mapped function at \p val in X1, with respect to the derivative information provided
             *   in \p derivative_orders, whose i^th element is the order of the derivative with respect to the i^th coordinate. The derivative
             *   values are returned in \p deriv.
             */
            virtual void getFunctionDerivative(const std::vector<FESystemUInt>& derivative_orders, const FESystem::Numerics::VectorBase<ValType>& val,
                                               FESystem::Numerics::VectorBase<ValType>& deriv) const;
            
            /*!
             *   Returns the value of the discrete interpolation functions.
             */
            virtual void getDiscreteFunction(const FESystem::Numerics::VectorBase<ValType>& val, FESystem::Numerics::VectorBase<ValType>& deriv) const;
            
            /*!
             *   Returns the derivative of the discrete interpolation functions. The derivative order with respect to the i^th coordinate is specified in the
             *   i^th element of \p derivative_order. The location at which this is to be evaluated is \p val. The dimension of \p deriv is equal the number
             *   of interpolation function values.
             */
            virtual void getDiscreteFunctionDerivative(const std::vector<FESystemUInt>& derivative_orders, const FESystem::Numerics::VectorBase<ValType>& val,
                                                       FESystem::Numerics::VectorBase<ValType>& deriv) const;
            
        protected:
                        
            /*!
             *   stores the state of initialization of the object
             */
            FESystemBoolean if_initialized;
            
            /*!
             *    dimension of the mapped space
             */
            FESystemUInt dimension;

            /*!
             *    order of Legendre function for each dimension
             */
            FESystemUInt order;

            /*!
             *   These are the Legendre functions associated with each dimension.
             */
            const FESystem::Functions::LegendreFunction<ValType> * legendre_function;
            
            /*!
             *    This is a helper class that allows the calculation of order of the polynomial to be used for shape function tensor product. This works
             *    like an odometer such that each increment will increase the polynomial order of the last dimension, and if the last dimension has reached 
             *    the highest value, it will increment the second last and reset the last to the lowest value.
             */
            class OrderCounter
            {
            public:
                OrderCounter(const FESystemUInt dim, const FESystemUInt o);
                
                void increment();
                
                FESystemUInt getOrderForDim(const FESystemUInt dim) const;
                
                void write(std::ostream& out) const;
                
            protected:
                
                const FESystemUInt dimension;
                
                const FESystemUInt order;
                
                std::vector<FESystemUInt> order_per_dim;
            };
            
        };
    }
}


#endif /* defined(__FESystem__LegendreFunctionMapping__) */
