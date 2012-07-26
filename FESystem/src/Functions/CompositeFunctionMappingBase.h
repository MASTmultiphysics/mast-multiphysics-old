//
//  CompositeFunctionMappingBase.h
//  FESystem
//
//  Created by Manav Bhatia on 4/14/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem_composite_function_mapping_base_h__
#define __fesystem_composite_function_mapping_base_h__

// FESystem includes
#include "Functions/FunctionMappingBase.h"

namespace FESystem
{
    namespace Functions
    {
        /*!
         *   This method deals with three mapped spaces: X1, which is a m1-dimensional space, X2, which is 
         *   a m2-dimensional space and X3, which is a m3-dimensional space, such that X2 = X2(X1), and 
         *   X3 = X3(X1) are known mappings. This provides the methods to calculate the mappings for X3 = X3(X2)
         *   inferred from the X2(X1) and X3(X1) mappings. Note that the method will provide the composite 
         *   Jacobians  dX3/dX2 and derivatives of a specified order. X2 here is called the sub-function. 
         *   This is useful for finite element derivative approximations in the physical coordinates (x, y, z) (the sub-function, X2, in this case). 
         *   For (xi, eta, phi) (X1 in this case) being the computational coordinates, and w(xi, eta, phi) (X3 in this case) as the original available 
         *   mapping, the dw/dx, dw/dy and dw/dz (and other higher derivative orders) can be calculated by this method. 
         */
        template <typename ValType>
        class CompositeFunctionMappingBase: FESystem::Functions::FunctionMappingBase<ValType>
        {
        public:
            /*!
             *   Constructor
             */
            CompositeFunctionMappingBase();
            
            virtual ~CompositeFunctionMappingBase();
            
            /*!
             *   This is not a valid function call for this class. The user may call the functions for the original and sub-function mappings
             */
            virtual std::pair<FESystemUInt, FESystemUInt> getDimensions() const; 

            /*!
             *   clears the data structures
             */
            void clear();
            
            /*!
             *   initializes the data structure for the two mappings. The X3(X1) is defined by \p original and X2(X1) is defined by 
             *   \p sub.  (Please see the class description for discussion of X1, X2 and X3). 
             */
            void reinit(const FESystem::Functions::FunctionMappingBase<ValType>& original, const FESystem::Functions::FunctionMappingBase<ValType>& sub);
            
            /*!
             *   Returns a constant reference to the X3(X1) mapping. 
             */
            const FESystem::Functions::FunctionMappingBase<ValType>& getOriginalFunctionMapping() const;

            /*!
             *   Returns a constant reference to the X2(X1) mapping.  
             */
            const FESystem::Functions::FunctionMappingBase<ValType>& getSubFunctionMapping() const;
            
            /*!
             *   This method is not defined for the composite mapping, and should not be used. 
             */
            virtual void map(const FESystem::Numerics::VectorBase<ValType>& vin, FESystem::Numerics::VectorBase<ValType>& vout) const;
            
            //            /*!
            //             *   Inverse maps the input vector \p vin to the original space vector \p vout. The user must ensure that the vector 
            //             *   dimensions are consistent
            //             */
            //            virtual void inverseMap(const FESystem::Numerics::VectorBase<ValType>& vin, FESystem::Numerics::VectorBase<ValType>& vout)=0;
            
            /*!
             *   Calculates and returns the d X3/ d X2 Jacobian in \p mat at the point \p vin defined in the X1 space. 
             */
            virtual void getMappingJacobian(const FESystem::Numerics::VectorBase<ValType>& vin, FESystem::Numerics::MatrixBase<ValType>& mat) const;
            
            /*!
             *   Calculates and returns the derivative of X3 in the X2 space at point \p val, where \p derivative_orders defines the derivative order 
             *   of each coordinate in X2 space.
             */
            virtual void getFunctionDerivative(const std::vector<FESystemUInt>& derivative_orders, const FESystem::Numerics::VectorBase<ValType>& val, FESystem::Numerics::VectorBase<ValType>& deriv) const;
            
            
        protected:
            
            /*!
             *   stores the state of initialization
             */
            FESystemBoolean if_initialized;

            /*!
             *   pointer to the X3(X1) mapping
             */
            const FESystem::Functions::FunctionMappingBase<ValType>* original_function_mapping;

            /*!
             *   pointer to the X2(X1) mapping
             */
            const FESystem::Functions::FunctionMappingBase<ValType>* sub_function_mapping;
        };
        
        
        /*!
         *   This defines a specialized class for discrete composite mapping functions where a separate interpolation functions exist for degrees of freedoms,
         *   such as the Lagrange interpolation function. The mapping function provides an additional method to return these interpolation functions
         *   in a vector form (as opposed to proving the interpolated quantity). 
         */
        template <typename ValType>
        class DiscreteCompositeMappingFunctionBase: public FESystem::Functions::CompositeFunctionMappingBase<ValType>
        {
        public: 
            DiscreteCompositeMappingFunctionBase();
            
            ~DiscreteCompositeMappingFunctionBase();
            
            /*!
             *   clears the data structures
             */
            void clear();
            
            /*!
             *   This is an invalid function call and an exception will be thrown. 
             */
            virtual void inverseMap(const FESystem::Numerics::VectorBase<ValType>& vin, FESystem::Numerics::VectorBase<ValType>& vout) const;

            /*!
             *   initializes the data structure for the two mappings. The X3(X1) is defined by \p original and X2(X1) is defined by 
             *   \p sub.  (Please see the class description for discussion of X1, X2 and X3). 
             */
            void reinit(const FESystem::Functions::DiscreteFunctionMappingBase<ValType>& original, const FESystem::Functions::DiscreteFunctionMappingBase<ValType>& sub);
            
            /*!
             *   Returns a constant reference to the X3(X1) mapping. 
             */
            const FESystem::Functions::DiscreteFunctionMappingBase<ValType>& getOriginalFunctionMapping() const;
            
            /*!
             *   Returns a constant reference to the X2(X1) mapping.  
             */
            const FESystem::Functions::DiscreteFunctionMappingBase<ValType>& getSubFunctionMapping() const;

            /*!
             *   Returns the discrete interpolation functions. The location at which this is to be evaluated is \p val. 
             */
            virtual void getDiscreteFunction(const FESystem::Numerics::VectorBase<ValType>& vin, FESystem::Numerics::VectorBase<ValType>& vout) const;

            /*!
             *   Returns the derivative of the discrete interpolation functions. The derivative order with respect to the i^th coordinate is specified in the 
             *   i^th element of \p derivative_order. The location at which this is to be evaluated is \p val. The dimension of \p deriv is equal the number 
             *   of interpolation function values. 
             */
            virtual void getDiscreteFunctionDerivative(const std::vector<FESystemUInt>& derivative_orders, const FESystem::Numerics::VectorBase<ValType>& vin, 
                                                       FESystem::Numerics::VectorBase<ValType>& deriv) const;
            
            
        protected:
            
            /*!
             *   pointer to the X3(X1) mapping
             */
            const FESystem::Functions::DiscreteFunctionMappingBase<ValType>* original_discrete_function_mapping;
            
            /*!
             *   pointer to the X2(X1) mapping
             */
            const FESystem::Functions::DiscreteFunctionMappingBase<ValType>* sub_discrete_function_mapping;
            
        };

    }
}


#endif //__fesystem_composite_function_mapping_base_h__

