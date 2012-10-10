//
//  LagrangeFunctionMapping.h
//  FESystem
//
//  Created by Manav Bhatia on 4/13/12.
//  Copyright (c) 2012. All rights reserved.
//


#ifndef __fesystem_lagrange_function_mapping_h__
#define __fesystem_lagrange_function_mapping_h__


// FESystem includes
#include "Functions/FunctionMappingBase.h"


namespace FESystem
{
    // Forward declerations
    namespace Functions {template <typename ValType> class LagrangeFunction;}
    namespace Utility {template <typename ValType> class Table;}
    
    namespace Functions
    {
        /*!
         *   This class derives from the base class FunctionBase, and implements the Lagrange function 
         *   with abcissa and ordinate types as FESystemDoubles. 
         */
        template <typename ValType>
        class LagrangeFunctionMapping: public FESystem::Functions::DiscreteFunctionMappingBase<ValType>
        {
        public:
            /*!
             *    default constructor
             */
            LagrangeFunctionMapping();
            
            /*!
             *    destructor
             */
            virtual ~LagrangeFunctionMapping();
            
            /*!
             *   clears the dats structures
             */ 
            virtual void clear();
                        
            /*!
             *   This initializes the mapping using the given dimension, \p dim, for the mapped space, and the 
             *   Lagrange functions for the domain. One Lagrange function per domain is needed. 
             */ 
            void reinit(FESystemUInt dim, const std::vector<FESystem::Functions::LagrangeFunction<ValType>*>& funcs, const FESystem::Utility::Table<FESystemUInt>& pt_table);
                        
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
             *   This returns the ID of the point in the local axis. This is useful since the Lagrange function creates
             *   higher dimensional polynomials using tensor products, and each point has a unique number in each dimension 
             *   whose shape functions will be multiplied together to create the multidimensinal shape functions. 
             *   This must be set by the user through the reinit function. 
             */   
            FESystemUInt getDiscretePointNumberAlongLocalDimension(FESystemUInt p_id, FESystemUInt dim_id) const;
            
            /*!
             *   stores the state of initialization of the object
             */
            FESystemBoolean if_initialized;
            
            /*!
             *    dimension of the mapped space
             */
            FESystemUInt dimension;
            
            /*!
             *   These are the Lagrange functions associated with each dimension. 
             */
            std::vector<const FESystem::Functions::LagrangeFunction<ValType> *> lagrange_functions;
            
            /*!
             *   pointer to the table that stores the point to dimension data. This is a 2-D table with the number of rows equal to the number of 
             *   points and the number of columns equal to the number of dimensions.
             */
            const FESystem::Utility::Table<FESystemUInt>* point_id_table;
        };
        
    }
}



#endif  //__fesystem_lagrange_function_h__

