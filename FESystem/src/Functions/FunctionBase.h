//
//  FunctionBase.h
//  FESystem
//
//  Created by Manav Bhatia on 3/23/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem_function_base_h__
#define __fesystem_function_base_h__


// FESystem includes
#include "Base/FESystemTypes.h"



namespace FESystem
{
    namespace Numerics 
    {
        
        /*!
         *    this class defines a base class for one dimensional function types. The 
         *    class is defined as a template with two parameters: one as the abcissa type,
         *    and the other as an ordinate
         */
        template <typename AbscissaType, typename OrdinateType>
        class FunctionBase
        {
        public:
            
            /*!
             *    default constructor
             */
            FunctionBase();
            
            /*!
             *    destructor
             */
            virtual ~FunctionBase();
            
            /*!
             *   @returns value of the function at the given abcissa 
             *   @param abcissa value at which the function should be evaluated
             */
            virtual OrdinateType getFunctionValue(const AbscissaType& abcissa) const = 0;
            
            
            /*!
             *   @returns derivative of the function at the given abcissa 
             *   @param abcissa value at which the function should be evaluated 
             *   @param order derivative  
             */
            virtual OrdinateType getFunctionDerivative(const AbscissaType& abcissa, const FESystemUInt& order) const = 0;
                        
        protected:
            
        };
    }
}




template <typename AbcissaType, typename OrdinateType>
FESystem::Numerics::FunctionBase<AbcissaType, OrdinateType>::FunctionBase()
{

}




template <typename AbcissaType, typename OrdinateType>
FESystem::Numerics::FunctionBase<AbcissaType, OrdinateType>::~FunctionBase()
{

}





#endif // __fesystem_function_base_h__
