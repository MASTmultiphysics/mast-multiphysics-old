//
//  LegendreFunction.h
//  FESystem
//
//  Created by Manav Bhatia on 1/25/13.
//
//

#ifndef FESystem_LegendreFunction_h
#define FESystem_LegendreFunction_h

// FESystem includes
#include "Functions/FunctionBase.h"
#include "Base/FESystemExceptions.h"


namespace FESystem
{
    namespace Functions
    {
        /*!
         *   This class derives from the base class FunctionBase, and implements the Legendre function
         *   with abcissa and ordinate types as ValTypes.
         */
        template <typename ValType>
        class LegendreFunction: public FESystem::Numerics::FunctionBase<ValType, ValType>
        {
        public:
            /*!
             *    default constructor
             */
            LegendreFunction();
            
            /*!
             *    destructor
             */
            virtual ~LegendreFunction();
            
            
            /*!
             *   dummy function for LegendreFunction since the evaluation needs the point for which the
             *   function needs to be evaluated. Use the overloaded function instead
             */
            virtual ValType getFunctionValue(const ValType& abcissa) const;
            
            
            /*!
             *   @returns value of the function at the given abcissa
             *   @param abcissa value at which the function should be evaluated
             *   @param func_order specified the order of the function value to be returned
             */
            ValType getFunctionValue(const ValType& abcissa, const FESystemUInt& func_order) const;
            
            
            /*!
             *   dummy function for LegendreFunction since the evaluation needs the point for which the
             *   function needs to be evaluated. Use the overloaded function instead
             */
            virtual ValType getFunctionDerivative(const ValType& abcissa, const FESystemUInt& deriv_order) const;
            
            
            /*!
             *   @returns value of the function at the given abcissa
             *   @param abcissa value at which the function should be evaluated
             *   @param func_order specified the order of the function value to be returned
             *   @param deriv_order specified the derivative order to be returned
             */
            ValType getFunctionDerivative(const ValType& abcissa, const FESystemUInt& func_order, const FESystemUInt& deriv_order) const;
            
        protected:
            
        };
    }
}



template <typename ValType>
inline
FESystem::Functions::LegendreFunction<ValType>::LegendreFunction():
FESystem::Numerics::FunctionBase<ValType, ValType>()
{
}


template <typename ValType>
inline
FESystem::Functions::LegendreFunction<ValType>::~LegendreFunction()
{
}



template <typename ValType>
inline
ValType
FESystem::Functions::LegendreFunction<ValType>::getFunctionValue(const ValType& abcissa) const
{
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}


template <typename ValType>
inline
ValType
FESystem::Functions::LegendreFunction<ValType>::getFunctionValue(const ValType& abcissa, const FESystemUInt& func_order) const
{
    FESystemAssert0((abcissa >= -1.0) && (abcissa <= 1.0), FESystem::Exception::InvalidValue);
    
    switch (func_order)
    {
        case 0:
            return 1;
            break;
            
        case 1:
            return abcissa;
            break;

        default:
            return ((2*func_order-1)*abcissa*this->getFunctionValue(abcissa, func_order-1) - (func_order-1)*this->getFunctionValue(abcissa, func_order-2))/func_order;
            break;
    }
}


template <typename ValType>
inline
ValType
FESystem::Functions::LegendreFunction<ValType>::getFunctionDerivative(const ValType& abcissa, const FESystemUInt& order) const
{
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}



template <typename ValType>
inline
ValType
FESystem::Functions::LegendreFunction<ValType>::getFunctionDerivative(const ValType& abcissa, const FESystemUInt& func_order, const FESystemUInt& deriv_order) const
{
    FESystemAssert0((abcissa >= -1.0) && (abcissa <= 1.0), FESystem::Exception::InvalidValue);
    FESystemAssert0(deriv_order > 0, FESystem::Exception::InvalidValue);
    
    switch (func_order)
    {
        case 0:
            return 0.0;
            break;
            
        case 1:
        {
            switch (deriv_order)
            {
                case 1:
                    return 1.0;
                    break;
                    
                default:
                    return 0.0;
                    break;
            }
        }
            break;
            
        default:
        {
            ValType val = (2*func_order-1)*abcissa*this->getFunctionDerivative(abcissa, func_order-1, deriv_order) - (func_order-1)*this->getFunctionDerivative(abcissa, func_order-2, deriv_order);
            if (deriv_order == 1)
                val += (2*func_order-1)*this->getFunctionValue(abcissa, func_order-1);
            else // deriv_order > 1
                val += (2*func_order-1)*this->getFunctionDerivative(abcissa, func_order-1, deriv_order-1);
            return val/(1.0*func_order);
        }
            break;
    }
}



#endif // FESystem_LegendreFunction_h
