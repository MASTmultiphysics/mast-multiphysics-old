//
//  constant_function.h
//  MAST
//
//  Created by Manav Bhatia on 1/23/14.
//  Copyright (c) 2014 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_constant_function_h__
#define __MAST_constant_function_h__


// MAST includes
#include "Numerics/function_base.h"


namespace MAST {
    
    /*!
     *    This is a function that does not change.
     */
    template <typename ValType>
    class ConstantFunction: public FunctionValue<ValType> {
    public:
        
        ConstantFunction(const std::string& nm,
                         const ValType& val):
        FunctionValue<ValType>(nm),
        _val(val)
        { }
        
        /*!
         *   returns the type of this function: CONSTANT_FUNCTION
         */
        virtual MAST::FunctionType type() const {
            
            return MAST::CONSTANT_FUNCTION;
        }
        
        /*!
         *   returns the value of this function
         */
        virtual ValType operator() () const
        {   return _val; }
        
        
        /*!
         *   returns the sensitivity of this function
         */
        virtual ValType partial_derivative (const MAST::SensitivityParameters& p) const {
            this->total_derivative(p);
        }
        
        
        /*!
         *   returns the sensitivity of this function. This is the same as partial
         *   sensitivity for a constant function.
         */
        virtual ValType total_derivative (const MAST::SensitivityParameters& p) const {
            // only first order sensitivities are calculated at this point
            libmesh_assert_equal_to(p.total_order(), 1);
            const MAST::FunctionBase& f = p.get_first_order_derivative_parameter();
            
            if (&f == this)
                return 1.;
            else
                return 0.;
        }
        
        
        
        /*!
         *    Returns the pointer to value of this function.
         */
        virtual ValType* ptr()
        {   return &_val; }
        
        
        
        /*!
         *  returns true if the function depends on the provided value
         */
        virtual bool depends_on(Real* val) const {
            return (val == &_val);
        }
        
        
        
        /*!
         *  returns false since a constant function does not depend on any
         *  function.
         */
        virtual bool depends_on(const FunctionBase& f) const {
            if (&f == this)
                return true;
            else
                return false;
        }
        
        
        
        /*!
         *  sets the value of this function
         */
        void operator =(const ValType& val)
        {   _val = val; }
        
        
    protected:
        
        ValType _val;
    };
    
}

#endif // __MAST_constant_function_h__
