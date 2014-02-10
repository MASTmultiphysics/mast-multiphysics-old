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
    class ConstantFunction: public FieldFunction<ValType> {
    public:
        
        ConstantFunction(const std::string& nm,
                         const ValType& val):
        FieldFunction<ValType>(nm),
        _val(val)
        { }
        
        ConstantFunction(const MAST::ConstantFunction<ValType>& f):
        FieldFunction<ValType>(f),
        _val(f._val)
        { }

        
        /*!
         *   @returns a clone of the function
         */
        virtual std::auto_ptr<MAST::FieldFunction<ValType> > clone() const {
            return std::auto_ptr<MAST::FieldFunction<ValType> >
            (new MAST::ConstantFunction<ValType>(*this));
        }

        /*!
         *   returns the value of this function
         */
        virtual void operator() (const Point& p, const Real t, ValType& v) const
        {   v = _val; }
        
        
        /*!
         *   returns the sensitivity of this function
         */
        virtual void partial (const MAST::FieldFunctionBase& f,
                              const Point& p, const Real t,
                              ValType& v) const {
            this->total(f, p, t, v);
        }
        
        
        /*!
         *   returns the sensitivity of this function. This is the same as partial
         *   sensitivity for a constant function.
         */
        virtual void total (const MAST::FieldFunctionBase& f,
                            const Point& p, const Real t,
                            ValType& v) const {

            if (this == &f)
                v = 1.;
            else
                v = 0.;
        }
        
        
        
        /*!
         *    Returns the pointer to value of this function.
         */
        virtual ValType* ptr()
        {   return &_val; }
        
        
        
        /*!
         *  returns false since a constant function does not depend on any
         *  function.
         */
        virtual bool depends_on(const MAST::FieldFunctionBase& f) const {
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
