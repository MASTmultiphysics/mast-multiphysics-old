//
//  constant_field_function.h
//  MAST
//
//  Created by Manav Bhatia on 1/22/14.
//  Copyright (c) 2014 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_constant_field_function_h__
#define __MAST_constant_field_function_h__



// C++ includes
#include <vector>
#include <string>

// MAST includes
#include "Numerics/function_base.h"
#include "Numerics/constant_function.h"


namespace MAST {
    
    template <typename ValType>
    class ConstantFieldFunction: public MAST::FieldFunction<ValType> {
        
    public:
        ConstantFieldFunction(const std::string& nm,
                              const ValType& v):
        MAST::FieldFunction<ValType>(nm),
        _output(v),
        _func(NULL)
        { }

        
        ConstantFieldFunction(const std::string& nm, const MAST::ConstantFunction<ValType>& f):
        MAST::FieldFunction<ValType>(nm),
        _output(f()),
        _func(&f)
        { }

        
        virtual ~ConstantFieldFunction()
        { }
        
        
        /*!
         *    returns the function kind
         */
        virtual MAST::FunctionType type() const {
            return MAST::CONSTANT_FIELD_FUNCTION;
        }


        /*!
         *  returns true if the function depends on the provided value
         */
        virtual bool depends_on(Real* val) const {
            if (!_func)
                return false;
            else
                return _func->depends_on(val);
        }

        
        /*!
         *  returns true if the function depends on the provided value
         */
        virtual bool depends_on(const MAST::FunctionBase& f) const {
            if (!_func)
                return false;
            else if (_func == &f)
                return true;
            else
                return false;
        }

        
        /*!
         *    Returns the value of this function.
         */
        virtual void operator() (const Point& p, const Real t, ValType& v) const {
            v = _output;
        }
        
        
        
        /*!
         *    Returns the partial derivative of this function with respect to
         *    the sensitivity parameter \par p
         */
        virtual void partial_derivative (const MAST::SensitivityParameters& par,
                                         const Point& p, const Real t, ValType& v) const {
            return this->total_derivative(par, p, t, v);
        }
        
        
        
        /*!
         *    Returns the total derivative of this function with respect to
         *    the sensitivity parameter \par p
         */
        virtual void total_derivative (const MAST::SensitivityParameters& par,
                                       const Point& p, const Real t, ValType& v) const {
            // if the function was not specified, then the sensitivity is zero
            if (!_func) {
                v = 0.;
                return;
            }
            else {
                switch (par.total_order()) {
                    case 1:
                        if (&(par.get_first_order_derivative_parameter()) == _func) {
                            v = 1.;
                            return;
                        }
                        break;
                        
                    default:
                        libmesh_error();
                        break;
                }
            }
        }
        

        
        /*!
         *  sets the value of this function
         */
        void operator =(const ValType& val)
        {   _output = val; }

    protected:
        
        /*!
         *    constant output value
         */
        ValType _output;
        
        /*!
         *   function value that provides this
         */
        const MAST::ConstantFunction<ValType> *_func;
    };
}




#endif // __MAST_constant_field_function_h__
