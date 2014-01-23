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

// libMesh includes
#include "libmesh/const_function.h"


namespace MAST {
    
    template <typename ValType>
    class ConstantFieldFunction: public MAST::FieldFunction<ValType> {
        
    public:
        ConstantFieldFunction(const std::string& nm, const ValType& o):
        MAST::FieldFunction<ValType>(nm),
        _function(NULL),
        _function_total_derivative(NULL) {
            
            _function = new ConstFunction(o);
        }
        
        
        virtual ~ConstantFieldFunction() {
            
            // delete the objects before exiting
            if (_function)
                delete _function;
            
            if (_function_total_derivative)
                delete _function_total_derivative;
        }
        
        
        
        /*!
         *    initialize the data structure for sensitivity analysis
         */
        void init_for_sens(const ValType& o) {
            // prepare the localized solution
            if (_function_total_derivative)
                delete _function_total_derivative;
            
            _function_total_derivative = new ConstFunction(o);
        }
        
        
        /*!
         *    Returns the value of this function.
         */
        virtual void operator() (const Point& p, const Real t, ValType& v) const {
            (*_function)(p, t, v);
        }
        
        
        
        /*!
         *    Returns the partial derivative of this function with respect to
         *    the sensitivity parameter \par p
         */
        virtual void partial_derivative (const MAST::SensitivityParameters& par,
                                         const Point& p, const Real t, ValType& v) const {
            // there are no partial derivatives here. A system only calculates
            // total derivatives
            libMesh::err
            << "Error: ConstantFieldFunction does not provide partial derivative. "
            << std::endl << "Exiting." << std::endl;
            libmesh_error(false);
        }
        
        
        
        /*!
         *    Returns the total derivative of this function with respect to
         *    the sensitivity parameter \par p
         */
        virtual void total_derivative (const MAST::SensitivityParameters& par,
                                       const Point& p, const Real t, ValType& v) const {
            // make sure that this has been initialized
            libmesh_assert(_function_total_derivative);
            (*_function_total_derivative)(p, t, v);
        }
        
        
    protected:
        
        /*!
         *   mesh function object that provides the interpolation
         */
        ConstFunction *_function;
        
        /*!
         *   total derivative mesh function
         */
        ConstFunction *_function_total_derivative;
    };
}




#endif // __MAST_constant_field_function_h__
