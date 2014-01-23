//
//  sensitivity_parameters.h
//  MAST
//
//  Created by Manav Bhatia on 11/8/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_sensitivity_parameters_h__
#define __MAST_sensitivity_parameters_h__

// C++ includes
#include <map>

// MAST includes
#include "Base/MAST_data_types.h"

namespace MAST {
    
    // Forward declerations
    class FunctionBase;
    
    /*!
     *    class provides an interface to store and query parameter set for
     *    sensitivity calculations.
     */
    class SensitivityParameters {
    public:
       
        /*!
         *   data type of the map used to store the function and its derivative
         *   order
         */
        typedef std::map<const MAST::FunctionBase*, unsigned int> ParameterMap;

        /*!
         *    public constructor
         */
        SensitivityParameters():
        _order(0)
        { }
        
        /*!
         *    virtual destructor
         */
        virtual ~SensitivityParameters()
        { }

        /*!
         *    adds parameter and its derivative order. \p o should be greater 
         *    than 0.
         */
        void add_parameter(const MAST::FunctionBase* f, unsigned int o) {
            libmesh_assert(!this->contains(f));
            libmesh_assert_greater(o, 0);
            _parameters[f] = o;
            _order += o;
        }
        
        /*!
         *    checks if the specified parameter is included for sensitivity
         */
        bool contains(const MAST::FunctionBase* f) const {
            return (_parameters.find(f) != _parameters.end());
        }
        
        /*!
         *    returns the derivative order wrt to the specified parameter
         */
        unsigned int parameter_order(const MAST::FunctionBase* f) const;
        
        /*!
         *    returns the total derivative order for this set
         */
        unsigned int total_order() const {
            return _order;
        }

        /*!
         *    returns the parameter map
         */
        bool shape_sensitivity() const;
        
        /*!
         *  number of parameters
         */
        unsigned int n_params() const {
            return (unsigned int)this->_parameters.size();
        }
        
        /*!
         *    returns the parameter map
         */
        const SensitivityParameters::ParameterMap& get_map() const {
            return this->_parameters;
        }
        
        
        /*!
         *    this is a utility function to return the first order sensitivity 
         *    parameter. This should only be used when total_order() == 1
         */
        const MAST::FunctionBase& get_first_order_derivative_parameter() const {
            libmesh_assert_equal_to(this->total_order(), 1);
            return *_parameters.begin()->first;
        }
        
    protected:

        /*!
         *    total derivative order of this parameter set
         */
        unsigned int _order;
        
        /*!
         *    sensitivity parameter set
         */
        SensitivityParameters::ParameterMap _parameters;
    };
}


#endif // __MAST_sensitivity_parameters_h__
