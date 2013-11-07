//
//  temperature_function.h
//  MAST
//
//  Created by Manav Bhatia on 10/23/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_temperature_function_h__
#define __MAST_temperature_function_h__

// libMesh includes
#include "libmesh/point.h"

// MAST includes
#include "Numerics/function_base.h"



namespace MAST
{
    class Temperature: public MAST::FunctionValue<Real> {
    public:

        Temperature():
        MAST::FunctionValue<Real>("Temperature")
        { }
        
        /*!
         *    virtual destructor
         */
        virtual ~Temperature() { }
        
        /*!
         *    initialize the value at the specified point
         */
        virtual void initialize(const Point& p) = 0;
        
        /*!
         *    returns the temperature value
         */
        virtual Real operator() () const = 0;
        
        /*!
         *    returns the reference temperature
         */
        virtual Real reference() const = 0;
    };
    
    
    class ConstantTemperature: public MAST::Temperature {
    public:
        /*!
         *   constructor
         */
        ConstantTemperature():
        MAST::Temperature(),
        _temperature(0.),
        _reference(0.)
        { }

        /*!
         *    virtual destructor
         */
        virtual ~ConstantTemperature() { }

        /*!
         *    sets the value of temperature and the reference
         */
        void  set_temperature(Real temp, Real ref_temp) {
            _temperature = temp;
            _reference = ref_temp;
        }
        
        /*!
         *    initialize the value at the specified point
         */
        virtual void initialize(const Point& p) { }
        
        /*!
         *    returns the temperature value
         */
        virtual Real operator() () const {
            return _temperature;
        }
        
        /*!
         *    returns the reference temperature
         */
        virtual Real reference() const {
            return _reference;
        }
        
    protected:
        
        /*!
         *    temperature and reference value
         */
        Real _temperature, _reference;
    };
}








#endif  // __MAST_temperature_function_h__
