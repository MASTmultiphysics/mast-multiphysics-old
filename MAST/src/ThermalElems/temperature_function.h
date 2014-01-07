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
        
        /*!
         *    returns the function kind
         */
        virtual MAST::FunctionType type() const {
            return MAST::CONSTANT_FUNCTION;
        }
        
        /*!
         *  returns true if the function depends on the provided value
         */
        virtual bool depends_on(Real* val) const {
            return (val == &_temperature);
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
         *   returns the sensitivity of this function
         */
        virtual Real sensitivity (const MAST::SensitivityParameters& p) const
        {
            // only first order sensitivities are calculated at this point
            libmesh_assert_equal_to(p.total_order(), 1);
            
            const MAST::SensitivityParameters::ParameterMap& p_map = p.get_map();
            MAST::SensitivityParameters::ParameterMap::const_iterator it, end;
            it = p_map.begin(); end = p_map.end();
            
            const MAST::FunctionBase& f = *(it->first);
            
            if (this->depends_on(f))
                return 1.;
            else
                return 0.;
        }
        
        /*!
         *    Returns the pointer to value of this function.
         */
        virtual Real* ptr() {
            return &_temperature;
        };
        
        /*!
         *  sets the value of this function
         */
        void operator =(const Real& val)
        { _temperature = val;}

    protected:
        
        /*!
         *    temperature and reference value
         */
        Real _temperature, _reference;
    };
}








#endif  // __MAST_temperature_function_h__
