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
    class Temperature: public MAST::FieldFunction<Real> {
    public:

        Temperature():
        MAST::FieldFunction<Real>("Temperature")
        { }
        
        Temperature(const MAST::Temperature& f):
        MAST::FieldFunction<Real>(f)
        { }
        
        /*!
         *    virtual destructor
         */
        virtual ~Temperature() { }
        
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

        ConstantTemperature(const MAST::ConstantTemperature& f):
        MAST::Temperature(f),
        _temperature(f._temperature),
        _reference(f._reference)
        { }
        
        /*!
         *   @returns a clone of the function
         */
        virtual std::auto_ptr<MAST::FieldFunction<Real> > clone() const {
            return std::auto_ptr<MAST::FieldFunction<Real> >
            (new MAST::ConstantTemperature(*this));
        }

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
         *    Returns the value of this function.
         */
        virtual void operator() (const Point& p, const Real t, Real& v) const {
            v = _temperature;
        }
        
        /*!
         *    Returns the partial derivative of this function with respect to
         *    the sensitivity parameter \par p
         */
        virtual void partial (const MAST::FieldFunctionBase& f,
                              const Point& p, const Real t, Real& v) const {
            libmesh_error();
        }
        
        /*!
         *    Returns the total derivative of this function with respect to
         *    the sensitivity parameter \par p
         */
        virtual void total (const MAST::FieldFunctionBase& f,
                                       const Point& p, const Real t, Real& v) const {
            libmesh_error();
        }
        
        /*!
         *    returns the reference temperature
         */
        virtual Real reference() const {
            return _reference;
        }
        
        /*!
         *  returns true if the function depends on the provided value
         */
        virtual bool depends_on(Real* val) const {
            return (val == &_temperature);
        }
        
    protected:
        
        /*!
         *    temperature and reference value
         */
        Real _temperature, _reference;
    };
}








#endif  // __MAST_temperature_function_h__
