//
//  function_base.h
//  MAST
//
//  Created by Manav Bhatia on 10/15/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_function_base_h__
#define __MAST_function_base_h__

// C++ includes
#include <string>
#include <map>
#include <set>
#include <memory>


//  MAST includes
#include "Base/MAST_data_types.h"
#include "Numerics/sensitivity_parameters.h"

// libMesh includes
#include "libmesh/point.h"

namespace MAST
{
    
    
    enum FunctionType {
        CONSTANT_FUNCTION,
        MULTILINEAR_FUNCTION_1D
    };
    
    
    enum FunctionAttributeType {
        SHAPE_FUNCTION,
        MATERIAL_FUNCTION,
        ELEMENT_SIZE_FUNCTION
    };
    
    class FunctionBase {
    public:
        
        /*!
         *   initializes the parameter to the given name
         */
        FunctionBase(const std::string& nm ):
        _name(nm)
        { }
        
        /*!
         *   virtual destructor
         */
        virtual ~FunctionBase() { }
        
        /*!
         *    returns the function kind
         */
        virtual MAST::FunctionType type() const = 0;
        
        /*!
         *   returns the name of this function
         */
        const std::string& name() const {
            return _name;
        }
        
        /*!
         *  returns true if the function depends on the provided value
         */
        virtual bool depends_on(Real* val) const = 0;

        /*!
         *  returns true if the function depends on the provided value
         */
        virtual bool depends_on(const MAST::FunctionBase& f) const = 0;

        /*!
         *  returns true if the function depends on the provided set 
         *  of parameters
         */
        virtual bool depends_on(const MAST::SensitivityParameters& p) const;

        
        /*!
         *  checks if the function has an attribute
         */
        bool has_attribute(MAST::FunctionAttributeType a) const {
            return _attributes.count(a);
        }
        
        /*!
         *   returns the value of this function
         */
        template <typename ValType>
        ValType operator() () const;

        
        /*!
         *   returns the partial derivative of this function with respect to the
         *   parameters specified in the argument
         */
        template <typename ValType>
        ValType partial_derivative (const MAST::SensitivityParameters& p) const;


        /*!
         *   returns the total derivative of this function with respect to the
         *   parameters specified in the argument
         */
        template <typename ValType>
        ValType total_derivative (const MAST::SensitivityParameters& p) const;

        
        /*!
         *   returns the value of this function at Point \par p and time \par t.
         *   This is used for FieldFunction objects. The value is returned in 
         *   object \par v.
         */
        template <typename ValType>
        void operator() (const Point& p, const Real t, ValType& v) const;
        
        
        /*!
         *   returns the partial derivative of this function with respect to the
         *   sensitivity parameter \par par at Point \par p and time \par t.
         *   This is used for FieldFunction objects. The value is returned in
         *   object \par v.
         */
        template <typename ValType>
        void partial_derivative (const MAST::SensitivityParameters& par,
                                 const Point& p, const Real t,
                                 ValType& v) const;

        /*!
         *   returns the total derivative of this function with respect to the
         *   sensitivity parameter \par par at Point \par p and time \par t.
         *   This is used for FieldFunction objects. The value is returned in
         *   object \par v.
         */
        template <typename ValType>
        void total_derivative (const MAST::SensitivityParameters& par,
                               const Point& p, const Real t,
                               ValType& v) const;

        /*!
         *    adds a parameter on which this function is dependent
         */
        void add_function(const FunctionBase& f){
            // make sure it does not already exist
            libmesh_assert(!_function_parameters.count(f.name()));
            _function_parameters.insert(std::pair<std::string, const FunctionBase*>
                                        (f.name(), &f));
        }
        
    protected:
        
        /*!
         *    name of this parameter
         */
        std::string _name;
        
        /*!
         *   maps of functions that \p this function depends on
         */
        std::map<std::string, const FunctionBase*> _function_parameters;
        
        /*!
         *   attributes of this function
         */
        std::set<MAST::FunctionAttributeType> _attributes;
    };
    

    /*!
     *    This creates a base class for functions that depend on specific 
     *    parameters and provide sensitivity operations with respect to the 
     *    parameters.
     */
    template <typename ValType>
    class FunctionValue: public MAST::FunctionBase {
    public:
        FunctionValue(const std::string& nm):
        MAST::FunctionBase(nm)
        { }
        
        /*!
         *    Returns the value of this function.
         */
        virtual ValType operator() () const = 0;
        
        /*!
         *    Returns the partial derivative of this function with respect to
         *    the sensitivity parameter \par p
         */
        virtual ValType partial_derivative (const MAST::SensitivityParameters& p) const = 0;

        /*!
         *    Returns the total derivative of this function with respect to
         *    the sensitivity parameter \par p
         */
        virtual ValType total_derivative (const MAST::SensitivityParameters& p) const = 0;

        /*!
         *    Returns the pointer to value of this function.
         */
        virtual ValType* ptr() = 0;

        /*!
         *  sets the value of this function
         */
        virtual void operator =(const ValType& val) = 0;
    };

    

    /*!
     *    This creates the base class for functions that have a saptial and 
     *    temporal dependence, and provide sensitivity operations with respect 
     *    to the functions and parameters.
     */
    template <typename ValType>
    class FieldFunction: public MAST::FunctionBase {
    public:
        FieldFunction(const std::string& nm):
        MAST::FunctionBase(nm)
        { }
        
        /*!
         *    Returns the value of this function.
         */
        virtual void operator() (const Point& p, const Real t, ValType& v) const = 0;
        
        /*!
         *    Returns the partial derivative of this function with respect to
         *    the sensitivity parameter \par p
         */
        virtual void partial_derivative (const MAST::SensitivityParameters& par,
                                         const Point& p, const Real t, ValType& v) const = 0;
        
        /*!
         *    Returns the total derivative of this function with respect to
         *    the sensitivity parameter \par p
         */
        virtual void total_derivative (const MAST::SensitivityParameters& par,
                                       const Point& p, const Real t, ValType& v) const = 0;
    };

    
    
    template <typename ValType>
    ValType
    MAST::FunctionBase::operator() () const {
        return dynamic_cast<FunctionValue<ValType>&>(*this)();
    }

    

    template <typename ValType>
    void
    MAST::FunctionBase::operator() (const Point& p,
                                    const Real t,
                                    ValType& v) const {
        return dynamic_cast<FieldFunction<ValType>&>(*this)();
    }

    
    
    
    template <typename ValType>
    ValType
    MAST::FunctionBase::partial_derivative (const MAST::SensitivityParameters& par) const {
        return dynamic_cast<MAST::FunctionValue<ValType>&>(*this).partial_derivative(par);
    }

    
    template <typename ValType>
    ValType
    MAST::FunctionBase::total_derivative (const MAST::SensitivityParameters& par) const {
        return dynamic_cast<MAST::FunctionValue<ValType>&>(*this).total_derivative(par);
    }

    
    template <typename ValType>
    void
    MAST::FunctionBase::partial_derivative (const MAST::SensitivityParameters& par,
                                            const Point& p,
                                            const Real t,
                                            ValType& v) const {
        dynamic_cast<MAST::FieldFunction<ValType>&>(*this).
        partial_derivative(par, p, t, v);
    }

    
    template <typename ValType>
    void
    MAST::FunctionBase::total_derivative (const MAST::SensitivityParameters& par,
                                          const Point& p,
                                          const Real t,
                                          ValType& v) const {
        dynamic_cast<MAST::FieldFunction<ValType>&>(*this).
        total_derivative(par, p, t, v);
    }

    
    /*!
     *    This is a function that does not change.
     */
    template <typename ValType>
    class ConstantFunction: public FunctionValue<ValType> {
    public:
        
        ConstantFunction(const std::string& nm):
        FunctionValue<ValType>(nm),
        _val(ValType())
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
        virtual ValType operator() () const {
            
            return _val; }


        /*!
         *   returns the sensitivity of this function
         */
        virtual ValType partial_derivative (const MAST::SensitivityParameters& p) const {
            
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
         *   returns the sensitivity of this function. This is the same as partial 
         *   sensitivity for a constant function.
         */
        virtual ValType total_derivative (const MAST::SensitivityParameters& p) const {
            return this->partial_derivative(p);
        }

        
        /*!
         *    Returns the pointer to value of this function.
         */
        virtual ValType* ptr() {
            return &_val;
        };

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
        void operator =(const ValType& val) {
            
            _val = val; }
        
    protected:
        
        ValType _val;
    };
 
    
    /*!
     *    creates a function and returns it as a smart pointer
     */
    template <typename ValType>
    std::auto_ptr<MAST::FunctionValue<ValType> >
    build_function(const std::string& nm, MAST::FunctionType type) {
        
        std::auto_ptr<MAST::FunctionValue<ValType> > rval;
        
        switch (type) {
            case MAST::CONSTANT_FUNCTION:
                rval.reset(new MAST::ConstantFunction<ValType>(nm));
                break;
                
            default:
                // should not get here
                libmesh_assert(false);
                break;
        }
        
        return rval;
    }
    
}

#endif // __MAST_function_base_h__
