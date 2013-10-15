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
#include <memory>

//  MAST includes
#include "Base/MAST_data_types.h"

namespace MAST
{

    enum FunctionType {
        CONSTANT_FUNCTION,
        MULTILINEAR_FUNCTION_1D
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
         *   returns the value of this function
         */
        template <typename ValType>
        ValType operator() ();
        
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
    };
    
    
    template <typename ValType>
    class FunctionValue: public MAST::FunctionBase {
    public:
        FunctionValue(const std::string& nm):
        MAST::FunctionBase(nm)
        { }
        
        /*!
         *    Returns the value of this function.
         */
        virtual ValType operator() () = 0;
    };
    
    
    template <typename ValType>
    inline ValType
    FunctionBase::operator() () {
        return dynamic_cast<FunctionValue<ValType>&>(*this)();
    }
    
    
    
    /*!
     *    This is a function that does not change
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
        virtual ValType operator() ()
        { return _val; }
        
        /*!
         *  sets the value of this function
         */
        void operator =(ValType& val)
        { _val = val;}
        
    protected:
        
        ValType _val;
    };
 
    
    /*!
     *    creates a function and returns it as a smart pointer
     */
    template <typename ValType>
    std::auto_ptr<MAST::FunctionBase> build_function(const std::string& nm,
                                                     MAST::FunctionType type) {
        
        std::auto_ptr<MAST::FunctionBase> rval;
        
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
