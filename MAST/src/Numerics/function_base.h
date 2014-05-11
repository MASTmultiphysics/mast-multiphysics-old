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


// libMesh includes
#include "libmesh/point.h"

namespace MAST
{
    
    
    class FieldFunctionBase {
    public:
        
        /*!
         *   initializes the parameter to the given name
         */
        FieldFunctionBase(const std::string& nm ):
        _name(nm),
        _master(NULL)
        { }

        FieldFunctionBase(const MAST::FieldFunctionBase& f):
        _name(f._name),
        _master(f.master())
        { }

        
        /*!
         *   virtual destructor
         */
        virtual ~FieldFunctionBase() { }
        
        
        /*!
         *   returns the name of this function
         */
        const std::string& name() const {
            return _name;
        }
        
        
        const MAST::FieldFunctionBase* master() const {
            if (_master)     // this function has a master
                return _master;
            else             // this is the master
                return this;
        }
        
        
        /*!
         *  returns true if the function depends on the provided value
         */
        virtual bool depends_on(const MAST::FieldFunctionBase& f) const {
            if ((_functions.count(f.master())) || // one of the functions is the master
                (f.master() == this->master()))   // this function is the same
                return true;
            
            // check with all functions if they are dependent
            std::set<const MAST::FieldFunctionBase*>::const_iterator
            it = _functions.begin(), end = _functions.end();
            
            for ( ; it != end; it++)
                if ((*it)->depends_on(*f.master()))
                    return true;
            
            // if it gets here, then there is no dependency
            return false;
        }

        
        /*!
         *  @returns true if the function is a shape parameter. False by
         *  default. This should be reimplemneted in a new function
         *  that is a shape function.
         */
        virtual bool is_shape_parameter() const {
            return false;
        }
        
        /*!
         *   returns the value of this function at libMesh::Point \par p and time \par t.
         *   This is used for FieldFunction objects. The value is returned in 
         *   object \par v.
         */
        template <typename ValType>
        void operator() (const libMesh::Point& p, const libMesh::Real t, ValType& v) const;
        
        
        /*!
         *   returns the partial derivative of this function with respect to the
         *   sensitivity parameter \par par at libMesh::Point \par p and time \par t.
         *   This is used for FieldFunction objects. The value is returned in
         *   object \par v.
         */
        template <typename ValType>
        void partial (const MAST::FieldFunctionBase& f,
                      const libMesh::Point& p, const libMesh::Real t,
                      ValType& v) const;
        
        /*!
         *   returns the total derivative of this function with respect to the
         *   sensitivity parameter \par par at libMesh::Point \par p and time \par t.
         *   This is used for FieldFunction objects. The value is returned in
         *   object \par v.
         */
        template <typename ValType>
        void total (const MAST::FieldFunctionBase& f,
                    const libMesh::Point& p, const libMesh::Real t,
                    ValType& v) const;
        
    protected:
        
        /*!
         *    name of this parameter
         */
        std::string _name;
        
        /*!
         *    pointer to the master function that this is a copy of. A NULL 
         *    pointer implies that this is the master function
         */
        const MAST::FieldFunctionBase* _master;
        
        /*!
         *   set of functions that \p this function depends on
         */
        std::set<const MAST::FieldFunctionBase*> _functions;
    };
    
    
    
    /*!
     *    This creates the base class for functions that have a saptial and 
     *    temporal dependence, and provide sensitivity operations with respect 
     *    to the functions and parameters.
     */
    template <typename ValType>
    class FieldFunction: public MAST::FieldFunctionBase {
    public:
        FieldFunction(const std::string& nm):
        MAST::FieldFunctionBase(nm)
        { }

        FieldFunction(const MAST::FieldFunction<ValType>& f):
        MAST::FieldFunctionBase(f)
        { }

        /*!
         *   @returns a clone of the function
         */
        virtual std::auto_ptr<MAST::FieldFunction<ValType> > clone() const = 0;

        /*!
         *    Returns the value of this function.
         */
        virtual void operator() (const libMesh::Point& p, const libMesh::Real t, ValType& v) const = 0;
        
        /*!
         *    Returns the partial derivative of this function with respect to
         *    the sensitivity parameter \par p
         */
        virtual void partial (const MAST::FieldFunctionBase& f,
                              const libMesh::Point& p, const libMesh::Real t, ValType& v) const = 0;
        
        /*!
         *    Returns the total derivative of this function with respect to
         *    the sensitivity parameter \par p
         */
        virtual void total (const MAST::FieldFunctionBase& f,
                            const libMesh::Point& p, const libMesh::Real t, ValType& v) const = 0;
    };

    
    
    template <typename ValType>
    void
    MAST::FieldFunctionBase::operator() (const libMesh::Point& p,
                                         const libMesh::Real t,
                                         ValType& v) const {
        return dynamic_cast<FieldFunction<ValType>&>(*this)();
    }

    
    template <typename ValType>
    void
    MAST::FieldFunctionBase::partial (const MAST::FieldFunctionBase& f,
                                      const libMesh::Point& p,
                                      const libMesh::Real t,
                                      ValType& v) const {
        dynamic_cast<MAST::FieldFunction<ValType>&>(*this).partial(f, p, t, v);
    }

    
    template <typename ValType>
    void
    MAST::FieldFunctionBase::total (const MAST::FieldFunctionBase& f,
                                               const libMesh::Point& p,
                                               const libMesh::Real t,
                                               ValType& v) const {
        dynamic_cast<MAST::FieldFunction<ValType>&>(*this).total(f, p, t, v);
    }

    
}

#endif // __MAST_function_base_h__
