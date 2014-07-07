/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2014  Manav Bhatia
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

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
        _val(new ValType)
        { *_val = val;}
        
        ConstantFunction(const MAST::ConstantFunction<ValType>& f):
        FieldFunction<ValType>(f),
        _val(f._val)
        { }

        ~ConstantFunction() {
            // since the master function owns the pointer, it should not be deleted
            // if the function is not a master function.
            if (this == this->master())
                delete _val;
        }
        
        
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
        virtual void operator() (const libMesh::Point& p, const Real t, ValType& v) const
        {   v = *_val; }
        
        
        /*!
         *   returns the sensitivity of this function
         */
        virtual void partial (const MAST::FieldFunctionBase& f,
                              const libMesh::Point& p, const Real t,
                              ValType& v) const {
            this->total(f, p, t, v);
        }
        
        
        /*!
         *   returns the sensitivity of this function. This is the same as partial
         *   sensitivity for a constant function.
         */
        virtual void total (const MAST::FieldFunctionBase& f,
                            const libMesh::Point& p, const Real t,
                            ValType& v) const {
            // non-scalar values have a zero sensitivity
            (*this)(p, t, v);
            v.zero();
        }

        
        /*!
         *    Returns the pointer to value of this function.
         */
        virtual ValType* ptr()
        {   return _val; }
        
        
        
        /*!
         *  returns false since a constant function does not depend on any
         *  function.
         */
        virtual bool depends_on(const MAST::FieldFunctionBase& f) const {
            if (f.master() == this->master())
                return true;
            else
                return false;
        }
        
        
        
        /*!
         *  sets the value of this function
         */
        void operator =(const ValType& val)
        {   *_val = val; }
        
        
    protected:
        
        ValType* _val;
    };
    
    

    template < >
    inline void MAST::ConstantFunction<Real>::total (const MAST::FieldFunctionBase& f,
                                                     const libMesh::Point& p, const Real t,
                                                     Real& v) const {
        if (this->master() == f.master())
            v = 1.;
        else
            v = 0.;
    }
    
}

#endif // __MAST_constant_function_h__
