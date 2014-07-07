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

#ifndef __MAST_property_card_base_h__
#define __MAST_property_card_base_h__


// C++ includes
#include <vector>

// MAST includes
#include "Numerics/function_base.h"

namespace MAST {
    
    /*!
     *   provides a methods to store property values
     */
    class PropertyCardBase {
        
    public:
        PropertyCardBase()
        { }
        
        /*!
         *   destructor deletes the function pointers
         */
        virtual ~PropertyCardBase()
        { }
        

        /*!
         *   checks if the card contains the specified property value
         */
        bool contains(const std::string& nm) const {
            std::map<std::string, MAST::FieldFunctionBase*>::const_iterator it =
            _properties.find(nm);
            
            // make sure that this funciton exists
            return (it != _properties.end());
        }

        
        
        /*!
         *    adds the function to this card and returns a reference to it.
         */
        void add(MAST::FieldFunctionBase& f) {
            // make sure that this funciton does not already exist
            libmesh_assert(!_properties.count(f.name()));
            
            bool success = _properties.insert(std::pair<std::string, MAST::FieldFunctionBase*>
                                              (f.name(), &f)).second;
            libmesh_assert(success);
        }
        
        
        /*!
         *   returns a constant reference to the specified function
         */
        template <typename ValType>
        const ValType& get(const std::string& nm) const {
            std::map<std::string, MAST::FieldFunctionBase*>::const_iterator it =
            _properties.find(nm);
            
            // make sure that this funciton exists
            libmesh_assert(it != _properties.end());
            
            return dynamic_cast<ValType&>(*(it->second));
        }
        
        
        /*!
         *   returns a writable reference to the specified function
         */
        template <typename ValType>
        ValType& get(const std::string& nm) {
            std::map<std::string, MAST::FieldFunctionBase*>::iterator it =
            _properties.find(nm);
            
            // make sure that this funciton exists
            libmesh_assert(it != _properties.end());
            
            return dynamic_cast<ValType&>(*(it->second));
        }
        
        
        /*!
         *  returns true if the property card depends on the function \p f
         */
        virtual bool depends_on(const MAST::FieldFunctionBase& f) const {
            // check with all the properties to see if any one of them is
            // dependent on the provided parameter, or is the parameter itself
            std::map<std::string, MAST::FieldFunctionBase*>::const_iterator
            it = _properties.begin(), end = _properties.end();
            for ( ; it!=end; it++) {
                if (it->second->depends_on(f))
                    return true;
            }
            
            // if it gets here, then there is no dependency
            return false;
        }

        
    protected:
        
        /*!
         *    map of the functions in this card
         */
        std::map<std::string, MAST::FieldFunctionBase*> _properties;
    };
    
}


#endif // __MAST_property_card_base_h__
