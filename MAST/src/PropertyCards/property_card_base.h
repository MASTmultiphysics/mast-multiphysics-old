//
//  property_card_base.h
//  MAST
//
//  Created by Manav Bhatia on 10/15/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_property_card_base_h__
#define __MAST_property_card_base_h__


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
        virtual ~PropertyCardBase() {
            // delete the functions in this card
            std::map<std::string, MAST::FunctionBase*>::iterator
            it = _properties.begin(), end = _properties.end();
            
            for ( ; it != end; it++)
                delete it->second;
        }
        
        
        /*!
         *    adds the function to this card and returns a reference to it.
         */
        template <typename ValType>
        MAST::FunctionBase& add(const std::string& nm, MAST::FunctionType type) {
            // make sure that this funciton does not already exist
            libmesh_assert(!_properties.count(nm));
            
            FunctionBase* ptr = MAST::build_function<ValType>(nm, type).release();
            
            bool success = _properties.insert(std::pair<std::string, MAST::FunctionBase*>
                                              (nm, ptr)).second;
            libmesh_assert(success);
            
            return *ptr;
        }
        
        
        /*!
         *   returns a constant reference to the specified function
         */
        template <typename ValType>
        const MAST::FunctionValue<ValType>& get(const std::string& nm) const {
            std::map<std::string, MAST::FunctionBase*>::const_iterator it =
            _properties.find(nm);
            
            // make sure that this funciton exists
            libmesh_assert(it != _properties.end());
            
            return dynamic_cast<MAST::FunctionValue<ValType>&>(*(it->second));
        }
        
        
        /*!
         *   returns a writable reference to the specified function
         */
        MAST::FunctionBase& get(const std::string& nm) {
            std::map<std::string, MAST::FunctionBase*>::iterator it =
            _properties.find(nm);
            
            // make sure that this funciton exists
            libmesh_assert(it != _properties.end());
            
            return *(it->second);
        }
        
        
        
    protected:
        
        /*!
         *    map of the functions in this card
         */
        std::map<std::string, MAST::FunctionBase*> _properties;
    };
    
}


#endif // __MAST_property_card_base_h__
