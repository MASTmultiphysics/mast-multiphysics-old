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
         *   checks if the card contains the specified property value
         */
        bool contains(const std::string& nm) const {
            std::map<std::string, MAST::FunctionBase*>::const_iterator it =
            _properties.find(nm);
            
            // make sure that this funciton exists
            return (it != _properties.end());
        }

        
        
        /*!
         *    adds the function to this card and returns a reference to it.
         */
        template <typename ValType>
        MAST::FunctionValue<ValType>& add(const std::string& nm,
                                          MAST::FunctionType type) {
            // make sure that this funciton does not already exist
            libmesh_assert(!_properties.count(nm));
            
            MAST::FunctionValue<ValType>* ptr =
            MAST::build_function<ValType>(nm, type).release();
            
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
        template <typename ValType>
        MAST::FunctionValue<ValType>& get(const std::string& nm) {
            std::map<std::string, MAST::FunctionBase*>::iterator it =
            _properties.find(nm);
            
            // make sure that this funciton exists
            libmesh_assert(it != _properties.end());
            
            return dynamic_cast<MAST::FunctionValue<ValType>&>(*(it->second));
        }
        
        
        /*!
         *  returns true if alteast one property in this property card 
         *  depends on this parameter
         */
        virtual bool depends_on(Real* val) const {
            // check with all the properties to see if any one of them is
            // dependent on the provided parameter
            bool rval = false;
            std::map<std::string, MAST::FunctionBase*>::const_iterator
            it = _properties.begin(), end = _properties.end();
            for ( ; it!=end; it++) {
                rval = it->second->depends_on(val);
                if (rval)
                    break;
            }
            
            return rval;
        }

        
        /*!
         *  returns true if the property card depends on the function \p f
         */
        virtual bool depends_on(const FunctionBase& f) const {
            // check with all the properties to see if any one of them is
            // dependent on the provided parameter, or is the parameter itself
            bool rval = false;
            std::map<std::string, MAST::FunctionBase*>::const_iterator
            it = _properties.begin(), end = _properties.end();
            for ( ; it!=end; it++) {
                rval = (it->second->depends_on(f) ||
                        it->second == &f);
                if (rval)
                    break;
            }
            
            return rval;
        }

        
        /*!
         *  returns true if the property card depends on the functions in \p p
         */
        virtual bool depends_on(const MAST::SensitivityParameters& p) const {
            
            const MAST::SensitivityParameters::ParameterMap& p_map = p.get_map();
            MAST::SensitivityParameters::ParameterMap::const_iterator it, end;
            it = p_map.begin(); end = p_map.end();
            
            std::vector<bool> rvals(p_map.size());
            std::fill(rvals.begin(), rvals.end(), false);
            unsigned int i=0;
            bool rval = true;
            
            for ( ; it != end; it++) {
                const MAST::FunctionBase& f = *(it->first);
                // check with all the properties to see if any one of them is
                // dependent on the provided parameter, or is the parameter itself
                std::map<std::string, MAST::FunctionBase*>::const_iterator
                it = _properties.begin(), end = _properties.end();
                for ( ; it!=end; it++) {
                    rvals[i] = (it->second->depends_on(f) ||
                            it->second == &f);
                    if (rvals[i])
                        break;
                }

                rval = rval && rvals[i];

                if (!rval)
                    return false;
                
                // increment the counter
                i++;
            }
            
            return rval;
        }

        
    protected:
        
        /*!
         *    map of the functions in this card
         */
        std::map<std::string, MAST::FunctionBase*> _properties;
    };
    
}


#endif // __MAST_property_card_base_h__
