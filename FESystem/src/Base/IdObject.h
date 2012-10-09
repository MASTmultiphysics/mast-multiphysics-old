//
//  IdObject.h
//  FESystem
//
//  Created by Manav Bhatia on 4/19/12.
//  Copyright (c) 2012 . All rights reserved.
//

#ifndef __fesystem_id_object_h__
#define __fesystem_id_object_h__

// C++ includes
#include <set>

// FESystem includes
#include "Base/FESystemTypes.h"


namespace FESystem
{
    namespace Base
    {
        /*!
         *   ID object that allows setting internal and external IDs
         */
        class IdObject
        {
        public:
            /*!
             *   Constructor
             */
            IdObject():
            internal_id(0),
            external_id(0)
            {
                
            }
            
            ~IdObject()
            {
                
            }
            
            /*!
             *   sets the internal ID
             */
            void setInternalID(FESystemUInt i)
            {
                this->internal_id = i;
            }

            /*!
             *   sets the external ID
             */
            void setExternalID(FESystemUInt i)
            {
                this->external_id = i;
            }

            /*!
             *   returns the internal ID
             */
            FESystemUInt getInternalID() const
            {
                return this->internal_id;
            }

            /*!
             *   returns the external ID
             */
            FESystemUInt getExternalID() const
            {
                return this->external_id;
            }
            
            
            /*!
             *    Set a user-defined tag
             */
            void setTag(FESystemInt n)
            {
                FESystemAssert0(!this->tags.count(n), FESystem::Exception::InvalidValue);
                this->tags.insert(n);
            }

            /*!
             *    removes a user-defined tag
             */
            void unsetTag(FESystemInt n)
            {
                FESystemAssert0(this->tags.count(n), FESystem::Exception::InvalidValue);
                this->tags.erase(n);
            }

            
            /*!
             *    checks if a user defined tags is set
             */
            FESystemBoolean checkForTag(FESystemInt n)
            {
                return this->tags.count(n);
            }

        protected:
            
            /*!
             *   internal ID
             */
            FESystemUInt internal_id;
            
            /*!
             *   external ID
             */
            FESystemUInt external_id;
            
            /*!
             *   Integer tags that the user can use to identify various properties
             */
            std::set<FESystemInt> tags;
        };
    }
}


#endif // __fesystem_id_object_h__
