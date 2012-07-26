//
//  IdObject.h
//  FESystem
//
//  Created by Manav Bhatia on 4/19/12.
//  Copyright (c) 2012 . All rights reserved.
//

#ifndef __fesystem_id_object_h__
#define __fesystem_id_object_h__

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

        protected:
            
            /*!
             *   internal ID
             */
            FESystemUInt internal_id;
            
            /*!
             *   external ID
             */
            FESystemUInt external_id;            
        };
    }
}


#endif // __fesystem_id_object_h__
