//
//  temperature_function.h
//  MAST
//
//  Created by Manav Bhatia on 10/23/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_temperature_function_h__
#define __MAST_temperature_function_h__

// MAST includes
#include ""


namespace MAST
{
    class Temperature {
    public:

        /*!
         *    virtual destructor
         */
        virtual ~Temperature();
        
        /*!
         *    virtual destructor
         */
        virtual void initialize() = 0;
        
        /*!
         *    virtual destructor
         */
        virtual Real value() = 0;
        
        /*!
         *    virtual destructor
         */
        virtual Real reference() = 0;
    };
}



#endif  // __MAST_temperature_function_h__
