//
//  DegreeOfFreedomUnit.h
//  FESystem
//
//  Created by Manav Bhatia on 4/6/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem_degree_of_freedom_unit_h__
#define __fesystem_degree_of_freedom_unit_h__

// C++ includes
#include <vector>

// FESystem includes
#include "Base/FESystemTypes.h"

namespace FESystem
{
    namespace Base 
    {
        /*!
         *   
         */
        class DegreeOfFreedomUnit
        {
        public:
            /*!
             *   system id of the degree of freedom unit  
             */
            FESystemUInt system_id;
            
            /*!
             *   internal variable id of the degree of freedom unit
             */
            FESystemUInt variable_id;
            
            /*!
             *    current p_level of the degree of freedom unit
             */
            FESystemUInt p_level;
            
            /*!
             *    if the p-level for this degree of freedom unit can be modified
             */
            FESystemBoolean if_p_adaptable;
            
            /*!
             *    vector of global ids for the individual p-levels. In case the 
             *    unit does not use p-level polynomials, then this is only a 1 dimensional 
             *    vector, which stores the global dof id of the variable
             */
            std::vector<FESystemUInt> global_dof_id;
        };
    }
}


#endif  // __fesystem_degree_of_freedom_unit_h__ 

