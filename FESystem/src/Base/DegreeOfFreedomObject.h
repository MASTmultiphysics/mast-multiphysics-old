//
//  DegreeOfFreedomObject.h
//  FESystem
//
//  Created by Manav Bhatia on 4/6/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem_degree_of_freedom_object_h__
#define __fesystem_degree_of_freedom_object_h__

// C++ includes
#include <memory>

// FESystem includes
#include "Base/FESystemTypes.h"


namespace FESystem
{
    // Forward decleration
    namespace Utility {template <typename ValType> class AutoPtrTable;}
    
    namespace Base
    {
        // Forward declerations
        class DegreeOfFreedomUnit;
        
        /*!
         *   This class defines contains information about the degree of freedoms. The mesh entities like nodes and elements all 
         *   derive from this class. 
         */
        class DegreeOfFreedomObject
        {
        public:
            /*!
             *   Constructor
             */
            DegreeOfFreedomObject();
            
            ~DegreeOfFreedomObject();

            /*!
             *    Initializes the table to the required size of number of variables \p n_vars
             */
            void init(FESystemUInt n_vars);
            
            /*!
             *   returns the number of degree of freedom units in this degree of freedom object
             */
            FESystemUInt getNDegreeOfFreedomUnits() const;
            
            /*!
             *   Adds a new DegreeOfFreedomUnit for variable \p var_id to this DegreeOfFreedomObject. 
             *   It is an error to add new variables to the same system and variable ids.
             */
            FESystem::Base::DegreeOfFreedomUnit& addDegreeOfFreedomUnit(FESystemUInt var_id);
            
            /*!
             *   Returns the DegreeOfFreedomUnit for variable \p var_id in this DegreeOfFreedomObject. 
             *   An exception is thrown if the variable does not exist
             */
            FESystem::Base::DegreeOfFreedomUnit& getDegreeOfFreedomUnit(FESystemUInt var_id);

            /*!
             *   Returns the DegreeOfFreedomUnit for variable \p var_id to this DegreeOfFreedomObject. 
             *   An exception is thrown if the variable does not exist
             */
            const FESystem::Base::DegreeOfFreedomUnit& getDegreeOfFreedomUnit(FESystemUInt var_id) const;
            
        protected:
            
            /*!
             *   Table of variables for this degree of freedom object. The table is two dimensional, where the 
             *   first dimension is the system id, and the second dimension is the variable number for the system in 
             *   question
             */
            std::auto_ptr< FESystem::Utility::AutoPtrTable<FESystem::Base::DegreeOfFreedomUnit> > variable_table;

        };
    }
}



#endif  // __fesystem_degree_of_freedom_object_h__
