//
//  Structural1DElementBase.h
//  FESystem
//
//  Created by Manav Bhatia on 8/13/12.
//
//

#ifndef FESystem_Structural1DElementBase_h
#define FESystem_Structural1DElementBase_h


// FESystem includes
#include "Disciplines/Structure/StructuralElementBase.h"


namespace FESystem
{
    namespace Structures
    {
        
        class Structural1DElementBase: public FESystem::Structures::StructuralElementBase
        {
        public:
            Structural1DElementBase();
            
            virtual ~Structural1DElementBase();
            
        protected:

        };
    }
}

#endif
