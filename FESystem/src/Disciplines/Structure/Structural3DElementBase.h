//
//  Structural3DElementBase.h
//  FESystem
//
//  Created by Manav Bhatia on 8/13/12.
//
//

#ifndef FESystem_Structural3DElementBase_h
#define FESystem_Structural3DElementBase_h

// FESystem includes
#include "Disciplines/Structure/StructuralElementBase.h"


namespace FESystem
{
    namespace Structures
    {
        
        class Structural3DElementBase: public FESystem::Structures::StructuralElementBase
        {
        public:
            Structural3DElementBase();
            
            virtual ~Structural3DElementBase();
            
        protected:
            
        };
    }
}



#endif
