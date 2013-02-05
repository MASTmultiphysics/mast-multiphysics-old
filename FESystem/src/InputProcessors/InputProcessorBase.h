//
//  InputProcessorBase.h
//  FESystem
//
//  Created by Manav Bhatia on 2/4/13.
//
//

#ifndef __FESystem__InputProcessorBase__
#define __FESystem__InputProcessorBase__

// C++ includes
#include <istream>

// FESystem includes
#include "Base/FESystemTypes.h"

namespace FESystem
{
    // Forward declerations
    namespace Mesh {class MeshBase;}
    namespace Geometry {class Point;}
    
    namespace InputProcessor
    {
        class InputProcessorBase
        {
        public:
            
            virtual void readMeshFromInput(std::istream& input, const FESystemUInt mesh_dim, const FESystemBoolean if_local_physical_cs_same_as_global, FESystem::Geometry::Point& origin, FESystem::Mesh::MeshBase& mesh) = 0;
            
        protected:
            
        };
    }
}



#endif /* defined(__FESystem__InputProcessorBase__) */
