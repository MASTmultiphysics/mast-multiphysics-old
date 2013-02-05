//
//  GmshInputProcessor.h
//  FESystem
//
//  Created by Manav Bhatia on 2/4/13.
//
//

#ifndef __FESystem__GmshInputProcessor__
#define __FESystem__GmshInputProcessor__

// FESystem includes
#include "InputProcessors/InputProcessorBase.h"
#include "Mesh/ElementType.h"

namespace FESystem
{
    namespace InputProcessor
    {
        class GmshInputProcessor: public FESystem::InputProcessor::InputProcessorBase
        {
        public:
            GmshInputProcessor();
            
            virtual ~GmshInputProcessor();
            
            virtual void readMeshFromInput(std::istream& input, const FESystemUInt mesh_dim, const FESystemBoolean if_local_physical_cs_same_as_global, FESystem::Geometry::Point& origin, FESystem::Mesh::MeshBase& mesh);

        protected:

            void getGmshElemTypeNum(const FESystemUInt elem_type_num, FESystem::Mesh::ElementType& type);
            
        };
    }
}



#endif /* defined(__FESystem__GmshInputProcessor__) */
