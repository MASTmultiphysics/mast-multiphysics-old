//
//  Vertex.h
//  FESystem
//
//  Created by Manav Bhatia on 2/14/13.
//
//

#ifndef __FESystem__Vertex__
#define __FESystem__Vertex__

// FESystem includes
#include "Mesh/ElemBase.h"


namespace FESystem
{
    namespace Mesh
    {
        class VertexElem: public FESystem::Mesh::ElemBase
        {
        public:
            VertexElem(FESystemBoolean local_cs_same_as_global);
            
            virtual ~VertexElem();
            
        protected:
            
            
        };
    }
}


#endif /* defined(__FESystem__Vertex__) */
