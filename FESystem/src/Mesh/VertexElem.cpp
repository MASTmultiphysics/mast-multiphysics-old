//
//  Vertex.cpp
//  FESystem
//
//  Created by Manav Bhatia on 2/14/13.
//
//

// FEsystem includes
#include "Mesh/VertexElem.h"


FESystem::Mesh::VertexElem::VertexElem(FESystemBoolean local_cs_same_as_global):
FESystem::Mesh::ElemBase(1, VERTEX1, local_cs_same_as_global)
{
    
}



FESystem::Mesh::VertexElem::~VertexElem()
{
    
}
