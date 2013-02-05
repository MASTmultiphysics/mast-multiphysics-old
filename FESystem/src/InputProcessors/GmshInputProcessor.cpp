//
//  GmshInputProcessor.cpp
//  FESystem
//
//  Created by Manav Bhatia on 2/4/13.
//
//


// FESystem includes
#include "InputProcessors/GmshInputProcessor.h"
#include "Mesh/MeshBase.h"
#include "Mesh/ElemBase.h"
#include "Mesh/Node.h"
#include "Base/FESystemExceptions.h"



FESystem::InputProcessor::GmshInputProcessor::GmshInputProcessor()
{
    
}



FESystem::InputProcessor::GmshInputProcessor::~GmshInputProcessor()
{
    
}



void
FESystem::InputProcessor::GmshInputProcessor::readMeshFromInput(std::istream& input, const FESystemUInt mesh_dim, const FESystemBoolean if_local_physical_cs_same_as_global, FESystem::Geometry::Point& origin, FESystem::Mesh::MeshBase& mesh)
{
    const FESystem::Geometry::CoordinateSystemBase& global_cs = origin.getCoordinateSystem();

    // read the header
    std::string tmp, str_needed;
    FESystemUInt n0, n1, n2;
    FESystemDouble d0, d1, d2;
    
    str_needed = "$MeshFormat"; input >> tmp; FESystemAssert2(tmp == str_needed, FESystem::Exception::InvalidTag, tmp, str_needed);
    d1 = 2.2; input >> d0; FESystemAssert0(d1 == d0, FESystem::Exception::InvalidValue);
    input >> n0; input >> n0; // read and discard these values
    str_needed = "$EndMeshFormat"; input >> tmp; FESystemAssert2(tmp == str_needed, FESystem::Exception::InvalidTag, tmp, str_needed);

    // read nodes
    str_needed = "$Nodes"; input >> tmp; FESystemAssert2(tmp == str_needed, FESystem::Exception::InvalidTag, tmp, str_needed);
    input >> n0;
    std::auto_ptr<std::vector<FESystem::Mesh::Node*> > nodes;
    nodes.reset(mesh.createNodes(n0, global_cs).release()); // create the nodes
    FESystem::Mesh::Node* node_p;
    
    for (FESystemUInt i=0; i<n0; i++)
    {
        input >> n1; // node id
        input >> d0; input >> d1; input >> d2; // coordinates
        node_p = (*nodes)[i];
        node_p->setVal(0, d0);
        node_p->setVal(1, d1);
        node_p->setVal(2, d2);
        node_p->setExternalID(n1);
    }
    str_needed = "$EndNodes"; input >> tmp; FESystemAssert2(tmp == str_needed, FESystem::Exception::InvalidTag, tmp, str_needed);

    // read elements
    str_needed = "$Elements"; input >> tmp; FESystemAssert2(tmp == str_needed, FESystem::Exception::InvalidTag, tmp, str_needed);
    input >> n0;
    // read the location of the stream so that this can be rewound to read the boundary elements
    std::streamoff pos = input.tellg();
    
    FESystem::Mesh::ElemBase* elem_p;
    FESystem::Mesh::ElementType elem_type;
    
    for (FESystemUInt i=0; i<n0; i++)
    {
        input >> n1; input >> n2; // id and type
        this->getGmshElemTypeNum(n2, elem_type);
        
        if (FESystem::Mesh::getElementDimension(elem_type) == mesh_dim) // this is an element and should be added to the mesh
        {
            elem_p = (*mesh.createElements(1, elem_type, if_local_physical_cs_same_as_global))[0];
            elem_p->setExternalID(n1); // set the external id
            
            // read the element tags
            input >> n1; // number of tags, but only the first tag is of interest
            for (FESystemUInt j=0; j<n1-1; j++) // the last tag is that of the elementary entity, which is not of interest
            {
                input >> n2;
                elem_p->setTag(n2);
            }
            input >> n2; // last tag is ignored
            
            // read the nodes that connect this element
            for (FESystemUInt j=0; j<FESystem::Mesh::getNNodesForElement(elem_type); j++)
            {
                input >> n1; // node id
                elem_p->setNode(j, *((*nodes)[n1-1]));
            }
        }
        else
        {
            // read and discard the remaining data
            input >> n1; // number of tags
            for (FESystemUInt j=0; j<n1; j++)
                input >> n2; // read and discard
            for (FESystemUInt j=0; j<FESystem::Mesh::getNNodesForElement(elem_type); j++)
                input >> n2; // read and discard
        }
    }
    
    // now reread the boundary elements and set the tags for each of the boundaries
    input.seekg(pos);
    std::vector<FESystemUInt> tags;
    std::set<FESystem::Mesh::Node*> boundary_nodes;
    
    for (FESystemUInt i=0; i<n0; i++)
    {
        input >> n1; input >> n1; // id and type; id of the boundary element is ignored
        this->getGmshElemTypeNum(n1, elem_type);
        
        if (FESystem::Mesh::getElementDimension(elem_type) == mesh_dim-1) // this is a boundary, and its tag should be set for the specific boundary of the element to which it is connected
        {
            // read the element tags
            input >> n1; // number of tags, but only the first tag is of interest
            tags.resize(n1);
            for (FESystemUInt j=0; j<n1-1; j++) // the last tag is that of the elementary entity, which is not of interest
            {
                input >> n2; tags[j] = n2;
            }
            input >> n2; // last tag is ignored
            
            // read the nodes that connect this element
            boundary_nodes.clear();
            for (FESystemUInt j=0; j<FESystem::Mesh::getNNodesForElement(elem_type); j++)
            {
                input >> n1; // node id
                boundary_nodes.insert((*nodes)[n1-1]);
            }
            // now tell the mesh to set tags for element with boundary condition
            mesh.setTagsForBoundaryWithNodes(boundary_nodes, tags);
        }
        else
            break;
    }
    //str_needed = "$EndElements"; input >> tmp; FESystemAssert2(tmp == str_needed, FESystem::Exception::InvalidTag, tmp, str_needed);

    mesh.reinit();
}



void
FESystem::InputProcessor::GmshInputProcessor::getGmshElemTypeNum(const FESystemUInt elem_type_num, FESystem::Mesh::ElementType& type)
{
    switch (elem_type_num)
    {
        case 1:
            type = FESystem::Mesh::EDGE2;
            break;
            
        case 8:
            type = FESystem::Mesh::EDGE3;
            break;
            
        case 2:
            type = FESystem::Mesh::TRI3;
            break;
            
        case 9:
            type = FESystem::Mesh::TRI6;
            break;
            
        case 3:
            type = FESystem::Mesh::QUAD4;
            break;
            
        case 10:
            type = FESystem::Mesh::QUAD9;
            break;
            
        case 4:
            type = FESystem::Mesh::TET4;
            break;
            
        case 5:
            type = FESystem::Mesh::HEX8;
            break;
            
        case 6:
            type = FESystem::Mesh::PRISM6;
            break;
            
        case 7:
            type = FESystem::Mesh::PYRAMID5;
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
            
    }
    
}
