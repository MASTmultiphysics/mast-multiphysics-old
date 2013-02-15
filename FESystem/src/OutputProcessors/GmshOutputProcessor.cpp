//
//  GmshOutputProcessor.cpp
//  FESystem
//
//  Created by Manav Bhatia on 4/18/12.
//  Copyright (c) 2012. All rights reserved.
//

// C++ includes
#include <iomanip>

// FESystem includes
#include "OutputProcessors/GmshOutputProcessor.h"
#include "Numerics/VectorBase.h"
#include "DegreesOfFreedom/DegreeOfFreedomMap.h"
#include "DegreesOfFreedom/DegreeOfFreedomUnit.h"
#include "Mesh/MeshBase.h"
#include "Mesh/ElemBase.h"
#include "Mesh/Node.h"
#include "Base/FESystemExceptions.h"




FESystem::OutputProcessor::GmshOutputProcessor::GmshOutputProcessor()
{
    
}

FESystem::OutputProcessor::GmshOutputProcessor::~GmshOutputProcessor()
{
    
}
            

void
FESystem::OutputProcessor::GmshOutputProcessor::writeMesh(std::ostream& output, const FESystem::Mesh::MeshBase& mesh, const FESystem::DegreeOfFreedom::DegreeOfFreedomMap& dof_map)
{
    std::vector<FESystem::Mesh::ElementType> elem_type_vec;
    
    // write the file header
    output << "$MeshFormat" << std::endl
    << 2.2 << "  " << 0 << "  " << sizeof(double) << std::endl
    << "$EndMeshFormat" << std::endl;

    // header
    output << "$Nodes " << std::endl << mesh.getNNodes() << std::endl;

    const std::vector<FESystem::Mesh::Node*> nodes = mesh.getNodes();    
    for (FESystemUInt i=0; i<nodes.size(); i++)
    {
        output << std::setw(10) << i+1; 
        for (FESystemUInt j=0; j<3; j++)
            output << std::setw(15) << nodes[i]->getVal(j);
        output << std::endl;
    }
    output << "$EndNodes" << std::endl;

    output << "$Elements" << std::endl << mesh.getNElements() << std::endl;
    const std::vector<FESystem::Mesh::ElemBase*> elems = mesh.getElements();    
    FESystemUInt elem_type_num, n_nodes_to_write;
    for (FESystemUInt i=0; i<elems.size(); i++)
    {
        this->getGmshElemTypeNum(elems[i]->getElementType(), elem_type_num, n_nodes_to_write);
        output << std::setw(10) << i+1 << std::setw(10) << elem_type_num << std::setw(5) << 0; 
        for (FESystemUInt j=0; j<n_nodes_to_write; j++)
            output << std::setw(10) << elems[i]->getNode(j).getInternalID()+1;
        output << std::endl;
    }
    output << "$EndElements" << std::endl;
    
}
            


void 
FESystem::OutputProcessor::GmshOutputProcessor::getGmshElemTypeNum(FESystem::Mesh::ElementType type, FESystemUInt& elem_type_num, FESystemUInt& n_nodes_to_write)
{
    switch (type)
    {
        case FESystem::Mesh::EDGE2:
            elem_type_num = 1;
            n_nodes_to_write = 2;
            break;

        case FESystem::Mesh::EDGE3:
            elem_type_num = 8;
            n_nodes_to_write = 3;
            break;

        case FESystem::Mesh::TRI3:
            elem_type_num = 2;
            n_nodes_to_write = 3;
            break;
            
        case FESystem::Mesh::TRI6:
        case FESystem::Mesh::TRI7:
            elem_type_num = 9;
            n_nodes_to_write = 6;
            break;

        case FESystem::Mesh::QUAD4:
            elem_type_num = 3;
            n_nodes_to_write = 4;
            break;
            
        case FESystem::Mesh::QUAD9:
            elem_type_num = 10;
            n_nodes_to_write = 9;
            break;
            
        case FESystem::Mesh::TET4:
            elem_type_num = 4;
            n_nodes_to_write = 4;
            break;
            
        case FESystem::Mesh::HEX8:
            elem_type_num = 5;
            n_nodes_to_write = 8;
            break;
            
        case FESystem::Mesh::PRISM6:
            elem_type_num = 6;
            n_nodes_to_write = 6;
            break;
            
        case FESystem::Mesh::PYRAMID5:
            elem_type_num = 7;
            n_nodes_to_write = 5;
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
            
    }

}


void
FESystem::OutputProcessor::GmshOutputProcessor::writeSolution(std::ostream& output, 
                                                              const std::string& data_name,
                                                              const FESystem::Mesh::MeshBase& mesh, 
                                                              const FESystem::DegreeOfFreedom::DegreeOfFreedomMap& dof_map,
                                                              const std::vector<FESystemUInt>& variables_to_write,
                                                              const FESystem::Numerics::VectorBase<FESystemDouble>& vec)
{
}




