//
//  VtkOutputProcessor.cpp
//  FESystem
//
//  Created by Manav Bhatia on 4/18/12.
//  Copyright (c) 2012. All rights reserved.
//

// C++ includes
#include <iomanip>

// FESystem includes
#include "OutputProcessors/VtkOutputProcessor.h"
#include "Numerics/VectorBase.h"
#include "Base/DegreeOfFreedomMap.h"
#include "Base/DegreeOfFreedomUnit.h"
#include "Mesh/MeshBase.h"
#include "Mesh/ElemBase.h"
#include "Mesh/Node.h"
#include "Base/FESystemExceptions.h"




FESystem::OutputProcessor::VtkOutputProcessor::VtkOutputProcessor():
n_previous_dofs(0)
{
    
}

FESystem::OutputProcessor::VtkOutputProcessor::~VtkOutputProcessor()
{
    
}


void
FESystem::OutputProcessor::VtkOutputProcessor::writeMesh(std::ostream& output, const FESystem::Mesh::MeshBase& mesh, const FESystem::Base::DegreeOfFreedomMap& dof_map)
{
    FESystemUInt n_nodes_to_write=0, elem_type_num=0;
    
    // write the file header
    output << "# vtk DataFile Version 2.0" << std::endl
    << "Analysis Data" << std::endl
    << "ASCII" << std::endl
    << "DATASET UNSTRUCTURED_GRID" << std::endl;
    
    const std::vector<FESystem::Mesh::Node*> nodes = mesh.getNodes();
    const std::vector<FESystem::Mesh::ElemBase*> elems = mesh.getElements();

    // write the point data
    output << "POINTS " << std::setw(20) << nodes.size() << " double" << std::endl;
    for (FESystemUInt i_node=0; i_node<nodes.size(); i_node++)
    {
        for (FESystemUInt i_dim=0; i_dim<3; i_dim++)
            output << std::setw(20) << nodes[i_node]->getVal(i_dim);
        output << std::endl;
    }
    
    // write the cell data
    // calculate the total number of nodal data for VTK input file
    FESystemUInt n_data = elems.size();
    for (FESystemUInt i_elem=0; i_elem<elems.size(); i_elem++)
    {
        this->getVtkElemTypeNum(elems[i_elem]->getElementType(), elem_type_num, n_nodes_to_write);
        n_data += n_nodes_to_write;
    }
    
    output << "CELLS " << std::setw(20) << elems.size() << " " << std::setw(20) <<  n_data << std::endl;
    for (FESystemUInt i_elem=0; i_elem<elems.size(); i_elem++)
    {
        this->getVtkElemTypeNum(elems[i_elem]->getElementType(), elem_type_num, n_nodes_to_write);
        output << std::setw(10) << n_nodes_to_write;
        for (FESystemUInt i_node=0; i_node<n_nodes_to_write; i_node++)
            output << std::setw(10) << elems[i_elem]->getNode(i_node).getInternalID();
        output << std::endl;
    }
    
    
    // write the cell types
    output << "CELL_TYPES " << std::setw(20) << elems.size() << std::endl;
    for (FESystemUInt i_elem=0; i_elem<elems.size(); i_elem++)
    {
        this->getVtkElemTypeNum(elems[i_elem]->getElementType(), elem_type_num, n_nodes_to_write);
        output << std::setw(20) << elem_type_num << std::endl;    
    }
    
    this->n_previous_dofs = 0;
}



void
FESystem::OutputProcessor::VtkOutputProcessor::writeSolution(std::ostream& output, 
                                                              const std::string& data_name,
                                                              const FESystem::Mesh::MeshBase& mesh, 
                                                              const FESystem::Base::DegreeOfFreedomMap& dof_map,
                                                              const std::vector<FESystemUInt>& variables_to_write,
                                                              const FESystem::Numerics::VectorBase<FESystemDouble>& vec)
{    
    const std::vector<FESystem::Mesh::Node*> nodes = mesh.getNodes();
    
    if (this->n_previous_dofs == 0)
        output << "POINT_DATA " << nodes.size() << std::endl;

    if (variables_to_write.size() == 3)
    {    
        output<< "VECTORS " << data_name << " double" << std::endl;
        
        for (FESystemUInt i_node=0; i_node<nodes.size(); i_node++)
        {
            for (FESystemUInt i_var=0; i_var<variables_to_write.size(); i_var++)
                output << std::setw(20) << vec.getVal(nodes[i_node]->getDegreeOfFreedomUnit(variables_to_write[i_var]).global_dof_id[0]) << " ";
            output << std::endl;
        }
        
        this->n_previous_dofs = 3;    
    }
    else 
    {    
        for (FESystemUInt i=0; i<variables_to_write.size(); i++)
        {
            output<< "SCALARS " << data_name << i << " double" << std::endl;
            output<< "LOOKUP_TABLE default" << std::endl;
            
            for (FESystemUInt i_node=0; i_node<nodes.size(); i_node++)
                output << std::setw(20) << vec.getVal(nodes[i_node]->getDegreeOfFreedomUnit(variables_to_write[i]).global_dof_id[0]) << std::endl;
        }
        
        this->n_previous_dofs = variables_to_write.size(); 
    }

}



void
FESystem::OutputProcessor::VtkOutputProcessor::getVtkElemTypeNum(FESystem::Mesh::ElementType type, FESystemUInt& elem_type_num, FESystemUInt& n_nodes_to_write)
{
    switch (type)
    {
        case FESystem::Mesh::EDGE2:
            elem_type_num = 3;
            n_nodes_to_write = 2;
            break;

        case FESystem::Mesh::EDGE3:
            elem_type_num = 21;
            n_nodes_to_write = 3;
            break;

        case FESystem::Mesh::TRI3:
            elem_type_num = 5;
            n_nodes_to_write = 3;
            break;

        case FESystem::Mesh::TRI6:
        case FESystem::Mesh::TRI7:
            elem_type_num = 22;
            n_nodes_to_write = 6;
            break;

        case FESystem::Mesh::QUAD4:
            elem_type_num = 9;
            n_nodes_to_write = 4;
            break;

        case FESystem::Mesh::QUAD9:
            elem_type_num = 23;
            n_nodes_to_write = 9;
            break;

        case FESystem::Mesh::TET4:
            elem_type_num = 10;
            n_nodes_to_write = 4;
            break;
            
        case FESystem::Mesh::HEX8:
            elem_type_num = 12;
            n_nodes_to_write = 8;
            break;
            
        case FESystem::Mesh::PRISM6:
            elem_type_num = 13;
            n_nodes_to_write = 6;
            break;
            
        case FESystem::Mesh::PYRAMID5:
            elem_type_num = 14;
            n_nodes_to_write = 5;
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
            
    }
}



