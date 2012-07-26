//
//  TecplotOutputProcessor.cpp
//  FESystem
//
//  Created by Manav Bhatia on 4/18/12.
//  Copyright (c) 2012. All rights reserved.
//

// C++ includes
#include <iomanip>

// FESystem includes
#include "OutputProcessors/TecplotOutputProcessor.h"
#include "Numerics/VectorBase.h"
#include "Base/DegreeOfFreedomMap.h"
#include "Base/DegreeOfFreedomUnit.h"
#include "Mesh/MeshBase.h"
#include "Mesh/ElemBase.h"
#include "Mesh/Node.h"
#include "Base/FESystemExceptions.h"



FESystem::OutputProcessor::TecplotOutputProcessor::TecplotOutputProcessor():
previous_var_zone_number(1)
{
    
}


FESystem::OutputProcessor::TecplotOutputProcessor::~TecplotOutputProcessor()
{
    
}


void
FESystem::OutputProcessor::TecplotOutputProcessor::writeMesh(std::ostream& output, const FESystem::Mesh::MeshBase& mesh, const FESystem::Base::DegreeOfFreedomMap& dof_map)
{
    std::vector<FESystem::Mesh::ElementType> elem_type_map;
    this->populateElemTypeVector(mesh, elem_type_map);
    
    FESystemUInt n_dofs = dof_map.getNVariables();
    output  << "TITLE = \" FESystem Output\" " << std::endl; 
    output << "VARIABLES = \"X\",  \"Y\",  \"Z\"";    // list all the variables: first three are x,y,z and then the variables of this analysis. 
    for (FESystemUInt i_var=0; i_var<n_dofs; i_var++)
        output << ",  \"V" << i_var << "\" ";
    
    std::string zone_name, elem_name;
    
    const std::vector<FESystem::Mesh::Node*>& nodes = mesh.getNodes();
    const std::vector<FESystem::Mesh::ElemBase*>& elems = mesh.getElements();
    
    // now create the zones based on the element types in the mesh
    for (FESystemUInt i_types=0; i_types<elem_type_map.size(); i_types++)
    {
        this->getFiniteElementZoneName(elem_type_map[i_types], zone_name, elem_name);
        
        output << "\n ZONE  T=\"Mesh\", NODES=" << mesh.getNNodes() << ",  ELEMENTS=" << mesh.getNElemsOfType(elem_type_map[i_types]) 
        << ", DATAPACKING=POINT,  ZONETYPE="<< elem_name;// << ", PASSIVEVARLIST=(["<< 4 <<"-"<< 3+n_dofs <<"])";
        
        if (i_types > 0)
            output << ", VARSHARELIST=[1-"<< 3+n_dofs <<"]=1" << std::endl; // later zones use the node locations from the first zone
        else
        {
            // write the entire nodal data for the first zone
            output << std::endl;
            for (FESystemUInt i_node=0; i_node<nodes.size(); i_node++)
            {
                for (FESystemUInt i_dim=0; i_dim<3; i_dim++)
                    output << std::setw(20) << nodes[i_node]->getVal(i_dim) << "  ";
                // write dummy values for all variables
                for (FESystemUInt i_var=0; i_var<n_dofs; i_var++)
                    output << "  0.0";
                output << std::endl;
            }
        }
        
        // now, write the connectivity for this element zone
        for (FESystemUInt i_elem=0; i_elem<elems.size(); i_elem++)
            if (elems[i_elem]->getElementType() == elem_type_map[i_types]) 
            {
                for (FESystemUInt i_node=0; i_node<elems[i_elem]->getNNodes(); i_node++)
                    output << std::setw(10) << 1+elems[i_elem]->getNode(i_node).getInternalID();
                output << std::endl;
            }
    }
    
}



void
FESystem::OutputProcessor::TecplotOutputProcessor::writeSolution(std::ostream& output, 
                                                                 const std::string& data_name,
                                                                 const FESystem::Mesh::MeshBase& mesh, 
                                                                 const FESystem::Base::DegreeOfFreedomMap& dof_map,
                                                                 const std::vector<FESystemUInt>& variables_to_write,
                                                                 const FESystem::Numerics::VectorBase<FESystemDouble>& vec)
{
    
    std::vector<FESystem::Mesh::ElementType> elem_type_map;
    this->populateElemTypeVector(mesh, elem_type_map);
    
    FESystemUInt n_dofs = dof_map.getNVariables();
    FESystemAssert0(variables_to_write.size() == n_dofs, FESystem::Exception::InvalidValue); // this output processor writes all the data

    std::string zone_name, elem_name;
    
    const std::vector<FESystem::Mesh::Node*>& nodes = mesh.getNodes();
    
    // now create the zones based on the element types in the mesh
    for (FESystemUInt i_types=0; i_types<elem_type_map.size(); i_types++)
    {
        this->getFiniteElementZoneName(elem_type_map[i_types], zone_name, elem_name);
        
        output << "\n ZONE  T=\""<<data_name<<"\", NODES=" << mesh.getNNodes() << ",  ELEMENTS=" << mesh.getNElemsOfType(elem_type_map[i_types]) 
        << ", DATAPACKING=POINT,  ZONETYPE="<< elem_name << ", VARSHARELIST=(["<< 1 <<"-"<< 3<<"]=1), CONNECTIVITYSHAREZONE="<<i_types+1;
        
        if (i_types > 0)
            output << ", VARSHARELIST=[4-"<<3+n_dofs<<"]="<<this->previous_var_zone_number<< std::endl; // later zones use the node locations from the first zone
        else
        {
            // write the entire nodal data for the first zone
            output << std::endl;
            for (FESystemUInt i_node=0; i_node<nodes.size(); i_node++)
            {
                for (FESystemUInt i_var=0; i_var<n_dofs; i_var++)
                    output << std::setw(20) << vec.getVal(nodes[i_node]->getDegreeOfFreedomUnit(i_var).global_dof_id[0]) << "  ";
                output << std::endl;
            }
            this->previous_var_zone_number += elem_type_map.size();
        }
        
    }
    
}



void 
FESystem::OutputProcessor::TecplotOutputProcessor::getFiniteElementZoneName(FESystem::Mesh::ElementType elem_type, 
                                                                            std::string& zone_name,
                                                                            std::string& elem_name)
{
    switch (elem_type)
    {
        case FESystem::Mesh::EDGE2:
        {
            zone_name = "EDGE ELEMS";
            elem_name = "FELINESEG";
        }
            break;
            
        case FESystem::Mesh::TRI3:
        {
            zone_name = "TRI ELEMS";
            elem_name = "FETRIANGLE";
        }
            break;
            
        case FESystem::Mesh::QUAD4:
        {
            zone_name = "QUAD ELEMS";
            elem_name = "FEQUADRILATERAL";
        }
            break;            
            
            //tecplot treats all linear 3-D elements as bricks. 
        case FESystem::Mesh::HEX8:
        case FESystem::Mesh::PRISM6:
        case FESystem::Mesh::PYRAMID5:
        {
            zone_name = "HEX8 ELEMS";
            elem_name = "FEBRICK";
        }
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
    }
}


void
FESystem::OutputProcessor::TecplotOutputProcessor::populateElemTypeVector(const FESystem::Mesh::MeshBase& mesh,
                                                                          std::vector<FESystem::Mesh::ElementType>& elem_type_vec)
{
    std::vector<FESystem::Mesh::ElementType> vec;
    
    vec.push_back(FESystem::Mesh::EDGE2);
    vec.push_back(FESystem::Mesh::TRI3);
    vec.push_back(FESystem::Mesh::QUAD4);
    vec.push_back(FESystem::Mesh::TET4);
    vec.push_back(FESystem::Mesh::HEX8);
    vec.push_back(FESystem::Mesh::PRISM6);
    vec.push_back(FESystem::Mesh::PYRAMID5);
    
    // no second order elements are added
    vec.push_back(FESystem::Mesh::EDGE3);
    vec.push_back(FESystem::Mesh::TRI7);
    vec.push_back(FESystem::Mesh::QUAD9);
    vec.push_back(FESystem::Mesh::TET15);
    vec.push_back(FESystem::Mesh::HEX26);
    vec.push_back(FESystem::Mesh::PRISM21);
    vec.push_back(FESystem::Mesh::PYRAMID19);

    FESystemUInt n_elem_of_type = 0;

    // iterate over the number of elements and if a certain type exists in the mesh, add it to the vector to be returned
    for (FESystemUInt i_type=0; i_type<vec.size(); i_type++)
    {
        if (mesh.getNElemsOfType(vec[i_type]))
            elem_type_vec.push_back(vec[i_type]);
    }
}


