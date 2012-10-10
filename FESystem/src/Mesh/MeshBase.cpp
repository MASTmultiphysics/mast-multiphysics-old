
/*
 *  MeshBase.cpp
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/15/09.
 *  Copyright 2009 . All rights reserved.
 *
 */


// FESystem includes
#include "Mesh/MeshBase.h"
#include "Mesh/Node.h"
#include "Mesh/ElemBase.h"


FESystem::Mesh::MeshBase::MeshBase():
if_initialized(false)
{
    
}
    


FESystem::Mesh::MeshBase::~MeshBase()
{
    // iterate over all the nodes  and delete them
    {
        std::vector<FESystem::Mesh::Node*>::iterator it, end;
        it = this->nodes.begin();
        end = this->nodes.end();
        
        for ( ; it != end; it++)
            if ( *it != NULL)
                delete *it;
    }
    
    // iterate over all the elements and delete them
    {
        std::vector<FESystem::Mesh::ElemBase*>::iterator it, end;
        it = this->elements.begin();
        end = this->elements.end();
        
        for ( ; it != end; it++)
            if ( *it != NULL)
                delete *it;
    }
    
}
 


void
FESystem::Mesh::MeshBase::reinit()
{
    this->if_initialized = false;

    // node maps
    this->internal_id_to_node_map.clear();
    this->external_id_to_node_map.clear();

    // elem maps
    this->internal_id_to_elem_map.clear();
    this->external_id_to_elem_map.clear();
    
    // go over all nodes and elements and give them internal IDs and process all the external IDs by 
    // populating the map. 
    FESystemBoolean insert_success;
    FESystemUInt id_num = 0;
    {
        FESystem::Mesh::Node* n_ptr = NULL;
        for (FESystemUInt i_node=0; i_node<this->nodes.size(); i_node++)
        {
            n_ptr = this->nodes[i_node];
            n_ptr->setInternalID(id_num);
            // add node internal ID
            insert_success = this->internal_id_to_node_map.insert(std::map<FESystemUInt, FESystem::Mesh::Node*>::value_type(id_num, n_ptr)).second;
            FESystemAssert1(insert_success, FESystem::Exception::InvalidID, id_num);
            
            // add node external ID
            insert_success = this->external_id_to_node_map.insert(std::map<FESystemUInt, FESystem::Mesh::Node*>::value_type(n_ptr->getExternalID(), n_ptr)).second;
            FESystemAssert1(insert_success, FESystem::Exception::InvalidID, n_ptr->getExternalID());
            id_num++;
        }
    }
    
    // now the elements
    id_num = 0;
    {
        FESystem::Mesh::ElemBase* e_ptr = NULL;
        std::pair<std::map<FESystemUInt, FESystem::Mesh::ElemBase*>::iterator, FESystemBoolean > elem_insert_return;
        for (FESystemUInt i_elem=0; i_elem<this->elements.size(); i_elem++)
        {
            e_ptr = this->elements[i_elem];
            e_ptr->setInternalID(id_num);
            // add node internal ID
            insert_success = this->internal_id_to_elem_map.insert(std::map<FESystemUInt, FESystem::Mesh::ElemBase*>::value_type(id_num, e_ptr)).second;
            FESystemAssert1(insert_success, FESystem::Exception::InvalidID, id_num);
            
            // add node external ID
            insert_success = this->external_id_to_elem_map.insert(std::map<FESystemUInt, FESystem::Mesh::ElemBase*>::value_type(e_ptr->getExternalID(), e_ptr)).second;
            FESystemAssert1(insert_success, FESystem::Exception::InvalidID, e_ptr->getExternalID());
            id_num++;
        }
    }

    this->if_initialized = true;
}




FESystem::Mesh::Node& 
FESystem::Mesh::MeshBase::getNodeFromInternalID(FESystemUInt i)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    std::map<FESystemUInt, FESystem::Mesh::Node*>::iterator
    it = this->internal_id_to_node_map.find(i), 
    end = this->internal_id_to_node_map.end();
    
    FESystemAssert1(it != end, FESystem::Exception::InvalidID, i);
    return *(it->second);
}


FESystem::Mesh::Node& 
FESystem::Mesh::MeshBase::getNodeFromExternalID(FESystemUInt i)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    std::map<FESystemUInt, FESystem::Mesh::Node*>::iterator
    it = this->external_id_to_node_map.find(i), 
    end = this->external_id_to_node_map.end();
    
    FESystemAssert1(it != end, FESystem::Exception::InvalidID, i);
    return *(it->second);
}


const FESystem::Mesh::Node& 
FESystem::Mesh::MeshBase::getNodeFromInternalID(FESystemUInt i) const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    std::map<FESystemUInt, FESystem::Mesh::Node*>::const_iterator
    it = this->internal_id_to_node_map.find(i), 
    end = this->internal_id_to_node_map.end();
    
    FESystemAssert1(it != end, FESystem::Exception::InvalidID, i);
    return *(it->second);
}


const FESystem::Mesh::Node& 
FESystem::Mesh::MeshBase::getNodeFromExternalID(FESystemUInt i) const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    std::map<FESystemUInt, FESystem::Mesh::Node*>::const_iterator
    it = this->external_id_to_node_map.find(i), 
    end = this->external_id_to_node_map.end();
    
    FESystemAssert1(it != end, FESystem::Exception::InvalidID, i);
    return *(it->second);
}


FESystem::Mesh::ElemBase& 
FESystem::Mesh::MeshBase::getElemFromInternalID(FESystemUInt i)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    std::map<FESystemUInt, FESystem::Mesh::ElemBase*>::iterator
    it = this->internal_id_to_elem_map.find(i), 
    end = this->internal_id_to_elem_map.end();
    
    FESystemAssert1(it != end, FESystem::Exception::InvalidID, i);
    return *(it->second);    
}


FESystem::Mesh::ElemBase& 
FESystem::Mesh::MeshBase::getElemFromExternalID(FESystemUInt i)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    std::map<FESystemUInt, FESystem::Mesh::ElemBase*>::iterator
    it = this->external_id_to_elem_map.find(i), 
    end = this->external_id_to_elem_map.end();
    
    FESystemAssert1(it != end, FESystem::Exception::InvalidID, i);
    return *(it->second);    
}


const FESystem::Mesh::ElemBase& 
FESystem::Mesh::MeshBase::getElemFromInternalID(FESystemUInt i) const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    std::map<FESystemUInt, FESystem::Mesh::ElemBase*>::const_iterator
    it = this->internal_id_to_elem_map.find(i), 
    end = this->internal_id_to_elem_map.end();
    
    FESystemAssert1(it != end, FESystem::Exception::InvalidID, i);
    return *(it->second);    
}


const FESystem::Mesh::ElemBase& 
FESystem::Mesh::MeshBase::getElemFromExternalID(FESystemUInt i) const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    std::map<FESystemUInt, FESystem::Mesh::ElemBase*>::const_iterator
    it = this->external_id_to_elem_map.find(i), 
    end = this->external_id_to_elem_map.end();
    
    FESystemAssert1(it != end, FESystem::Exception::InvalidID, i);
    return *(it->second);    
}



FESystemUInt 
FESystem::Mesh::MeshBase::getNNodes() const
{
    return this->nodes.size();
}
    
    



FESystemUInt 
FESystem::Mesh::MeshBase::getNElements() const
{
    return this->elements.size();
}
    


FESystemUInt 
FESystem::Mesh::MeshBase::getNElemsOfType(FESystem::Mesh::ElementType t) const
{
    FESystemUInt n=0;
    for (FESystemUInt i_elem=0; i_elem<this->elements.size(); i_elem++)
        if (this->elements[i_elem]->getElementType() == t)
            n++;
    return n;
}


std::vector<FESystem::Mesh::Node*>&
FESystem::Mesh::MeshBase::getNodes() 
{
    return this->nodes;
}


const std::vector<FESystem::Mesh::Node*>&
FESystem::Mesh::MeshBase::getNodes() const
{
    return this->nodes;
}



std::vector<FESystem::Mesh::ElemBase*>&
FESystem::Mesh::MeshBase::getElements() 
{
    return this->elements;
}



const std::vector<FESystem::Mesh::ElemBase*>&
FESystem::Mesh::MeshBase::getElements() const
{
    return this->elements;
}


FESystem::Mesh::Node&
FESystem::Mesh::MeshBase::createNode(const FESystem::Geometry::CoordinateSystemBase& cs)
{
    this->nodes.push_back(new FESystem::Mesh::Node(cs));
    return **(this->nodes.end());
}
    
    


std::auto_ptr<std::vector<FESystem::Mesh::Node*> >
FESystem::Mesh::MeshBase::createNodes(FESystemUInt n_nodes, const FESystem::Geometry::CoordinateSystemBase& cs)
{
    std::auto_ptr<std::vector<FESystem::Mesh::Node*> > rval(new std::vector<FESystem::Mesh::Node*>(n_nodes));
    for (FESystemUInt i=0; i<n_nodes; i++)
    {
        this->nodes.push_back(new FESystem::Mesh::Node(cs));
        (*rval)[i] = *(this->nodes.rbegin());
    }
    
    return rval;
}
    


std::auto_ptr<std::vector<FESystem::Mesh::ElemBase*> > 
FESystem::Mesh::MeshBase::createElements(FESystemUInt n_elems, FESystem::Mesh::ElementType elem_type, FESystemBoolean if_local_physical_cs_same_as_global)
{
    std::auto_ptr<std::vector<FESystem::Mesh::ElemBase*> > rval(new std::vector<FESystem::Mesh::ElemBase*>(n_elems));
    for (FESystemUInt i=0; i<n_elems; i++)
    {
        this->elements.push_back(FESystem::Mesh::ElementCreate(*this, elem_type, if_local_physical_cs_same_as_global).release());
        (*rval)[i] = *(this->elements.rbegin());
    }
    
    return rval;
}


