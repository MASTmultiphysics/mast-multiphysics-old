
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
#include "Base/FESystemExceptions.h"


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


void
FESystem::Mesh::MeshBase::setTagsForBoundaryWithNodes(std::set<FESystem::Mesh::Node*>& n, const std::vector<FESystemUInt>& tags)
{
    // iterate over the nodes and the elements to which this boundary is connected to. Note that there can be atmost one element to which
    // this boundary is connected
    
    FESystemUInt n_elems_sharing_given_boundary = 0;
    FESystem::Mesh::ElemBase* elem_p = NULL;
    
    std::set<FESystem::Mesh::ElemBase*>::const_iterator it, end;
    FESystemUInt b_id = 0;
    FESystemBoolean if_elem_has_given_boundary = false;
    
    std::set<FESystem::Mesh::Node*>::iterator n_it = n.begin(), n_end = n.end();
    
    for ( ; n_it != n_end; n_it++)
    {
        // get element that the node is connected to
        const std::set<FESystem::Mesh::ElemBase*>& elems = (*n_it)->getElementConnectivitySet();
        it = elems.begin(); end = elems.end();
        for ( ; it != end; it++)
        {
            (*it)->getIDOfBoundaryWithNodes(n, if_elem_has_given_boundary, b_id);
            if (if_elem_has_given_boundary && (elem_p != *it))
            {
                FESystemAssert0(n_elems_sharing_given_boundary == 0, FESystem::Exception::InvalidValue); // only one element can have this boundary
                n_elems_sharing_given_boundary++;
                elem_p = *it; // the code should get here only once
                for (FESystemUInt j=0; j<tags.size(); j++)
                    (*it)->setTagForBoundary(b_id, tags[j]);
            }
        }
    }
}



FESystem::Mesh::Node&
FESystem::Mesh::MeshBase::createNode(const FESystem::Geometry::CoordinateSystemBase& cs)
{
    this->nodes.push_back(new FESystem::Mesh::Node(cs));
    return **(this->nodes.end());
}
    
    


FESystemUInt
FESystem::Mesh::getNNodesForElement(FESystem::Mesh::ElementType type)
{
    switch (type)
    {
        case FESystem::Mesh::EDGE2:
            return 2;
            break;

        case FESystem::Mesh::EDGE3:
            return 3;
            break;

        case FESystem::Mesh::EDGE5:
            return 5;
            break;

        case FESystem::Mesh::QUAD4:
            return 4;
            break;

        case FESystem::Mesh::QUAD8:
            return 8;
            break;

        case FESystem::Mesh::QUAD9:
            return 9;
            break;

        case FESystem::Mesh::TRI3:
            return 3;
            break;

        case FESystem::Mesh::TRI6:
            return 6;
            break;

        case FESystem::Mesh::TRI7:
            return 7;
            break;

        case FESystem::Mesh::HEX8:
            return 8;
            break;

        case FESystem::Mesh::HEX27:
            return 27;
            break;

        case FESystem::Mesh::TET4:
            return 4;
            break;

        case FESystem::Mesh::PRISM6:
            return 6;
            break;

        case FESystem::Mesh::PYRAMID5:
            return 5;
            break;

        default:
            FESystemAssert1(false, FESystem::Exception::EnumerationNotHandled, type);
            break;
    }
}



FESystemUInt
FESystem::Mesh::getElementDimension(FESystem::Mesh::ElementType type)
{
    switch (type)
    {
        case FESystem::Mesh::EDGE2:
        case FESystem::Mesh::EDGE3:
        case FESystem::Mesh::EDGE5:
            return 1;
            break;
            
        case FESystem::Mesh::QUAD4:
        case FESystem::Mesh::QUAD8:
        case FESystem::Mesh::QUAD9:
        case FESystem::Mesh::TRI3:
        case FESystem::Mesh::TRI6:
        case FESystem::Mesh::TRI7:
            return 2;
            break;
            
        case FESystem::Mesh::HEX8:
        case FESystem::Mesh::HEX27:
        case FESystem::Mesh::TET4:
        case FESystem::Mesh::PRISM6:
        case FESystem::Mesh::PYRAMID5:
            return 3;
            break;
            
        default:
            FESystemAssert1(false, FESystem::Exception::EnumerationNotHandled, type);
            break;
    }
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



