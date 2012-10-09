
/*
 *  Node.cpp
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/27/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

// FESystem includes
#include "Mesh/Node.h"
#include "Mesh/ElemBase.h"


FESystem::Mesh::Node::Node(const FESystem::Geometry::CoordinateSystemBase& cs):
FESystem::Geometry::Point(cs),
FESystem::Base::DegreeOfFreedomObject()
{
    
}

FESystem::Mesh::Node::~Node()
{
    
}
    


void
FESystem::Mesh::Node::addElementConnection(FESystem::Mesh::ElemBase& el)
{
    // if the element is not already in the set, add it
    if (this->element_connections.count(&el) == 0)
        this->element_connections.insert(&el);
}


FESystemBoolean
FESystem::Mesh::Node::ifConnectedToElement(FESystem::Mesh::ElemBase& el)
{
    // if the element is not already in the set, add it
    if (this->element_connections.count(&el) == 0)
        return false;
    else
        return true;
}


const std::set<FESystem::Mesh::ElemBase*>&
FESystem::Mesh::Node::getElementConnectivitySet() const
{
    return this->element_connections;
}


FESystemUInt
FESystem::Mesh::Node::getNConnectedElements() const
{
    return this->element_connections.size();
}


FESystemUInt
FESystem::Mesh::Node::getNConnectedElementsOfDim(FESystemUInt dim) const
{
    std::set<FESystem::Mesh::ElemBase*>::const_iterator el_it = this->element_connections.begin(), el_end = this->element_connections.end();
    
    FESystemUInt n = 0;
    for ( ; el_it != el_end; el_it++)
        if ((*el_it)->getDimension() == dim)
            n++;

    return n;
}

