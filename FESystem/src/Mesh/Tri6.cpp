//
//  Tri6.cpp
//  FESystem
//
//  Created by Manav Bhatia on 5/21/12.
//  Copyright (c) 2012. All rights reserved.
//

// FESystem includes
#include "Mesh/Tri6.h"
#include "Mesh/Quad9.h"
#include "Mesh/Node.h"
#include "Geom/Point.h"
#include "Numerics/DenseMatrix.h"


// initialization of the pointer here
std::auto_ptr<FESystem::Numerics::MatrixBase<FESystemDouble> > FESystem::Mesh::Tri6::tri6_nondegenerate_to_degenerate_element_mapping(NULL);


FESystem::Mesh::Tri6::Tri6():
FESystem::Mesh::TriElemBase(6, FESystem::Mesh::TRI6),
center_node_for_parent_nondegenerate_elem(NULL)
{    
    // initialize the mapping matrix
    if (FESystem::Mesh::Tri6::tri6_nondegenerate_to_degenerate_element_mapping.get() == NULL)
    {
        FESystem::Mesh::Tri6::tri6_nondegenerate_to_degenerate_element_mapping.reset(new FESystem::Numerics::DenseMatrix<FESystemDouble>);
        FESystem::Mesh::Tri6::tri6_nondegenerate_to_degenerate_element_mapping->resize(6,9);

        FESystem::Mesh::Tri6::tri6_nondegenerate_to_degenerate_element_mapping->setVal(0, 0,  1.0);
        FESystem::Mesh::Tri6::tri6_nondegenerate_to_degenerate_element_mapping->setVal(0, 8,-0.125);

        FESystem::Mesh::Tri6::tri6_nondegenerate_to_degenerate_element_mapping->setVal(1, 1,  1.0);
        FESystem::Mesh::Tri6::tri6_nondegenerate_to_degenerate_element_mapping->setVal(1, 8,-0.125);

        // this does not have any contribution from the 9th node of QUAD9, because the corner nodes (#3, #4) give a -.025 contribution each, 
        // and the middle node (#7) gives a 0.5 contribution
        FESystem::Mesh::Tri6::tri6_nondegenerate_to_degenerate_element_mapping->setVal(2, 2,  1.0); // last node is repeated
        FESystem::Mesh::Tri6::tri6_nondegenerate_to_degenerate_element_mapping->setVal(2, 3,  1.0); // last node is repeated
        FESystem::Mesh::Tri6::tri6_nondegenerate_to_degenerate_element_mapping->setVal(2, 6,  1.0); // last node is repeated
        
        FESystem::Mesh::Tri6::tri6_nondegenerate_to_degenerate_element_mapping->setVal(3, 4,  1.0);
        FESystem::Mesh::Tri6::tri6_nondegenerate_to_degenerate_element_mapping->setVal(3, 8,  0.25);

        FESystem::Mesh::Tri6::tri6_nondegenerate_to_degenerate_element_mapping->setVal(4, 5,  1.0);
        FESystem::Mesh::Tri6::tri6_nondegenerate_to_degenerate_element_mapping->setVal(4, 8,  0.5);

        FESystem::Mesh::Tri6::tri6_nondegenerate_to_degenerate_element_mapping->setVal(5, 7,  1.0); 
        FESystem::Mesh::Tri6::tri6_nondegenerate_to_degenerate_element_mapping->setVal(5, 8,  0.5);
    }
}



FESystem::Mesh::Tri6::~Tri6()
{
    
}


FESystemUInt
FESystem::Mesh::Tri6::getGeometryOrder() const
{
    return 2;
}


const FESystem::Numerics::MatrixBase<FESystemDouble>& 
FESystem::Mesh::Tri6::getParentToDegenerateElemMappingMatrix() const
{
    return  *(FESystem::Mesh::Tri6::tri6_nondegenerate_to_degenerate_element_mapping);
}


void
FESystem::Mesh::Tri6::clearParentNondegenerateElement()
{
    if (this->parent_nondegenerate_elem != NULL)
        delete this->parent_nondegenerate_elem;
    this->parent_nondegenerate_elem = NULL;
    
    if (this->center_node_for_parent_nondegenerate_elem != NULL)
        delete this->center_node_for_parent_nondegenerate_elem;
    this->center_node_for_parent_nondegenerate_elem = NULL;
}


void 
FESystem::Mesh::Tri6::initializeParentNondegenerateElement()
{
    FESystemAssert0(this->parent_nondegenerate_elem == NULL, FESystem::Exception::InvalidState);
    this->parent_nondegenerate_elem = new FESystem::Mesh::Quad9();
    FESystemAssert0(this->center_node_for_parent_nondegenerate_elem == NULL, FESystem::Exception::InvalidState);
    
    // create the center node for the parent element; use the same coordinate system as the first node    
    this->center_node_for_parent_nondegenerate_elem = new FESystem::Mesh::Node(this->getNode(0).getCoordinateSystem());    

    this->parent_nondegenerate_elem = new FESystem::Mesh::Quad9();
    this->parent_nondegenerate_elem->setNode(0, this->getNode(0));
    this->parent_nondegenerate_elem->setNode(1, this->getNode(1));
    this->parent_nondegenerate_elem->setNode(2, this->getNode(2));
    this->parent_nondegenerate_elem->setNode(3, this->getNode(2));
    
    this->parent_nondegenerate_elem->setNode(4, this->getNode(3));
    this->parent_nondegenerate_elem->setNode(5, this->getNode(4));
    this->parent_nondegenerate_elem->setNode(6, this->getNode(2)); // same as the top node
    this->parent_nondegenerate_elem->setNode(7, this->getNode(5));
    
    for (FESystemUInt i=0; i<8; i++)
        this->center_node_for_parent_nondegenerate_elem->add(1.0/8.0, this->parent_nondegenerate_elem->getNode(i));

    
    this->parent_nondegenerate_elem->setNode(8, *(this->center_node_for_parent_nondegenerate_elem));
}


