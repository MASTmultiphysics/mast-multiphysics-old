//
//  Tri7.cpp
//  FESystem
//
//  Created by Manav Bhatia on 5/17/12.
//  Copyright (c) 2012 . All rights reserved.
//

// FESystem includes
#include "Mesh/Tri7.h"
#include "Mesh/Quad9.h"
#include "Mesh/Node.h"
#include "Geom/Point.h"
#include "Numerics/DenseMatrix.h"


// initialization of the pointer here
std::auto_ptr<FESystem::Numerics::MatrixBase<FESystemDouble> > FESystem::Mesh::Tri7::tri7_nondegenerate_to_degenerate_element_mapping(NULL);


FESystem::Mesh::Tri7::Tri7():
FESystem::Mesh::TriElemBase(7, FESystem::Mesh::TRI7)
{    
    // initialize the mapping matrix
    if (FESystem::Mesh::Tri7::tri7_nondegenerate_to_degenerate_element_mapping.get() == NULL)
    {
        FESystem::Mesh::Tri7::tri7_nondegenerate_to_degenerate_element_mapping.reset(new FESystem::Numerics::DenseMatrix<FESystemDouble>);
        FESystem::Mesh::Tri7::tri7_nondegenerate_to_degenerate_element_mapping->resize(7,9);
        FESystem::Mesh::Tri7::tri7_nondegenerate_to_degenerate_element_mapping->setVal(0, 0, 1.0);
        FESystem::Mesh::Tri7::tri7_nondegenerate_to_degenerate_element_mapping->setVal(1, 1, 1.0);
        FESystem::Mesh::Tri7::tri7_nondegenerate_to_degenerate_element_mapping->setVal(2, 2, 1.0); // last node is repeated
        FESystem::Mesh::Tri7::tri7_nondegenerate_to_degenerate_element_mapping->setVal(2, 3, 1.0); // last node is repeated
        FESystem::Mesh::Tri7::tri7_nondegenerate_to_degenerate_element_mapping->setVal(2, 6, 1.0); // last node is repeated

        FESystem::Mesh::Tri7::tri7_nondegenerate_to_degenerate_element_mapping->setVal(3, 4, 1.0);
        FESystem::Mesh::Tri7::tri7_nondegenerate_to_degenerate_element_mapping->setVal(4, 5, 1.0);
        FESystem::Mesh::Tri7::tri7_nondegenerate_to_degenerate_element_mapping->setVal(5, 7, 1.0); 

        FESystem::Mesh::Tri7::tri7_nondegenerate_to_degenerate_element_mapping->setVal(6, 8, 1.0);
    }
}



FESystem::Mesh::Tri7::~Tri7()
{
    
}



FESystemUInt
FESystem::Mesh::Tri7::getGeometryOrder() const
{
    return 2;
}


const FESystem::Numerics::MatrixBase<FESystemDouble>& 
FESystem::Mesh::Tri7::getParentToDegenerateElemMappingMatrix() const
{
    return  *(FESystem::Mesh::Tri7::tri7_nondegenerate_to_degenerate_element_mapping);
}



void 
FESystem::Mesh::Tri7::initializeParentNondegenerateElement()
{
    FESystemAssert0(this->parent_nondegenerate_elem == NULL, FESystem::Exception::InvalidState);
    this->parent_nondegenerate_elem = new FESystem::Mesh::Quad9();
    this->parent_nondegenerate_elem->setNode(0, this->getNode(0));
    this->parent_nondegenerate_elem->setNode(1, this->getNode(1));
    this->parent_nondegenerate_elem->setNode(2, this->getNode(2));
    this->parent_nondegenerate_elem->setNode(3, this->getNode(2));

    this->parent_nondegenerate_elem->setNode(4, this->getNode(3));
    this->parent_nondegenerate_elem->setNode(5, this->getNode(4));
    this->parent_nondegenerate_elem->setNode(6, this->getNode(2)); // same as the top node
    this->parent_nondegenerate_elem->setNode(7, this->getNode(5));

    this->parent_nondegenerate_elem->setNode(8, this->getNode(6));
}


