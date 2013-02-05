
/*
 *  Tri3.cpp
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/27/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

// FESystem includes
#include "Mesh/Tri3.h"
#include "Mesh/Quad4.h"
#include "Mesh/Node.h"
#include "Geom/Point.h"
#include "Numerics/DenseMatrix.h"


// initialization of the pointer here
std::auto_ptr<FESystem::Numerics::MatrixBase<FESystemDouble> > FESystem::Mesh::Tri3::tri3_nondegenerate_to_degenerate_element_mapping(NULL);


FESystem::Mesh::Tri3::Tri3(FESystemBoolean local_cs_same_as_global):
FESystem::Mesh::TriElemBase(3, FESystem::Mesh::TRI3, local_cs_same_as_global)
{    
    // initialize the mapping matrix
    if (FESystem::Mesh::Tri3::tri3_nondegenerate_to_degenerate_element_mapping.get() == NULL)
    {
        FESystem::Mesh::Tri3::tri3_nondegenerate_to_degenerate_element_mapping.reset(new FESystem::Numerics::DenseMatrix<FESystemDouble>);
        FESystem::Mesh::Tri3::tri3_nondegenerate_to_degenerate_element_mapping->resize(3,4);
        FESystem::Mesh::Tri3::tri3_nondegenerate_to_degenerate_element_mapping->setVal(0, 0, 1.0);
        FESystem::Mesh::Tri3::tri3_nondegenerate_to_degenerate_element_mapping->setVal(1, 1, 1.0);
        FESystem::Mesh::Tri3::tri3_nondegenerate_to_degenerate_element_mapping->setVal(2, 2, 1.0); // last node is repeated
        FESystem::Mesh::Tri3::tri3_nondegenerate_to_degenerate_element_mapping->setVal(2, 3, 1.0); // last node is repeated
    }
}



FESystem::Mesh::Tri3::~Tri3()
{
    
}


FESystemUInt
FESystem::Mesh::Tri3::getGeometryOrder() const
{
    return 1;
}


const FESystem::Numerics::MatrixBase<FESystemDouble>& 
FESystem::Mesh::Tri3::getParentToDegenerateElemMappingMatrix() const
{
    return  *(FESystem::Mesh::Tri3::tri3_nondegenerate_to_degenerate_element_mapping);
}



void
FESystem::Mesh::Tri3::clearParentNondegenerateElement()
{
    if (this->parent_nondegenerate_elem != NULL)
        delete this->parent_nondegenerate_elem;
    this->parent_nondegenerate_elem = NULL;
}


void 
FESystem::Mesh::Tri3::initializeParentNondegenerateElement()
{
    FESystemAssert0(this->parent_nondegenerate_elem == NULL, FESystem::Exception::InvalidState);
    this->parent_nondegenerate_elem = new FESystem::Mesh::Quad4(this->if_local_physical_cs_same_as_global);
    this->parent_nondegenerate_elem->setToParentNondegenerateElement();
    this->parent_nondegenerate_elem->setNode(0, this->getNode(0));
    this->parent_nondegenerate_elem->setNode(1, this->getNode(1));
    this->parent_nondegenerate_elem->setNode(2, this->getNode(2));
    this->parent_nondegenerate_elem->setNode(3, this->getNode(2));
}

