
/*
 *  Tet4.cpp
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/27/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

// FESystem includes
#include "Mesh/Tet4.h"
#include "Mesh/Hex8.h"
#include "Mesh/Node.h"
#include "Geom/Point.h"
#include "Numerics/DenseMatrix.h"


// initialization of the pointer here
std::auto_ptr<FESystem::Numerics::MatrixBase<FESystemDouble> > FESystem::Mesh::Tet4::tet4_nondegenerate_to_degenerate_element_mapping(NULL);


FESystem::Mesh::Tet4::Tet4(FESystemBoolean local_cs_same_as_global):
FESystem::Mesh::TetElemBase(4, FESystem::Mesh::TET4, local_cs_same_as_global)
{
    // initialize the mapping matrix
    if (FESystem::Mesh::Tet4::tet4_nondegenerate_to_degenerate_element_mapping.get() == NULL)
    {
        FESystem::Mesh::Tet4::tet4_nondegenerate_to_degenerate_element_mapping.reset(new FESystem::Numerics::DenseMatrix<FESystemDouble>);
        FESystem::Mesh::Tet4::tet4_nondegenerate_to_degenerate_element_mapping->resize(4,8);
        FESystem::Mesh::Tet4::tet4_nondegenerate_to_degenerate_element_mapping->setVal(0, 0, 1.0);
        FESystem::Mesh::Tet4::tet4_nondegenerate_to_degenerate_element_mapping->setVal(1, 1, 1.0);
        FESystem::Mesh::Tet4::tet4_nondegenerate_to_degenerate_element_mapping->setVal(2, 2, 1.0); // top node is condensed from top 4 nodes of the hex8
        FESystem::Mesh::Tet4::tet4_nondegenerate_to_degenerate_element_mapping->setVal(2, 3, 1.0);
        FESystem::Mesh::Tet4::tet4_nondegenerate_to_degenerate_element_mapping->setVal(2, 6, 1.0);
        FESystem::Mesh::Tet4::tet4_nondegenerate_to_degenerate_element_mapping->setVal(2, 7, 1.0);
        FESystem::Mesh::Tet4::tet4_nondegenerate_to_degenerate_element_mapping->setVal(3, 4, 1.0); // front node is condensed from the bottom front-two nodes
        FESystem::Mesh::Tet4::tet4_nondegenerate_to_degenerate_element_mapping->setVal(3, 5, 1.0);
    }
}



FESystem::Mesh::Tet4::~Tet4()
{
    
}


FESystemUInt
FESystem::Mesh::Tet4::getGeometryOrder() const
{
    return 1;
}


const FESystem::Numerics::MatrixBase<FESystemDouble>&
FESystem::Mesh::Tet4::getParentToDegenerateElemMappingMatrix() const
{
    return  *(FESystem::Mesh::Tet4::tet4_nondegenerate_to_degenerate_element_mapping);
}



void
FESystem::Mesh::Tet4::clearParentNondegenerateElement()
{
    if (this->parent_nondegenerate_elem != NULL)
        delete this->parent_nondegenerate_elem;
    this->parent_nondegenerate_elem = NULL;
}


void
FESystem::Mesh::Tet4::initializeParentNondegenerateElement()
{
    FESystemAssert0(this->parent_nondegenerate_elem == NULL, FESystem::Exception::InvalidState);
    this->parent_nondegenerate_elem = new FESystem::Mesh::Hex8(this->if_local_physical_cs_same_as_global);
    this->parent_nondegenerate_elem->setToParentNondegenerateElement();
    this->parent_nondegenerate_elem->setNode(0, this->getNode(0));
    this->parent_nondegenerate_elem->setNode(1, this->getNode(1));
    this->parent_nondegenerate_elem->setNode(2, this->getNode(2));
    this->parent_nondegenerate_elem->setNode(3, this->getNode(2));
    this->parent_nondegenerate_elem->setNode(4, this->getNode(3));
    this->parent_nondegenerate_elem->setNode(5, this->getNode(3));
    this->parent_nondegenerate_elem->setNode(6, this->getNode(2));
    this->parent_nondegenerate_elem->setNode(7, this->getNode(2));
}


