//
//  Quad8.cpp
//  FESystem
//
//  Created by Manav Bhatia on 5/21/12.
//  Copyright (c) 2012. All rights reserved.
//

#include "Mesh/Quad8.h"
#include "Mesh/Quad9.h"
#include "Mesh/Node.h"
#include "Geom/Point.h"
#include "Geom/RectangularCoordinateSystem.h"
#include "Numerics/LocalVector.h"
#include "Numerics/DenseMatrix.h"


// initialization of the pointer here
std::auto_ptr<FESystem::Numerics::MatrixBase<FESystemDouble> > FESystem::Mesh::Quad8::quad8_nondegenerate_to_degenerate_element_mapping(NULL);


FESystem::Mesh::Quad8::Quad8(FESystemBoolean local_cs_same_as_global):
FESystem::Mesh::QuadElemBase(8, FESystem::Mesh::QUAD8, local_cs_same_as_global),
center_node_for_parent_nondegenerate_elem(NULL)
{
    // initialize the mapping matrix
    if (FESystem::Mesh::Quad8::quad8_nondegenerate_to_degenerate_element_mapping.get() == NULL)
    {
        FESystem::Mesh::Quad8::quad8_nondegenerate_to_degenerate_element_mapping.reset(new FESystem::Numerics::DenseMatrix<FESystemDouble>);
        FESystem::Mesh::Quad8::quad8_nondegenerate_to_degenerate_element_mapping->resize(8,9);
        // corner nodes
        for (FESystemUInt i=0; i<4; i++)
        {
            FESystem::Mesh::Quad8::quad8_nondegenerate_to_degenerate_element_mapping->setVal(i,i,  1.0);
            FESystem::Mesh::Quad8::quad8_nondegenerate_to_degenerate_element_mapping->setVal(i,8,-0.25);
        }
        // midside nodes
        for (FESystemUInt i=0; i<4; i++)
        {
            FESystem::Mesh::Quad8::quad8_nondegenerate_to_degenerate_element_mapping->setVal(i+4,i+4,  1.0);
            FESystem::Mesh::Quad8::quad8_nondegenerate_to_degenerate_element_mapping->setVal(i+4,  8,  0.5);
        }
    }
}



FESystem::Mesh::Quad8::~Quad8()
{
    if (this->center_node_for_parent_nondegenerate_elem != NULL)
        delete this->center_node_for_parent_nondegenerate_elem;
}



FESystemUInt
FESystem::Mesh::Quad8::getGeometryOrder() const
{
    return 2;
}



const FESystem::Utility::Table<FESystemUInt>& 
FESystem::Mesh::Quad8::getPointDimensionIDTable() const
{
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}




const FESystem::Numerics::MatrixBase<FESystemDouble>& 
FESystem::Mesh::Quad8::getParentToDegenerateElemMappingMatrix() const
{
    return  *(FESystem::Mesh::Quad8::quad8_nondegenerate_to_degenerate_element_mapping);
}


bool
FESystem::Mesh::Quad8::ifDegerateElement() const
{
    return true;
}



const FESystem::Mesh::ElemBase& 
FESystem::Mesh::Quad8::getParentNondegenerateElem() const
{
    return *(this->parent_nondegenerate_elem);
}


void
FESystem::Mesh::Quad8::clearParentNondegenerateElement()
{
    if (this->parent_nondegenerate_elem != NULL)
        delete this->parent_nondegenerate_elem;
    this->parent_nondegenerate_elem = NULL;
    
    if (this->center_node_for_parent_nondegenerate_elem != NULL)
        delete this->center_node_for_parent_nondegenerate_elem;
    this->center_node_for_parent_nondegenerate_elem = NULL;
}


void
FESystem::Mesh::Quad8::initializeParentNondegenerateElement()
{
    FESystemAssert0(this->parent_nondegenerate_elem == NULL, FESystem::Exception::InvalidState);
    this->parent_nondegenerate_elem = new FESystem::Mesh::Quad9(this->if_local_physical_cs_same_as_global);
    FESystemAssert0(this->center_node_for_parent_nondegenerate_elem == NULL, FESystem::Exception::InvalidState);
 
    this->parent_nondegenerate_elem->setToParentNondegenerateElement();

    // create the center node for the parent element; use the same coordinate system as the first node
    this->center_node_for_parent_nondegenerate_elem = new FESystem::Mesh::Node(this->getNode(0).getCoordinateSystem()); 
    for (FESystemUInt i=0; i<8; i++)
        this->center_node_for_parent_nondegenerate_elem->add(1.0/8.0, this->getNode(i));
    
    for (FESystemUInt i=0; i<8; i++)
        this->parent_nondegenerate_elem->setNode(i, this->getNode(i));    
    this->parent_nondegenerate_elem->setNode(8, *(this->center_node_for_parent_nondegenerate_elem));
}



FESystemUInt 
FESystem::Mesh::Quad8::getNNodesAlongDimension(FESystemUInt n) const
{
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}



FESystemUInt 
FESystem::Mesh::Quad8::getLocalNodeIDALongDim(FESystemUInt n_id, FESystemUInt d) const
{
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}




void
FESystem::Mesh::Quad8::getLocalComputationalCoordinateForLocalNodeID(FESystemUInt id, FESystem::Geometry::Point& p) const
{
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}



