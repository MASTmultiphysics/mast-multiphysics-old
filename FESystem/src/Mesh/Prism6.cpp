//
//  Prism6.cpp
//  FESystem
//
//  Created by Manav Bhatia on 1/24/13.
//
//

#include "Mesh/Prism6.h"
#include "Mesh/Hex8.h"
#include "Geom/Point.h"
#include "Geom/RectangularCoordinateSystem.h"
#include "Mesh/Node.h"
#include "Numerics/LocalVector.h"
#include "Numerics/DenseMatrix.h"
#include "Utils/Table.h"


// initialization of the pointer here
std::auto_ptr<FESystem::Numerics::MatrixBase<FESystemDouble> > FESystem::Mesh::Prism6::prism6_nondegenerate_to_degenerate_element_mapping(NULL);


FESystem::Mesh::Prism6::Prism6(FESystemBoolean local_cs_same_as_global):
FESystem::Mesh::PrismElemBase(6, FESystem::Mesh::PRISM6, local_cs_same_as_global)
{
    // initialize the mapping matrix
    if (FESystem::Mesh::Prism6::prism6_nondegenerate_to_degenerate_element_mapping.get() == NULL)
    {
        FESystem::Mesh::Prism6::prism6_nondegenerate_to_degenerate_element_mapping.reset(new FESystem::Numerics::DenseMatrix<FESystemDouble>);
        FESystem::Mesh::Prism6::prism6_nondegenerate_to_degenerate_element_mapping->resize(6,8);
        FESystem::Mesh::Prism6::prism6_nondegenerate_to_degenerate_element_mapping->setVal(0, 0, 1.0);
        FESystem::Mesh::Prism6::prism6_nondegenerate_to_degenerate_element_mapping->setVal(1, 1, 1.0);
        FESystem::Mesh::Prism6::prism6_nondegenerate_to_degenerate_element_mapping->setVal(2, 2, 1.0);
        FESystem::Mesh::Prism6::prism6_nondegenerate_to_degenerate_element_mapping->setVal(3, 3, 1.0);
        FESystem::Mesh::Prism6::prism6_nondegenerate_to_degenerate_element_mapping->setVal(4, 4, 1.0); // front node is condensed from the front nodes
        FESystem::Mesh::Prism6::prism6_nondegenerate_to_degenerate_element_mapping->setVal(4, 5, 1.0);
        FESystem::Mesh::Prism6::prism6_nondegenerate_to_degenerate_element_mapping->setVal(5, 6, 1.0);
        FESystem::Mesh::Prism6::prism6_nondegenerate_to_degenerate_element_mapping->setVal(5, 7, 1.0);
    }
}



FESystem::Mesh::Prism6::~Prism6()
{
    
}


FESystemUInt
FESystem::Mesh::Prism6::getGeometryOrder() const
{
    return 1;
}


const FESystem::Numerics::MatrixBase<FESystemDouble>&
FESystem::Mesh::Prism6::getParentToDegenerateElemMappingMatrix() const
{
    return  *(FESystem::Mesh::Prism6::prism6_nondegenerate_to_degenerate_element_mapping);
}



void
FESystem::Mesh::Prism6::clearParentNondegenerateElement()
{
    if (this->parent_nondegenerate_elem != NULL)
        delete this->parent_nondegenerate_elem;
    this->parent_nondegenerate_elem = NULL;
}


void
FESystem::Mesh::Prism6::initializeParentNondegenerateElement()
{
    FESystemAssert0(this->parent_nondegenerate_elem == NULL, FESystem::Exception::InvalidState);
    this->parent_nondegenerate_elem = new FESystem::Mesh::Hex8(this->if_local_physical_cs_same_as_global);
    this->parent_nondegenerate_elem->setToParentNondegenerateElement();
    this->parent_nondegenerate_elem->setNode(0, this->getNode(0));
    this->parent_nondegenerate_elem->setNode(1, this->getNode(1));
    this->parent_nondegenerate_elem->setNode(2, this->getNode(2));
    this->parent_nondegenerate_elem->setNode(3, this->getNode(3));
    this->parent_nondegenerate_elem->setNode(4, this->getNode(4));
    this->parent_nondegenerate_elem->setNode(5, this->getNode(4));
    this->parent_nondegenerate_elem->setNode(6, this->getNode(5));
    this->parent_nondegenerate_elem->setNode(7, this->getNode(5));
}


