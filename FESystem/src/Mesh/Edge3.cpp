//
//  Edge3.cpp
//  FESystem
//
//  Created by Manav Bhatia on 5/17/12.
//  Copyright (c) 2012. All rights reserved.
//

// FESystem includes
#include "Mesh/Edge3.h"
#include "Utils/Table.h"
#include "Mesh/Node.h"
#include "Geom/RectangularCoordinateSystem.h"
#include "Numerics/DenseMatrix.h"


// initialization of the pointer here
std::auto_ptr<FESystem::Utility::Table<FESystemUInt> > FESystem::Mesh::Edge3::edge3_node_dim_ID_table(NULL);
std::auto_ptr<FESystem::Numerics::MatrixBase<FESystemDouble> > FESystem::Mesh::Edge3::edge3_nondegenerate_to_degenerate_element_mapping(NULL);



FESystem::Mesh::Edge3::Edge3(FESystemBoolean local_cs_same_as_global):
FESystem::Mesh::EdgeElemBase(3, FESystem::Mesh::EDGE3, local_cs_same_as_global) // specify that two nodes are used to define this element
{
    this->y_unit_vec = new FESystem::Numerics::LocalVector<FESystemDouble>;
    this->y_unit_vec->resize(3);
    
    // if the table has not been initialized, do that here. 
    if (FESystem::Mesh::Edge3::edge3_node_dim_ID_table.get() == NULL)
    {
        FESystem::Mesh::Edge3::edge3_node_dim_ID_table.reset(new FESystem::Utility::Table<FESystemUInt>);
        // now initialize the data
        std::vector<FESystemUInt> el(2);
        el[0] = this->getNNodes();
        el[1] = this->getDimension();
        FESystemUInt id = 0;
        FESystem::Mesh::Edge3::edge3_node_dim_ID_table->reinit(el);
        for (FESystemUInt i_node=0; i_node<this->getNNodes(); i_node++)
            for (FESystemUInt i_dim=0; i_dim<this->getDimension(); i_dim++)
            {
                el[0] = i_node;
                el[1] = i_dim;
                id = this->getLocalNodeIDALongDim(i_node, i_dim);
                FESystem::Mesh::Edge3::edge3_node_dim_ID_table->setVal(el, id);
            }
    }
    
    // initialize the mapping matrix
    if (FESystem::Mesh::Edge3::edge3_nondegenerate_to_degenerate_element_mapping.get() == NULL)
    {
        FESystem::Mesh::Edge3::edge3_nondegenerate_to_degenerate_element_mapping.reset(new FESystem::Numerics::DenseMatrix<FESystemDouble>);
        FESystem::Mesh::Edge3::edge3_nondegenerate_to_degenerate_element_mapping->resize(3,3);
        FESystem::Mesh::Edge3::edge3_nondegenerate_to_degenerate_element_mapping->setToIdentity();
    }
}



FESystem::Mesh::Edge3::~Edge3()
{

}


FESystemUInt
FESystem::Mesh::Edge3::getGeometryOrder() const
{
    return 2;
}


void 
FESystem::Mesh::Edge3::setVectorForXYPlane(const FESystem::Numerics::VectorBase<FESystemDouble>& vec)
{
    // make sure that the dimensions are appropriate
    FESystemAssert2(vec.getSize() == this->y_unit_vec->getSize(), FESystem::Exception::DimensionsDoNotMatch, vec.getSize(), this->y_unit_vec->getSize());
    
    // make sure that vec has a finite length
    FESystemAssert0(vec.getL2Norm() > 0, FESystem::Exception::InvalidValue);
    
    this->y_unit_vec->copyVector(vec);
}



const FESystem::Utility::Table<FESystemUInt>& 
FESystem::Mesh::Edge3::getPointDimensionIDTable() const
{
    return  *(FESystem::Mesh::Edge3::edge3_node_dim_ID_table);
}


const FESystem::Numerics::MatrixBase<FESystemDouble>& 
FESystem::Mesh::Edge3::getParentToDegenerateElemMappingMatrix() const
{
    return  *(FESystem::Mesh::Edge3::edge3_nondegenerate_to_degenerate_element_mapping);
}


FESystemUInt 
FESystem::Mesh::Edge3::getNNodesAlongDimension(FESystemUInt n) const
{
    FESystemUInt nnodes = 0;
    
    switch (n) {
        case 0:
            nnodes = 3;
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
    }
    
    return nnodes;
}



FESystemUInt 
FESystem::Mesh::Edge3::getLocalNodeIDALongDim(FESystemUInt n_id, FESystemUInt d) const
{
    FESystemUInt id = 0;
    
    switch (d) {
        case 0:
            switch (n_id) {
                case 0:
                    id = 0;
                    break;
                    
                case 1:
                    id = 1;
                    break;
                    
                case 2:
                    id = 2;
                    break;

                default:
                    FESystemAssert0(false, FESystem::Exception::InvalidValue);
                    break;
            }
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
    }
    
    return id;
}




void
FESystem::Mesh::Edge3::getLocalComputationalCoordinateForLocalNodeID(FESystemUInt id, FESystem::Geometry::Point& p) const
{
    switch (id) {
        case 0:
            p.setVal(0, -1.0); // left edge node
            break;
            
        case 1:
            p.setVal(0,  1.0); // right edge node
            break;
            
        case 2:
            p.setVal(0,  0.0); // right edge node
            break;

        default:
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
    }
}


