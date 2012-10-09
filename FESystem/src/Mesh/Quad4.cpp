
/*
 *  Quad4.cpp
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/27/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

#include "Mesh/Quad4.h"
#include "Geom/Point.h"
#include "Geom/RectangularCoordinateSystem.h"
#include "Mesh/Node.h"
#include "Numerics/LocalVector.h"
#include "Numerics/DenseMatrix.h"


// initialization of the pointer here
std::auto_ptr<FESystem::Utility::Table<FESystemUInt> > FESystem::Mesh::Quad4::quad4_node_dim_ID_table(NULL);
std::auto_ptr<FESystem::Numerics::MatrixBase<FESystemDouble> > FESystem::Mesh::Quad4::quad4_nondegenerate_to_degenerate_element_mapping(NULL);


FESystem::Mesh::Quad4::Quad4(FESystemBoolean local_cs_same_as_global):
FESystem::Mesh::QuadElemBase(4, FESystem::Mesh::QUAD4, local_cs_same_as_global)
{
    // if the table has not been initialized, do that here. 
    if (FESystem::Mesh::Quad4::quad4_node_dim_ID_table.get() == NULL)
    {
        FESystem::Mesh::Quad4::quad4_node_dim_ID_table.reset(new FESystem::Utility::Table<FESystemUInt>);
        // now initialize the data
        std::vector<FESystemUInt> el(2);
        el[0] = this->getNNodes();
        el[1] = this->getDimension();
        FESystemUInt id = 0;
        FESystem::Mesh::Quad4::quad4_node_dim_ID_table->reinit(el);
        for (FESystemUInt i_node=0; i_node<this->getNNodes(); i_node++)
            for (FESystemUInt i_dim=0; i_dim<this->getDimension(); i_dim++)
            {
                el[0] = i_node;
                el[1] = i_dim;
                id = this->getLocalNodeIDALongDim(i_node, i_dim);
                FESystem::Mesh::Quad4::quad4_node_dim_ID_table->setVal(el, id);
            }
    }
    
    // initialize the mapping matrix
    if (FESystem::Mesh::Quad4::quad4_nondegenerate_to_degenerate_element_mapping.get() == NULL)
    {
        FESystem::Mesh::Quad4::quad4_nondegenerate_to_degenerate_element_mapping.reset(new FESystem::Numerics::DenseMatrix<FESystemDouble>);
        FESystem::Mesh::Quad4::quad4_nondegenerate_to_degenerate_element_mapping->resize(4,4);
        FESystem::Mesh::Quad4::quad4_nondegenerate_to_degenerate_element_mapping->setToIdentity();
    }
}



FESystem::Mesh::Quad4::~Quad4()
{
    
}



FESystemUInt
FESystem::Mesh::Quad4::getGeometryOrder() const
{
    return 1;
}



const FESystem::Utility::Table<FESystemUInt>& 
FESystem::Mesh::Quad4::getPointDimensionIDTable() const
{
    return  *(FESystem::Mesh::Quad4::quad4_node_dim_ID_table);
}


const FESystem::Numerics::MatrixBase<FESystemDouble>& 
FESystem::Mesh::Quad4::getParentToDegenerateElemMappingMatrix() const
{
    return  *(FESystem::Mesh::Quad4::quad4_nondegenerate_to_degenerate_element_mapping);
}



bool
FESystem::Mesh::Quad4::ifDegerateElement() const
{
    return false;
}



const FESystem::Mesh::ElemBase& 
FESystem::Mesh::Quad4::getParentNondegenerateElem() const
{
    return *this;
}



void
FESystem::Mesh::Quad4::clearParentNondegenerateElement()
{
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}


void
FESystem::Mesh::Quad4::initializeParentNondegenerateElement()
{
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}



FESystemUInt 
FESystem::Mesh::Quad4::getNNodesAlongDimension(FESystemUInt n) const
{
    FESystemUInt nnodes = 0;
    
    switch (n) {
        case 0:
        case 1:
            nnodes = 2;
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
    }
    
    return nnodes;
}



FESystemUInt 
FESystem::Mesh::Quad4::getLocalNodeIDALongDim(FESystemUInt n_id, FESystemUInt d) const
{
    FESystemUInt dim_node_num = 0;
    
    switch (n_id) {
        case 0:
            switch (d) {
                case 0:
                case 1:
                    dim_node_num = 0;
                    break;
                    
                default:
                    FESystemAssert0(false, FESystem::Exception::InvalidValue);
                    break;
            }
            break;
            
        case 1:
            switch (d) {
                case 0:
                    dim_node_num = 1;
                    break;

                case 1:
                    dim_node_num = 0;
                    break;
                    
                default:
                    FESystemAssert0(false, FESystem::Exception::InvalidValue);
                    break;
            }
            break;

        case 2:
            switch (d) {
                case 0:
                case 1:
                    dim_node_num = 1;
                    break;
                    
                default:
                    FESystemAssert0(false, FESystem::Exception::InvalidValue);
                    break;
            }
            break;

        case 3:
            switch (d) {
                case 0:
                    dim_node_num = 0;
                    break;

                case 1:
                    dim_node_num = 1;
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
    
    return dim_node_num;
}




void
FESystem::Mesh::Quad4::getLocalComputationalCoordinateForLocalNodeID(FESystemUInt id, FESystem::Geometry::Point& p) const
{
    switch (id) {
        case 0:
            p.setVal(0, -1.0); 
            p.setVal(1, -1.0); 
            break;
            
        case 1:
            p.setVal(0,  1.0); 
            p.setVal(1, -1.0); 
            break;
            
        case 2:
            p.setVal(0,  1.0); 
            p.setVal(1,  1.0); 
            break;
            
        case 3:
            p.setVal(0, -1.0); 
            p.setVal(1,  1.0); 
            break;

        default:
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
    }
}


