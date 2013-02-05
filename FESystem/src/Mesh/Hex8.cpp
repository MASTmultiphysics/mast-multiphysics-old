
/*
 *  Hex8.cpp
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/27/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

#include "Mesh/Hex8.h"
#include "Geom/Point.h"
#include "Geom/RectangularCoordinateSystem.h"
#include "Mesh/Node.h"
#include "Numerics/LocalVector.h"
#include "Numerics/DenseMatrix.h"
#include "Utils/Table.h"


// initialization of the pointer here
std::auto_ptr<FESystem::Utility::Table<FESystemUInt> > FESystem::Mesh::Hex8::hex8_node_dim_ID_table(NULL);
std::auto_ptr<FESystem::Numerics::MatrixBase<FESystemDouble> > FESystem::Mesh::Hex8::hex8_nondegenerate_to_degenerate_element_mapping(NULL);


FESystem::Mesh::Hex8::Hex8(FESystemBoolean local_cs_same_as_global):
FESystem::Mesh::HexElemBase(8, FESystem::Mesh::HEX8, local_cs_same_as_global)
{
    // if the table has not been initialized, do that here.
    if (FESystem::Mesh::Hex8::hex8_node_dim_ID_table.get() == NULL)
    {
        FESystem::Mesh::Hex8::hex8_node_dim_ID_table.reset(new FESystem::Utility::Table<FESystemUInt>);
        // now initialize the data
        std::vector<FESystemUInt> el(2);
        el[0] = this->getNNodes();
        el[1] = this->getDimension();
        FESystemUInt id = 0;
        FESystem::Mesh::Hex8::hex8_node_dim_ID_table->reinit(el);
        for (FESystemUInt i_node=0; i_node<this->getNNodes(); i_node++)
            for (FESystemUInt i_dim=0; i_dim<this->getDimension(); i_dim++)
            {
                el[0] = i_node;
                el[1] = i_dim;
                id = this->getLocalNodeIDALongDim(i_node, i_dim);
                FESystem::Mesh::Hex8::hex8_node_dim_ID_table->setVal(el, id);
            }
    }
    
    // initialize the mapping matrix
    if (FESystem::Mesh::Hex8::hex8_nondegenerate_to_degenerate_element_mapping.get() == NULL)
    {
        FESystem::Mesh::Hex8::hex8_nondegenerate_to_degenerate_element_mapping.reset(new FESystem::Numerics::DenseMatrix<FESystemDouble>);
        FESystem::Mesh::Hex8::hex8_nondegenerate_to_degenerate_element_mapping->resize(8,8);
        FESystem::Mesh::Hex8::hex8_nondegenerate_to_degenerate_element_mapping->setToIdentity();
    }
}



FESystem::Mesh::Hex8::~Hex8()
{
    
}



FESystemUInt
FESystem::Mesh::Hex8::getGeometryOrder() const
{
    return 1;
}



const FESystem::Utility::Table<FESystemUInt>&
FESystem::Mesh::Hex8::getPointDimensionIDTable() const
{
    return  *(FESystem::Mesh::Hex8::hex8_node_dim_ID_table);
}


const FESystem::Numerics::MatrixBase<FESystemDouble>&
FESystem::Mesh::Hex8::getParentToDegenerateElemMappingMatrix() const
{
    return  *(FESystem::Mesh::Hex8::hex8_nondegenerate_to_degenerate_element_mapping);
}



bool
FESystem::Mesh::Hex8::ifDegerateElement() const
{
    return false;
}



const FESystem::Mesh::ElemBase&
FESystem::Mesh::Hex8::getParentNondegenerateElem() const
{
    return *this;
}



void
FESystem::Mesh::Hex8::clearParentNondegenerateElement()
{
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}


void
FESystem::Mesh::Hex8::initializeParentNondegenerateElement()
{
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}



FESystemUInt
FESystem::Mesh::Hex8::getNNodesAlongDimension(FESystemUInt n) const
{
    FESystemUInt nnodes = 0;
    
    switch (n) {
        case 0:
        case 1:
        case 2:
            nnodes = 2;
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
    }
    
    return nnodes;
}



FESystemUInt
FESystem::Mesh::Hex8::getLocalNodeIDALongDim(FESystemUInt n_id, FESystemUInt d) const
{
    FESystemUInt dim_node_num = 0;
    
    switch (n_id) {
        case 0:
            switch (d) {
                case 0:
                case 1:
                case 2:
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
                case 2:
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
                    
                case 2:
                    dim_node_num = 0;
                    break;

                default:
                    FESystemAssert0(false, FESystem::Exception::InvalidValue);
                    break;
            }
            break;
            
        case 3:
            switch (d) {
                case 0:
                case 2:
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
            
        case 4:
            switch (d) {
                case 0:
                case 1:
                    dim_node_num = 0;
                    break;
                    
                case 2:
                    dim_node_num = 1;
                    break;

                default:
                    FESystemAssert0(false, FESystem::Exception::InvalidValue);
                    break;
            }
            break;
            
        case 5:
            switch (d) {
                case 0:
                case 2:
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
            
        case 6:
            switch (d) {
                case 0:
                case 1:
                case 2:
                    dim_node_num = 1;
                    break;
                    
                default:
                    FESystemAssert0(false, FESystem::Exception::InvalidValue);
                    break;
            }
            break;
            
        case 7:
            switch (d) {
                case 0:
                    dim_node_num = 0;
                    break;
                    
                case 1:
                case 2:
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
FESystem::Mesh::Hex8::getLocalComputationalCoordinateForLocalNodeID(FESystemUInt id, FESystem::Geometry::Point& p) const
{
    switch (id) {
        case 0:
            p.setVal(0, -1.0);
            p.setVal(1, -1.0);
            p.setVal(2, -1.0);
            break;
            
        case 1:
            p.setVal(0,  1.0);
            p.setVal(1, -1.0);
            p.setVal(2, -1.0);
            break;
            
        case 2:
            p.setVal(0,  1.0);
            p.setVal(1,  1.0);
            p.setVal(2, -1.0);
            break;
            
        case 3:
            p.setVal(0, -1.0);
            p.setVal(1,  1.0);
            p.setVal(2, -1.0);
            break;

        case 4:
            p.setVal(0, -1.0);
            p.setVal(1, -1.0);
            p.setVal(2,  1.0);
            break;
            
        case 5:
            p.setVal(0,  1.0);
            p.setVal(1, -1.0);
            p.setVal(2,  1.0);
            break;
            
        case 6:
            p.setVal(0,  1.0);
            p.setVal(1,  1.0);
            p.setVal(2,  1.0);
            break;
            
        case 7:
            p.setVal(0, -1.0);
            p.setVal(1,  1.0);
            p.setVal(2,  1.0);
            break;

        default:
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
    }
}


