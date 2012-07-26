//
//  Quad9.cpp
//  FESystem
//
//  Created by Manav Bhatia on 5/17/12.
//  Copyright (c) 2012 . All rights reserved.
//

#include "Mesh/Quad9.h"
#include "Geom/Point.h"
#include "Geom/RectangularCoordinateSystem.h"
#include "Mesh/Node.h"
#include "Numerics/LocalVector.h"
#include "Numerics/DenseMatrix.h"


// initialization of the pointer here
std::auto_ptr<FESystem::Utility::Table<FESystemUInt> > FESystem::Mesh::Quad9::quad9_node_dim_ID_table(NULL);
std::auto_ptr<FESystem::Numerics::MatrixBase<FESystemDouble> > FESystem::Mesh::Quad9::quad9_nondegenerate_to_degenerate_element_mapping(NULL);


FESystem::Mesh::Quad9::Quad9():
FESystem::Mesh::QuadElemBase(9, FESystem::Mesh::QUAD9)
{
    // if the table has not been initialized, do that here. 
    if (FESystem::Mesh::Quad9::quad9_node_dim_ID_table.get() == NULL)
    {
        FESystem::Mesh::Quad9::quad9_node_dim_ID_table.reset(new FESystem::Utility::Table<FESystemUInt>);
        // now initialize the data
        std::vector<FESystemUInt> el(2);
        el[0] = this->getNNodes();
        el[1] = this->getDimension();
        FESystemUInt id = 0;
        FESystem::Mesh::Quad9::quad9_node_dim_ID_table->reinit(el);
        for (FESystemUInt i_node=0; i_node<this->getNNodes(); i_node++)
            for (FESystemUInt i_dim=0; i_dim<this->getDimension(); i_dim++)
            {
                el[0] = i_node;
                el[1] = i_dim;
                id = this->getLocalNodeIDALongDim(i_node, i_dim);
                FESystem::Mesh::Quad9::quad9_node_dim_ID_table->setVal(el, id);
            }
    }
    
    // initialize the mapping matrix
    if (FESystem::Mesh::Quad9::quad9_nondegenerate_to_degenerate_element_mapping.get() == NULL)
    {
        FESystem::Mesh::Quad9::quad9_nondegenerate_to_degenerate_element_mapping.reset(new FESystem::Numerics::DenseMatrix<FESystemDouble>);
        FESystem::Mesh::Quad9::quad9_nondegenerate_to_degenerate_element_mapping->resize(9,9);
        FESystem::Mesh::Quad9::quad9_nondegenerate_to_degenerate_element_mapping->setToIdentity();
    }
}



FESystem::Mesh::Quad9::~Quad9()
{
    
}


FESystemUInt
FESystem::Mesh::Quad9::getGeometryOrder() const
{
    return 2;
}


const FESystem::Utility::Table<FESystemUInt>& 
FESystem::Mesh::Quad9::getPointDimensionIDTable() const
{
    return  *(FESystem::Mesh::Quad9::quad9_node_dim_ID_table);
}




const FESystem::Numerics::MatrixBase<FESystemDouble>& 
FESystem::Mesh::Quad9::getParentToDegenerateElemMappingMatrix() const
{
    return  *(FESystem::Mesh::Quad9::quad9_nondegenerate_to_degenerate_element_mapping);
}


bool
FESystem::Mesh::Quad9::ifDegerateElement() const
{
    return false;
}



const FESystem::Mesh::ElemBase& 
FESystem::Mesh::Quad9::getParentNondegenerateElem() const
{
    return *this;
}



void
FESystem::Mesh::Quad9::initializeParentNondegenerateElement()
{
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}



FESystemUInt 
FESystem::Mesh::Quad9::getNNodesAlongDimension(FESystemUInt n) const
{
    FESystemUInt nnodes = 0;
    
    switch (n) {
        case 0:
        case 1:
            nnodes = 3;
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
    }
    
    return nnodes;
}



FESystemUInt 
FESystem::Mesh::Quad9::getLocalNodeIDALongDim(FESystemUInt n_id, FESystemUInt d) const
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
            
        case 4:
            switch (d) {
                case 0:
                    dim_node_num = 2;
                    break;

                case 1:
                    dim_node_num = 0;
                    break;
                    
                default:
                    FESystemAssert0(false, FESystem::Exception::InvalidValue);
                    break;
            }
            break;
            
        case 5:
            switch (d) {
                case 0:
                    dim_node_num = 1;
                    break;
                    
                case 1:
                    dim_node_num = 2;
                    break;
                    
                default:
                    FESystemAssert0(false, FESystem::Exception::InvalidValue);
                    break;
            }
            break;
            
        case 6:
            switch (d) {
                case 0:
                    dim_node_num = 2;
                    break;

                case 1:
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
                    dim_node_num = 2;
                    break;
                    
                default:
                    FESystemAssert0(false, FESystem::Exception::InvalidValue);
                    break;
            }
            break;

        case 8:
            switch (d) {
                case 0:
                    dim_node_num = 2;
                    break;
                    
                case 1:
                    dim_node_num = 2;
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
FESystem::Mesh::Quad9::getLocalComputationalCoordinateForLocalNodeID(FESystemUInt id, FESystem::Geometry::Point& p) const
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
            
        case 4:
            p.setVal(0,  0.0); 
            p.setVal(1, -1.0); 
            break;
            
        case 5:
            p.setVal(0,  1.0); 
            p.setVal(1,  0.0); 
            break;
            
        case 6:
            p.setVal(0,  0.0); 
            p.setVal(1,  1.0); 
            break;
            
        case 7:
            p.setVal(0, -1.0); 
            p.setVal(1,  0.0); 
            break;

        case 8:
            p.setVal(0,  0.0); 
            p.setVal(1,  0.0);
            break;

        default:
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
    }
}



