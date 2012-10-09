
/*
 *  FaceElemBase.cpp
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/27/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

// FESystem includes
#include "Mesh/FaceElemBase.h"
#include "Mesh/Node.h"
#include "Numerics/DenseMatrix.h"
#include "Numerics/LocalVector.h"
#include "Geom/RectangularCoordinateSystem.h"

// initialization of the pointer here
std::auto_ptr<std::map<FESystemUInt, std::vector<FESystemUInt> > > FESystem::Mesh::QuadElemBase::quad_boundary_node_set(NULL);
std::auto_ptr<std::map<FESystemUInt, std::vector<FESystemUInt> > > FESystem::Mesh::TriElemBase::tri_boundary_node_set(NULL);



FESystem::Mesh::FaceElemBase::FaceElemBase(FESystemUInt nnodes, FESystem::Mesh::ElementType type, FESystemBoolean local_cs_same_as_global):
FESystem::Mesh::ElemBase(nnodes, type, local_cs_same_as_global)
{
    
}


FESystem::Mesh::FaceElemBase::~FaceElemBase()
{
    
}

FESystemUInt 
FESystem::Mesh::FaceElemBase::getDimension() const
{
    return 2;
}



void
FESystem::Mesh::FaceElemBase::clearLocalPhysicalCoordinateSystem()
{
    if (this->local_coordinate_system != NULL)
        delete this->local_coordinate_system;
    
    this->local_coordinate_system = NULL;
}


void
FESystem::Mesh::FaceElemBase::initializeLocalPhysicalCoordinateSystem()
{
    // make sure that the coordinate system has not already been set
    FESystemAssert0(this->local_coordinate_system == NULL, FESystem::Exception::InvalidState);
    // calculate the local coordinate system
    FESystem::Mesh::Node &n0 = this->getNode(0), &n1 = this->getNode(1), &n2 = this->getNode(2);
    
    FESystem::Numerics::DenseMatrix<FESystemDouble> mat;
    FESystem::Numerics::LocalVector<FESystemDouble> vec_x, vec_y, vec_z;
    
    mat.resize(n0.getSize(), n0.getSize());
    vec_x.resize(n0.getSize()); vec_y.resize(n0.getSize()); vec_z.resize(n0.getSize());
    vec_x.copyVector(n1);
    vec_x.add(-1.0, n0);
    vec_x.scaleToUnitLength();
    
    // this is only temporary initialization of the point
    vec_y.copyVector(n2);
    vec_y.add(-1.0, n0);
    
    vec_x.crossProduct(vec_y, vec_z); // z-axis direction
    vec_z.scaleToUnitLength();
    
    // now calculate the actual y-axis
    vec_z.crossProduct(vec_x, vec_y);
    vec_y.scaleToUnitLength();
    
    mat.setColumnVals(0, 0, n0.getSize()-1, vec_x);
    mat.setColumnVals(1, 0, n0.getSize()-1, vec_y);
    mat.setColumnVals(2, 0, n0.getSize()-1, vec_z);
    
    this->local_coordinate_system = new FESystem::Geometry::RectangularCoordinateSystem(n0, mat);
}


FESystem::Mesh::QuadElemBase::QuadElemBase(FESystemUInt nnodes, FESystem::Mesh::ElementType type, FESystemBoolean local_cs_same_as_global):
FESystem::Mesh::FaceElemBase(nnodes, type, local_cs_same_as_global)
{
    // initialize the boundary ID map
    if (FESystem::Mesh::QuadElemBase::quad_boundary_node_set.get() == NULL)
    {
        FESystem::Mesh::QuadElemBase::quad_boundary_node_set.reset(new std::map<FESystemUInt, std::vector<FESystemUInt> >());
        // now add the boundary nodes
        std::vector<FESystemUInt> b_vec;
        for (FESystemUInt i=0; i<4; i++)
        {
            b_vec.clear();
            switch (i) 
            {                    
                case 0:
                    b_vec.push_back(0);
                    b_vec.push_back(1);
                    break;

                case 1:
                    b_vec.push_back(1);
                    b_vec.push_back(2);
                    break;

                case 2:
                    b_vec.push_back(2);
                    b_vec.push_back(3);
                    break;

                case 3:
                    b_vec.push_back(3);
                    b_vec.push_back(0);
                    break;

                default:
                    FESystemAssert1(false, FESystem::Exception::InvalidID, i);
                    break;
            }
            FESystem::Mesh::QuadElemBase::quad_boundary_node_set->insert(std::map<FESystemUInt, std::vector<FESystemUInt> >::value_type(i, b_vec));
        }
    }

}


FESystem::Mesh::QuadElemBase::~QuadElemBase()
{
    
}





void 
FESystem::Mesh::QuadElemBase::getConstantCoordinateIDAndValueForBoundary(const FESystemUInt b_id, FESystemUInt& coord_id, FESystemDouble& coord_val) const
{
    switch (b_id) {
            
        case 0: // bottom boundary
            coord_id = 1;
            coord_val = -1.0; 
            break;
            
        case 1: // right boundary
            coord_id = 0;
            coord_val = 1.0; 
            break;
            
        case 2: // upper boundary
            coord_id = 1;
            coord_val = 1.0; 
            break;
            
        case 3: // left boundary
            coord_id = 0;
            coord_val = -1.0; 
            break;
            
        default:
            FESystemAssert1(false, FESystem::Exception::InvalidID, b_id);
            break;
    }
}


const std::map<FESystemUInt, std::vector<FESystemUInt> >&
FESystem::Mesh::QuadElemBase::getBoundaryIDAndBoundaryNodeMap() const
{
    return *FESystem::Mesh::QuadElemBase::quad_boundary_node_set;
}




FESystem::Mesh::TriElemBase::TriElemBase(FESystemUInt nnodes, FESystem::Mesh::ElementType type, FESystemBoolean local_cs_same_as_global):
FESystem::Mesh::FaceElemBase(nnodes, type, local_cs_same_as_global)
{
    // initialize the boundary ID map
    if (FESystem::Mesh::TriElemBase::tri_boundary_node_set.get() == NULL)
    {
        FESystem::Mesh::TriElemBase::tri_boundary_node_set.reset(new std::map<FESystemUInt, std::vector<FESystemUInt> >());
        // now add the boundary nodes
        std::vector<FESystemUInt> b_vec;
        for (FESystemUInt i=0; i<3; i++)
        {
            b_vec.clear();
            switch (i) 
            {                    
                case 0:
                    b_vec.push_back(0);
                    b_vec.push_back(1);
                    break;
                    
                case 1:
                    b_vec.push_back(1);
                    b_vec.push_back(2);
                    break;
                    
                case 2:
                    b_vec.push_back(2);
                    b_vec.push_back(0);
                    break;
                    
                default:
                    FESystemAssert1(false, FESystem::Exception::InvalidID, i);
                    break;
            }
            FESystem::Mesh::TriElemBase::tri_boundary_node_set->insert(std::map<FESystemUInt, std::vector<FESystemUInt> >::value_type(i, b_vec));
        }
    }
        
}



FESystem::Mesh::TriElemBase::~TriElemBase()
{
    if (this->parent_nondegenerate_elem != NULL)
        delete this->parent_nondegenerate_elem;
}


const FESystem::Utility::Table<FESystemUInt>& 
FESystem::Mesh::TriElemBase::getPointDimensionIDTable() const
{
    // this does not exist for the triangular element. 
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}


FESystemUInt 
FESystem::Mesh::TriElemBase::getNNodesAlongDimension(FESystemUInt n) const
{
    // this does not exist for the triangular element. 
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}



FESystemUInt 
FESystem::Mesh::TriElemBase::getLocalNodeIDALongDim(FESystemUInt n_id, FESystemUInt d) const
{
    // this does not exist for the triangular element. 
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}




void
FESystem::Mesh::TriElemBase::getLocalComputationalCoordinateForLocalNodeID(FESystemUInt id, FESystem::Geometry::Point& p) const
{
    // this does not exist for the triangular element. 
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}


bool
FESystem::Mesh::TriElemBase::ifDegerateElement() const
{
    return true;
}


const FESystem::Mesh::ElemBase& 
FESystem::Mesh::TriElemBase::getParentNondegenerateElem() const
{
    FESystemAssert0(this->parent_nondegenerate_elem != NULL, FESystem::Exception::InvalidState);
    
    return *(this->parent_nondegenerate_elem);
}


void 
FESystem::Mesh::TriElemBase::getConstantCoordinateIDAndValueForBoundary(const FESystemUInt b_id, FESystemUInt& coord_id, FESystemDouble& coord_val) const
{
    switch (b_id) {
            
        case 0: // bottom boundary
            coord_id = 1;
            coord_val = -1.0; 
            break;
            
        case 1: // right boundary
            coord_id = 0;
            coord_val = 1.0; 
            break;
            
        case 2: // left boundary
            coord_id = 0;
            coord_val = -1.0; 
            break;
            
        default:
            FESystemAssert1(false, FESystem::Exception::InvalidID, b_id);
            break;
    }
}


const std::map<FESystemUInt, std::vector<FESystemUInt> >& 
FESystem::Mesh::TriElemBase::getBoundaryIDAndBoundaryNodeMap() const
{
    return *FESystem::Mesh::TriElemBase::tri_boundary_node_set;
}


