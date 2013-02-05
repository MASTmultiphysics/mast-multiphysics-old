
/*
 *  VolumeElemBase.cpp
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/27/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

// FESystem includes
#include "Mesh/VolumeElemBase.h"
#include "Mesh/Node.h"
#include "Numerics/DenseMatrix.h"
#include "Numerics/LocalVector.h"
#include "Geom/RectangularCoordinateSystem.h"

// initialization of the pointer here
std::auto_ptr<std::map<FESystemUInt, std::vector<FESystemUInt> > > FESystem::Mesh::HexElemBase::hex_boundary_node_set(NULL);
std::auto_ptr<std::map<FESystemUInt, std::vector<FESystemUInt> > > FESystem::Mesh::PrismElemBase::prism_boundary_node_set(NULL);
std::auto_ptr<std::map<FESystemUInt, std::vector<FESystemUInt> > > FESystem::Mesh::TetElemBase::tet_boundary_node_set(NULL);
std::auto_ptr<std::map<FESystemUInt, std::vector<FESystemUInt> > > FESystem::Mesh::PyramidElemBase::pyramid_boundary_node_set(NULL);



FESystem::Mesh::VolumeElemBase::VolumeElemBase(FESystemUInt nnodes, FESystem::Mesh::ElementType type, FESystemBoolean local_cs_same_as_global):
FESystem::Mesh::ElemBase(nnodes, type, local_cs_same_as_global)
{
    
}


FESystem::Mesh::VolumeElemBase::~VolumeElemBase()
{
    
}

FESystemUInt
FESystem::Mesh::VolumeElemBase::getDimension() const
{
    return 3;
}



void
FESystem::Mesh::VolumeElemBase::clearLocalPhysicalCoordinateSystem()
{
    if (this->local_coordinate_system != NULL)
        delete this->local_coordinate_system;
    
    this->local_coordinate_system = NULL;
}


void
FESystem::Mesh::VolumeElemBase::initializeLocalPhysicalCoordinateSystem()
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



void
FESystem::Mesh::VolumeElemBase::calculateSurfaceNormal(FESystem::Numerics::VectorBase<FESystemDouble> &n_vec) const
{
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}




FESystem::Mesh::HexElemBase::HexElemBase(FESystemUInt nnodes, FESystem::Mesh::ElementType type, FESystemBoolean local_cs_same_as_global):
FESystem::Mesh::VolumeElemBase(nnodes, type, local_cs_same_as_global)
{
    // initialize the boundary ID map
    if (FESystem::Mesh::HexElemBase::hex_boundary_node_set.get() == NULL)
    {
        FESystem::Mesh::HexElemBase::hex_boundary_node_set.reset(new std::map<FESystemUInt, std::vector<FESystemUInt> >());
        // now add the boundary nodes
        std::vector<FESystemUInt> b_vec;
        for (FESystemUInt i=0; i<4; i++)
        {
            b_vec.clear();
            switch (i)
            {
                case 0: // back
                    b_vec.push_back(3);
                    b_vec.push_back(2);
                    b_vec.push_back(1);
                    b_vec.push_back(0);
                    break;
                    
                case 1: // left
                    b_vec.push_back(0);
                    b_vec.push_back(4);
                    b_vec.push_back(7);
                    b_vec.push_back(3);
                    break;
                    
                case 2: // front
                    b_vec.push_back(4);
                    b_vec.push_back(5);
                    b_vec.push_back(6);
                    b_vec.push_back(7);
                    break;
                    
                case 3: // right
                    b_vec.push_back(1);
                    b_vec.push_back(2);
                    b_vec.push_back(6);
                    b_vec.push_back(5);
                    break;
                    
                case 4: // top
                    b_vec.push_back(2);
                    b_vec.push_back(3);
                    b_vec.push_back(7);
                    b_vec.push_back(6);
                    break;

                case 5: // bottom
                    b_vec.push_back(0);
                    b_vec.push_back(1);
                    b_vec.push_back(5);
                    b_vec.push_back(4);
                    break;

                default:
                    FESystemAssert1(false, FESystem::Exception::InvalidID, i);
                    break;
            }
            FESystem::Mesh::HexElemBase::hex_boundary_node_set->insert(std::map<FESystemUInt, std::vector<FESystemUInt> >::value_type(i, b_vec));
        }
    }
    
}


FESystem::Mesh::HexElemBase::~HexElemBase()
{
    
}




FESystemUInt
FESystem::Mesh::HexElemBase::getNBoundaries() const
{
    return 6;
}




void
FESystem::Mesh::HexElemBase::getConstantCoordinateIDAndValueForBoundary(const FESystemUInt b_id, FESystemUInt& coord_id, FESystemDouble& coord_val) const
{
    switch (b_id) {
        case 0: // back
            coord_id = 2;
            coord_val = -1.0;
            break;
            
        case 1: // left
            coord_id = 0;
            coord_val = -1.0;
            break;
            
        case 2: // front
            coord_id = 2;
            coord_val = 1.0;
            break;
            
        case 3: // right
            coord_id = 0;
            coord_val = 1.0;
            break;

        case 4: // top
            coord_id = 1;
            coord_val =  1.0;
            break;

        case 5: // bottom
            coord_id = 1;
            coord_val = -1.0;
            break;

        default:
            FESystemAssert1(false, FESystem::Exception::InvalidID, b_id);
            break;
    }
}


const std::map<FESystemUInt, std::vector<FESystemUInt> >&
FESystem::Mesh::HexElemBase::getBoundaryIDAndBoundaryNodeMap() const
{
    return *FESystem::Mesh::HexElemBase::hex_boundary_node_set;
}




FESystem::Mesh::PrismElemBase::PrismElemBase(FESystemUInt nnodes, FESystem::Mesh::ElementType type, FESystemBoolean local_cs_same_as_global):
FESystem::Mesh::VolumeElemBase(nnodes, type, local_cs_same_as_global)
{
    // initialize the boundary ID map
    if (FESystem::Mesh::PrismElemBase::prism_boundary_node_set.get() == NULL)
    {
        FESystem::Mesh::PrismElemBase::prism_boundary_node_set.reset(new std::map<FESystemUInt, std::vector<FESystemUInt> >());
        // now add the boundary nodes
        std::vector<FESystemUInt> b_vec;
        for (FESystemUInt i=0; i<3; i++)
        {
            b_vec.clear();
            switch (i)
            {
                case 0: // back
                    b_vec.push_back(3);
                    b_vec.push_back(2);
                    b_vec.push_back(1);
                    b_vec.push_back(0);
                    break;
                    
                case 1: // left
                    b_vec.push_back(3);
                    b_vec.push_back(0);
                    b_vec.push_back(4);
                    b_vec.push_back(5);
                    break;
                    
                case 2: // right
                    b_vec.push_back(1);
                    b_vec.push_back(2);
                    b_vec.push_back(5);
                    b_vec.push_back(4);
                    break;
                    
                case 3: // top
                    b_vec.push_back(2);
                    b_vec.push_back(3);
                    b_vec.push_back(5);
                    break;

                case 4: // bottom
                    b_vec.push_back(0);
                    b_vec.push_back(1);
                    b_vec.push_back(4);
                    break;

                default:
                    FESystemAssert1(false, FESystem::Exception::InvalidID, i);
                    break;
            }
            FESystem::Mesh::PrismElemBase::prism_boundary_node_set->insert(std::map<FESystemUInt, std::vector<FESystemUInt> >::value_type(i, b_vec));
        }
    }
    
}



FESystem::Mesh::PrismElemBase::~PrismElemBase()
{
    if (this->parent_nondegenerate_elem != NULL)
        delete this->parent_nondegenerate_elem;
}



FESystemUInt
FESystem::Mesh::PrismElemBase::getNBoundaries() const
{
    return 5;
}



const FESystem::Utility::Table<FESystemUInt>&
FESystem::Mesh::PrismElemBase::getPointDimensionIDTable() const
{
    // this does not exist for the triangular element.
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}


FESystemUInt
FESystem::Mesh::PrismElemBase::getNNodesAlongDimension(FESystemUInt n) const
{
    // this does not exist for the triangular element.
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}



FESystemUInt
FESystem::Mesh::PrismElemBase::getLocalNodeIDALongDim(FESystemUInt n_id, FESystemUInt d) const
{
    // this does not exist for the triangular element.
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}




void
FESystem::Mesh::PrismElemBase::getLocalComputationalCoordinateForLocalNodeID(FESystemUInt id, FESystem::Geometry::Point& p) const
{
    // this does not exist for the triangular element.
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}


bool
FESystem::Mesh::PrismElemBase::ifDegerateElement() const
{
    return true;
}


const FESystem::Mesh::ElemBase&
FESystem::Mesh::PrismElemBase::getParentNondegenerateElem() const
{
    FESystemAssert0(this->parent_nondegenerate_elem != NULL, FESystem::Exception::InvalidState);
    
    return *(this->parent_nondegenerate_elem);
}


void
FESystem::Mesh::PrismElemBase::getConstantCoordinateIDAndValueForBoundary(const FESystemUInt b_id, FESystemUInt& coord_id, FESystemDouble& coord_val) const
{
    switch (b_id) {
            
        case 0: // back
            coord_id = 2;
            coord_val = -1.0;
            break;
            
        case 1: // left
            coord_id = 0;
            coord_val = -1.0;
            break;
            
        case 2: // right
            coord_id = 0;
            coord_val =  1.0;
            break;

        case 3: // top
            coord_id = 1;
            coord_val =  1.0;
            break;

        case 4: // bottom
            coord_id = 1;
            coord_val = -1.0;
            break;

        default:
            FESystemAssert1(false, FESystem::Exception::InvalidID, b_id);
            break;
    }
}


const std::map<FESystemUInt, std::vector<FESystemUInt> >&
FESystem::Mesh::PrismElemBase::getBoundaryIDAndBoundaryNodeMap() const
{
    return *FESystem::Mesh::PrismElemBase::prism_boundary_node_set;
}




FESystem::Mesh::TetElemBase::TetElemBase(FESystemUInt nnodes, FESystem::Mesh::ElementType type, FESystemBoolean local_cs_same_as_global):
FESystem::Mesh::VolumeElemBase(nnodes, type, local_cs_same_as_global)
{
    // initialize the boundary ID map
    if (FESystem::Mesh::TetElemBase::tet_boundary_node_set.get() == NULL)
    {
        FESystem::Mesh::TetElemBase::tet_boundary_node_set.reset(new std::map<FESystemUInt, std::vector<FESystemUInt> >());
        // now add the boundary nodes
        std::vector<FESystemUInt> b_vec;
        for (FESystemUInt i=0; i<3; i++)
        {
            b_vec.clear();
            switch (i)
            {
                case 0: // back
                    b_vec.push_back(2);
                    b_vec.push_back(1);
                    b_vec.push_back(0);
                    break;
                    
                case 1: // left
                    b_vec.push_back(2);
                    b_vec.push_back(0);
                    b_vec.push_back(3);
                    break;
                    
                case 2: // right
                    b_vec.push_back(1);
                    b_vec.push_back(2);
                    b_vec.push_back(3);
                    break;

                case 3: // bottom
                    b_vec.push_back(0);
                    b_vec.push_back(1);
                    b_vec.push_back(3);
                    break;
                    
                default:
                    FESystemAssert1(false, FESystem::Exception::InvalidID, i);
                    break;
            }
            FESystem::Mesh::TetElemBase::tet_boundary_node_set->insert(std::map<FESystemUInt, std::vector<FESystemUInt> >::value_type(i, b_vec));
        }
    }
    
}



FESystem::Mesh::TetElemBase::~TetElemBase()
{
    if (this->parent_nondegenerate_elem != NULL)
        delete this->parent_nondegenerate_elem;
}



FESystemUInt
FESystem::Mesh::TetElemBase::getNBoundaries() const
{
    return 4;
}


const FESystem::Utility::Table<FESystemUInt>&
FESystem::Mesh::TetElemBase::getPointDimensionIDTable() const
{
    // this does not exist for the triangular element.
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}


FESystemUInt
FESystem::Mesh::TetElemBase::getNNodesAlongDimension(FESystemUInt n) const
{
    // this does not exist for the triangular element.
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}



FESystemUInt
FESystem::Mesh::TetElemBase::getLocalNodeIDALongDim(FESystemUInt n_id, FESystemUInt d) const
{
    // this does not exist for the triangular element.
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}




void
FESystem::Mesh::TetElemBase::getLocalComputationalCoordinateForLocalNodeID(FESystemUInt id, FESystem::Geometry::Point& p) const
{
    // this does not exist for the triangular element.
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}


bool
FESystem::Mesh::TetElemBase::ifDegerateElement() const
{
    return true;
}


const FESystem::Mesh::ElemBase&
FESystem::Mesh::TetElemBase::getParentNondegenerateElem() const
{
    FESystemAssert0(this->parent_nondegenerate_elem != NULL, FESystem::Exception::InvalidState);
    
    return *(this->parent_nondegenerate_elem);
}


void
FESystem::Mesh::TetElemBase::getConstantCoordinateIDAndValueForBoundary(const FESystemUInt b_id, FESystemUInt& coord_id, FESystemDouble& coord_val) const
{
    switch (b_id) {
            
        case 0: // back
            coord_id = 2;
            coord_val = -1.0;
            break;
            
        case 1: // left
            coord_id = 0;
            coord_val = -1.0;
            break;
            
        case 2: // right
            coord_id = 0;
            coord_val =  1.0;
            break;

        case 3: // bottom
            coord_id = 1;
            coord_val = -1.0;
            break;

        default:
            FESystemAssert1(false, FESystem::Exception::InvalidID, b_id);
            break;
    }
}


const std::map<FESystemUInt, std::vector<FESystemUInt> >&
FESystem::Mesh::TetElemBase::getBoundaryIDAndBoundaryNodeMap() const
{
    return *FESystem::Mesh::TetElemBase::tet_boundary_node_set;
}


FESystem::Mesh::PyramidElemBase::PyramidElemBase(FESystemUInt nnodes, FESystem::Mesh::ElementType type, FESystemBoolean local_cs_same_as_global):
FESystem::Mesh::VolumeElemBase(nnodes, type, local_cs_same_as_global)
{
    // initialize the boundary ID map
    if (FESystem::Mesh::PyramidElemBase::pyramid_boundary_node_set.get() == NULL)
    {
        FESystem::Mesh::PyramidElemBase::pyramid_boundary_node_set.reset(new std::map<FESystemUInt, std::vector<FESystemUInt> >());
        // now add the boundary nodes
        std::vector<FESystemUInt> b_vec;
        for (FESystemUInt i=0; i<3; i++)
        {
            b_vec.clear();
            switch (i)
            {
                case 0: // back
                    b_vec.push_back(3);
                    b_vec.push_back(2);
                    b_vec.push_back(1);
                    b_vec.push_back(0);
                    break;
                    
                case 1: // left
                    b_vec.push_back(3);
                    b_vec.push_back(0);
                    b_vec.push_back(4);
                    break;
                    
                case 2: // right
                    b_vec.push_back(1);
                    b_vec.push_back(2);
                    b_vec.push_back(4);
                    break;

                case 3: // top
                    b_vec.push_back(2);
                    b_vec.push_back(3);
                    b_vec.push_back(4);
                    break;

                case 4: // bottom
                    b_vec.push_back(0);
                    b_vec.push_back(1);
                    b_vec.push_back(4);
                    break;

                default:
                    FESystemAssert1(false, FESystem::Exception::InvalidID, i);
                    break;
            }
            FESystem::Mesh::PyramidElemBase::pyramid_boundary_node_set->insert(std::map<FESystemUInt, std::vector<FESystemUInt> >::value_type(i, b_vec));
        }
    }
    
}



FESystem::Mesh::PyramidElemBase::~PyramidElemBase()
{
    if (this->parent_nondegenerate_elem != NULL)
        delete this->parent_nondegenerate_elem;
}



FESystemUInt
FESystem::Mesh::PyramidElemBase::getNBoundaries() const
{
    return 5;
}


const FESystem::Utility::Table<FESystemUInt>&
FESystem::Mesh::PyramidElemBase::getPointDimensionIDTable() const
{
    // this does not exist for the triangular element.
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}


FESystemUInt
FESystem::Mesh::PyramidElemBase::getNNodesAlongDimension(FESystemUInt n) const
{
    // this does not exist for the triangular element.
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}



FESystemUInt
FESystem::Mesh::PyramidElemBase::getLocalNodeIDALongDim(FESystemUInt n_id, FESystemUInt d) const
{
    // this does not exist for the triangular element.
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}




void
FESystem::Mesh::PyramidElemBase::getLocalComputationalCoordinateForLocalNodeID(FESystemUInt id, FESystem::Geometry::Point& p) const
{
    // this does not exist for the triangular element.
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}


bool
FESystem::Mesh::PyramidElemBase::ifDegerateElement() const
{
    return true;
}


const FESystem::Mesh::ElemBase&
FESystem::Mesh::PyramidElemBase::getParentNondegenerateElem() const
{
    FESystemAssert0(this->parent_nondegenerate_elem != NULL, FESystem::Exception::InvalidState);
    
    return *(this->parent_nondegenerate_elem);
}


void
FESystem::Mesh::PyramidElemBase::getConstantCoordinateIDAndValueForBoundary(const FESystemUInt b_id, FESystemUInt& coord_id, FESystemDouble& coord_val) const
{
    switch (b_id) {
            
        case 0: // back
            coord_id = 2;
            coord_val = -1.0;
            break;
            
        case 1: // left
            coord_id = 0;
            coord_val = -1.0;
            break;
            
        case 2: // right
            coord_id = 0;
            coord_val =  1.0;
            break;

        case 3: // top
            coord_id = 1;
            coord_val =  1.0;
            break;

        case 4: // bottom
            coord_id = 1;
            coord_val = -1.0;
            break;

        default:
            FESystemAssert1(false, FESystem::Exception::InvalidID, b_id);
            break;
    }
}


const std::map<FESystemUInt, std::vector<FESystemUInt> >&
FESystem::Mesh::PyramidElemBase::getBoundaryIDAndBoundaryNodeMap() const
{
    return *FESystem::Mesh::PyramidElemBase::pyramid_boundary_node_set;
}


