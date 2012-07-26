
/*
 *  EdgeElemBase.cpp
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/27/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

// FESystem includes
#include "Mesh/EdgeElemBase.h"
#include "Mesh/Node.h"
#include "Numerics/LocalVector.h"
#include "Numerics/DenseMatrix.h"
#include "Geom/RectangularCoordinateSystem.h"


// initialization of the pointer here
std::auto_ptr<std::map<FESystemUInt, std::vector<FESystemUInt> > > FESystem::Mesh::EdgeElemBase::edge_boundary_node_set(NULL);



FESystem::Mesh::EdgeElemBase::EdgeElemBase(FESystemUInt n_nodes, FESystem::Mesh::ElementType type):
FESystem::Mesh::ElemBase(n_nodes, type)
{
    this->y_unit_vec = new FESystem::Numerics::LocalVector<FESystemDouble>;
    this->y_unit_vec->resize(3);
}




FESystem::Mesh::EdgeElemBase::~EdgeElemBase()
{
    delete this->y_unit_vec;
}


FESystemUInt 
FESystem::Mesh::EdgeElemBase::getDimension() const
{
    return 1;
}



void 
FESystem::Mesh::EdgeElemBase::setVectorForXYPlane(const FESystem::Numerics::VectorBase<FESystemDouble>& vec)
{
    // make sure that the dimensions are appropriate
    FESystemAssert2(vec.getSize() == this->y_unit_vec->getSize(), FESystem::Exception::DimensionsDoNotMatch, vec.getSize(), this->y_unit_vec->getSize());
    
    // make sure that vec has a finite length
    FESystemAssert0(vec.getL2Norm() > 0, FESystem::Exception::InvalidValue);
    
    this->y_unit_vec->copyVector(vec);
}


bool 
FESystem::Mesh::EdgeElemBase::ifDegerateElement() const
{
    return false;
}




const FESystem::Mesh::ElemBase& 
FESystem::Mesh::EdgeElemBase::getParentNondegenerateElem() const
{
    return *this;
}



void 
FESystem::Mesh::EdgeElemBase::getConstantCoordinateIDAndValueForBoundary(const FESystemUInt b_id, FESystemUInt& coord_id, FESystemDouble& coord_val) const
{
    switch (b_id) {
            
        case 0: // left point
            coord_id = 0;
            coord_val = -1.0; 
            break;
            
        case 1: // right point
            coord_id = 0;
            coord_val = 1.0; 
            break;
            
        default:
            FESystemAssert1(false, FESystem::Exception::InvalidID, b_id);
            break;
    }
}


void 
FESystem::Mesh::EdgeElemBase::initializeParentNondegenerateElement()
{
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}


const std::map<FESystemUInt, std::vector<FESystemUInt> >& 
FESystem::Mesh::EdgeElemBase::getBoundaryIDAndBoundaryNodeMap() const
{
    return *FESystem::Mesh::EdgeElemBase::edge_boundary_node_set;
}



void
FESystem::Mesh::EdgeElemBase::initializeLocalPhysicalCoordinateSystem()
{
    // make sure that the coordinate system has not already been set
    FESystemAssert0(this->local_coordinate_system == NULL, FESystem::Exception::InvalidState);
    FESystemAssert0(this->y_unit_vec->getL2Norm() > 0, FESystem::Exception::InvalidValue);
    
    // calculate the local coordinate system
    FESystem::Mesh::Node &n0 = this->getNode(0), &n1 = this->getNode(1);
    FESystem::Numerics::DenseMatrix<FESystemDouble> mat;
    FESystem::Numerics::LocalVector<FESystemDouble> vec_x, vec_y, vec_z;    
    mat.resize(n0.getSize(), n0.getSize());
    vec_x.resize(n0.getSize()); vec_y.resize(n0.getSize()); vec_z.resize(n0.getSize());
    vec_x.copyVector(n1);
    vec_x.add(-1.0, n0);
    vec_x.scaleToUnitLength();
    
    vec_x.crossProduct(*(this->y_unit_vec), vec_z); // z-axis direction
    vec_z.scaleToUnitLength();
    
    // now calculate the actual y-axis
    vec_z.crossProduct(vec_x, vec_y);
    vec_y.scaleToUnitLength();
    
    mat.setColumnVals(0, 0, n0.getSize()-1, vec_x);
    mat.setColumnVals(1, 0, n0.getSize()-1, vec_y);
    mat.setColumnVals(2, 0, n0.getSize()-1, vec_z);
    
    this->local_coordinate_system = new FESystem::Geometry::RectangularCoordinateSystem(n0, mat);
}

