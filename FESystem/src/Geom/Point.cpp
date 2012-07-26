
/*
 *  Point.cpp
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/27/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

// FESystem includes
#include "Base/FESystemBase.h"
#include "Geom/Point.h"
#include "Geom/CoordinateSystemBase.h"
#include "Geom/RectangularCoordinateSystem.h"
#include "Numerics/DenseMatrix.h"
            

FESystem::Geometry::Point::Point(FESystemUInt dim):
FESystem::Numerics::LocalVector<FESystemDouble>(dim),
if_global_origin(true),
coordinate_system(NULL)
{
    if (if_global_origin)
    {
        // the basis is an identity matrix
        FESystem::Numerics::DenseMatrix<FESystemDouble> mat;
        mat.resize(dim, dim); mat.setToIdentity();        
        this->coordinate_system = new FESystem::Geometry::RectangularCoordinateSystem(*this, mat);
        this->zero(); // the value of the global origin is always zero
    }
}


FESystem::Geometry::Point::Point(const FESystem::Geometry::CoordinateSystemBase& cs):
FESystem::Numerics::LocalVector<FESystemDouble>(),
if_global_origin(false),
coordinate_system(&cs)
{
    this->resize(cs.getDimension());
    this->zero();
}


FESystem::Geometry::Point::~Point()
{
    if (if_global_origin) 
        delete this->coordinate_system;
}
            


const FESystem::Geometry::CoordinateSystemBase& 
FESystem::Geometry::Point::getCoordinateSystem() const
{
    return *(this->coordinate_system);
}
            

