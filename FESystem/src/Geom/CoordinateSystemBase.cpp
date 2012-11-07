
/*
 *  CoordinateSystemBase.cpp
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/27/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

// FESystem includes
#include "Geom/CoordinateSystemBase.h"
#include "Numerics/VectorBase.h"
#include "Numerics/MatrixBase.h"
#include "Geom/Point.h"
#include "Functions/FunctionMappingBase.h"


FESystem::Geometry::CoordinateSystemBase::CoordinateSystemBase(const FESystem::Geometry::Point& o):
origin(o),
coordinate_mapping(NULL)
{

}


FESystem::Geometry::CoordinateSystemBase::~CoordinateSystemBase()
{
    // nothing specific needs to be done. The coordinate basis is an auto_ptr and will destroy itself. 
}
    


const FESystem::Geometry::Point& 
FESystem::Geometry::CoordinateSystemBase::getOrigin() const
{
    return this->origin;
}


FESystemUInt 
FESystem::Geometry::CoordinateSystemBase::getDimension() const
{
    return this->origin.getSize();
}

void
FESystem::Geometry::CoordinateSystemBase::mapPointInSelf(const FESystem::Geometry::Point& p, FESystem::Numerics::VectorBase<FESystemDouble>& vec) const
{
    // TODO: revisit for ensuring hierarchy of coordinate systems

    FESystem::Numerics::LocalVector<FESystemDouble> vec2;
    vec2.resize(this->getDimension());
    vec2.copyVector(p);
    vec2.add(-1.0, this->getOrigin());
    
    this->coordinate_mapping->inverseMap(vec2, vec);
}


const FESystem::Functions::FunctionMappingBase<FESystemDouble>& 
FESystem::Geometry::CoordinateSystemBase::getFunctionMappingObject() const
{
    FESystemAssert0(this->coordinate_mapping != NULL, FESystem::Exception::InvalidState);
    return *(this->coordinate_mapping);
}
