//
//  QuadratureBase.cpp
//  FESystem
//
//  Created by Manav Bhatia on 4/9/12.
//  Copyright (c) 2012. All rights reserved.
//

// FESystem includes
#include "Quadrature/QuadratureBase.h"
#include "Geom/Point.h"

FESystem::Quadrature::QuadratureBase::QuadratureBase():
if_initialized(false),
dimension(0),
order(0),
origin(NULL)
{
    
}

FESystem::Quadrature::QuadratureBase::~QuadratureBase()
{
    // iterate over the points in the vector and delete them
    for (FESystemInt i=0; i<this->quadrature_points.size(); i++)
        if (this->quadrature_points[i] != NULL)
            delete this->quadrature_points[i];
    
    if (this->origin != NULL)
        delete this->origin;    
}
            

FESystemUInt 
FESystem::Quadrature::QuadratureBase::getDimention() const
{
    return this->dimension;
}
            

FESystemUInt 
FESystem::Quadrature::QuadratureBase::getOrder() const
{
    return this->order;
}
            
            

const std::vector<FESystem::Geometry::Point*>& 
FESystem::Quadrature::QuadratureBase::getQuadraturePoints() const
{
    return this->quadrature_points;
}


const std::vector<FESystemDouble>& 
FESystem::Quadrature::QuadratureBase::getQuadraturePointWeights() const
{
    return this->quadrature_point_weights;
}

void
FESystem::Quadrature::QuadratureBase::initializeOrigin(FESystemUInt dim)
{
    // if the origin is already initialized, delete it and create a object
    if (this->origin != NULL)
        delete this->origin;
    
    this->origin = new FESystem::Geometry::Point(dim);
}

