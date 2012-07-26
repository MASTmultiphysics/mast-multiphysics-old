
/*
 *  RectangularCoordinateSystem.cpp
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/27/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

// FESystem includes
#include "Geom/RectangularCoordinateSystem.h"
#include "Functions/LinearFunctionMapping.h"


FESystem::Geometry::RectangularCoordinateSystem::RectangularCoordinateSystem(const FESystem::Geometry::Point& origin, const FESystem::Numerics::MatrixBase<FESystemDouble>& basis):
FESystem::Geometry::CoordinateSystemBase(origin)
{
    FESystem::Functions::LinearFunctionMapping<FESystemDouble>* map = new FESystem::Functions::LinearFunctionMapping<FESystemDouble>;
    this->coordinate_mapping = map;
    map->reinit(basis);
}
    
    
FESystem::Geometry::RectangularCoordinateSystem::~RectangularCoordinateSystem()
{
    delete this->coordinate_mapping;
}
    


void
FESystem::Geometry::RectangularCoordinateSystem::mapCoordinates(const FESystem::Geometry::Point& p_in, 
                                                                     FESystem::Geometry::Point& p_mapped) const
{
    
}


void 
FESystem::Geometry::RectangularCoordinateSystem::inverseMapCoordinates(const FESystem::Geometry::Point& p_mapped, 
                                                                            FESystem::Geometry::Point& p_inv_mapped) const
{
    
}
