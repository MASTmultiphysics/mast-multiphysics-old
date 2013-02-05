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
#include "Mesh/ElemBase.h"
#include "FiniteElems/FiniteElementBase.h"


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


void
FESystem::Quadrature::QuadratureBase::clear()
{
    this->if_initialized = false;
    this->dimension = 0;
    this->order = 0;
    
    if (this->origin != NULL)
        delete this->origin;
    this->origin = NULL;
    
    for (FESystemInt i=0; i<this->quadrature_points.size(); i++)
        if (this->quadrature_points[i] != NULL)
        {
            delete this->quadrature_points[i];
            this->quadrature_points[i] = NULL;
        }
    
    this->quadrature_points.clear();
    this->quadrature_point_weights.clear();
}


FESystemUInt
FESystem::Quadrature::QuadratureBase::getDimention() const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    return this->dimension;
}


FESystemUInt
FESystem::Quadrature::QuadratureBase::getOrder() const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    return this->order;
}
            
            

const std::vector<FESystem::Geometry::Point*>& 
FESystem::Quadrature::QuadratureBase::getQuadraturePoints() const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    return this->quadrature_points;
}


const std::vector<FESystemDouble>& 
FESystem::Quadrature::QuadratureBase::getQuadraturePointWeights() const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    return this->quadrature_point_weights;
}

void
FESystem::Quadrature::QuadratureBase::initializeOrigin(FESystemUInt dim)
{
    FESystemAssert0(!this->if_initialized, FESystem::Exception::InvalidState);

    // if the origin is already initialized, delete it and create a object
    if (this->origin != NULL)
        delete this->origin;
    
    this->origin = new FESystem::Geometry::Point(dim);
}



void
FESystem::Quadrature::initializeQRuleForElem(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, FESystem::Quadrature::QuadratureBase& qrule)
{
    qrule.clear();

    switch (fe.getFiniteElementType())
    {
        case FESystem::FiniteElement::FE_LAGRANGE:
            // isoparametric elements have the order of geometry and interpolation to be the same
            qrule.init(elem.getDimension(), elem.getGeometryOrder());
            break;
            
        default:
            FESystemAssert1(false, FESystem::Exception::EnumerationNotHandled, fe.getFiniteElementType());
            break;
    }
}



void
FESystem::Quadrature::initializeQRuleForElemBoundary(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, FESystem::Quadrature::QuadratureBase& qrule)
{
    qrule.clear();
    
    switch (fe.getFiniteElementType())
    {
        case FESystem::FiniteElement::FE_LAGRANGE:
            // isoparametric elements have the order of geometry and interpolation to be the same
            qrule.init(elem.getDimension()-1, elem.getGeometryOrder());
            break;
            
        default:
            FESystemAssert1(false, FESystem::Exception::EnumerationNotHandled, fe.getFiniteElementType());
            break;
    }
}

