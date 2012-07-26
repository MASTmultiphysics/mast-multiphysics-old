//
//  TrapezoidQuadrature.cpp
//  FESystem
//
//  Created by Manav Bhatia on 4/9/12.
//  Copyright (c) 2012. All rights reserved.
//

// FESystem includes
#include "Quadrature/TrapezoidQuadrature.h"
#include "Base/FESystemExceptions.h"
#include "Base/FESystemBase.h"
#include "Geom/Point.h"


FESystem::Quadrature::TrapezoidQuadrature::TrapezoidQuadrature():
FESystem::Quadrature::QuadratureBase()
{
    
}


FESystem::Quadrature::TrapezoidQuadrature::~TrapezoidQuadrature()
{
    
}


void 
FESystem::Quadrature::TrapezoidQuadrature::init(FESystemUInt dim, FESystemUInt order)
{
    FESystemAssert0(!this->if_initialized, FESystem::Exception::InvalidState);
    
    this->initializeOrigin(dim);
    
    // the total area is 4, which is uniformly divided into integration points
    // order of integration is transformed into number of points using an a rule that N points 
    // provide an accurate integration of a polynomial of order N-1.
    
    this->order = order;
    this->dimension = dim;
    
    FESystemUInt n_pts_per_dim = (order+1),
    n_pts_total = pow(n_pts_per_dim,dim);
    
    this->quadrature_points.resize(n_pts_total);
    this->quadrature_point_weights.resize(n_pts_total);
    
    FESystem::Geometry::Point* pt=NULL;
    FESystemDouble dx = 2.0/n_pts_per_dim, weight = pow(dx,dim);
    FESystemUInt i_pt=0;
    
    // iterate on the number of dimensions and add the points
    switch (dim) {
        case 1:
            for (FESystemUInt i=0; i<n_pts_per_dim; i++)
            {
                pt = new FESystem::Geometry::Point(this->origin->getCoordinateSystem());
                pt->setVal(0,(i+0.5)*dx-1.0);
                this->quadrature_points[i_pt] = pt;
                this->quadrature_point_weights[i_pt] = weight;
                i_pt++;
                pt = NULL;
            }
            break;
            
        case 2:
            for (FESystemUInt i=0; i<n_pts_per_dim; i++)
                for (FESystemUInt j=0; j<n_pts_per_dim; j++)
                {
                    pt = new FESystem::Geometry::Point(this->origin->getCoordinateSystem());
                    pt->setVal(0,(i+0.5)*dx-1.0);
                    pt->setVal(1,(j+0.5)*dx-1.0);
                    this->quadrature_points[i_pt] = pt;
                    this->quadrature_point_weights[i_pt] = weight;
                    i_pt++;
                    pt = NULL;
                }
            break;

        case 3:
            for (FESystemUInt i=0; i<n_pts_per_dim; i++)
                for (FESystemUInt j=0; j<n_pts_per_dim; j++)
                    for (FESystemUInt k=0; k<n_pts_per_dim; k++)
                    {
                        pt = new FESystem::Geometry::Point(this->origin->getCoordinateSystem());
                        pt->setVal(0,(i+0.5)*dx-1.0);
                        pt->setVal(1,(j+0.5)*dx-1.0);
                        pt->setVal(2,(k+0.5)*dx-1.0);
                        this->quadrature_points[i_pt] = pt;
                        this->quadrature_point_weights[i_pt] = weight;
                        i_pt++;
                        pt = NULL;
                    }
            break;

        default:
            FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
            break;
    }
    
}
