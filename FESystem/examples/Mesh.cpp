//
//  MeshTest.cpp
//  FESystemApplication
//
//  Created by Manav Bhatia on 3/19/12.
//  Copyright (c) 2012. All rights reserved.
//

// FESystem include
#include "TestingIncludes.h"
#include <cmath>



void distributePoints(const FESystemUInt n_divs, const std::vector<FESystemDouble>& div_locations, const std::vector<FESystemUInt>& n_subdivs_in_div, const std::vector<FESystemDouble>& relative_mesh_size_at_div, std::vector<FESystemDouble>& points)
{
    FESystemAssert2(div_locations.size() == n_divs+1, FESystem::Exception::DimensionsDoNotMatch, div_locations.size(), n_divs+1);
    FESystemAssert2(relative_mesh_size_at_div.size() == n_divs+1, FESystem::Exception::DimensionsDoNotMatch, relative_mesh_size_at_div.size(), n_divs+1);
    FESystemAssert2(n_subdivs_in_div.size() == n_divs, FESystem::Exception::DimensionsDoNotMatch, n_subdivs_in_div.size(), n_divs);
    
    // calculate total number of points
    FESystemUInt n_total_points = 1;
    for (FESystemUInt i=0; i<n_divs; i++)
        n_total_points += n_subdivs_in_div[i];
    
    // resize the points vector and set the first and last points of each division
    points.resize(n_total_points);

    FESystemUInt n=1;
    points[0] = div_locations[0];
    for (FESystemUInt i=0; i<n_divs; i++)
    {
        n += n_subdivs_in_div[i];
        points[n-1] = div_locations[i+1];
    }

    n=1;
    FESystemDouble dx=0.0, growth_factor = 0.0;
    // now calculate the base mesh size, and calculate the nodal points
    for (FESystemUInt i=0; i<n_divs; i++)
    {
        growth_factor = pow(relative_mesh_size_at_div[i+1]/relative_mesh_size_at_div[i], 1.0/(n_subdivs_in_div[i]-1.0));
        if (fabs(growth_factor-1.0)>1.0e-10)
            dx = (div_locations[i+1]-div_locations[i]) * (1.0-growth_factor)/(1.0-pow(growth_factor, n_subdivs_in_div[i]));
        else
        {
            growth_factor = 1.0;
            dx = (div_locations[i+1]-div_locations[i]) / n_subdivs_in_div[i];
        }
        
        for (FESystemUInt n_pt=1; n_pt<n_subdivs_in_div[i]; n_pt++)
        {
            points[n+n_pt-1] = points[n+n_pt-2] + dx;
            dx *= growth_factor;
        }
        n += n_subdivs_in_div[i];
    }
}




void createLineMesh(FESystem::Mesh::ElementType elem_type, FESystem::Mesh::MeshBase& mesh, FESystem::Geometry::Point& origin,
                    const std::vector<FESystemDouble>& points, FESystemUInt& n_elem_nodes, MeshType m_type, FESystemBoolean local_cs_same_as_global)
{
    // create a nx x ny grid of nodes, and connect them by quad4 elements
    FESystemUInt n_nodes, n_elems;
    
    const FESystem::Geometry::CoordinateSystemBase& global_cs = origin.getCoordinateSystem();
    FESystem::Numerics::LocalVector<FESystemDouble> vec;
    vec.resize(3); vec.zero(); vec.setVal(1, 1.0);
    
    std::auto_ptr<std::vector<FESystem::Mesh::Node*> > nodes;
    std::auto_ptr<std::vector<FESystem::Mesh::ElemBase*> > elems;
    
    // set the location of individual nodes
    switch (elem_type)
    {
        case FESystem::Mesh::EDGE2:
        {
            n_nodes=points.size();
            n_elems=(n_nodes-1);
            n_elem_nodes = 2;
            
            nodes.reset(mesh.createNodes(n_nodes, global_cs).release());
            elems.reset(mesh.createElements(n_elems, elem_type, local_cs_same_as_global).release());
            
            FESystemUInt id=0;
            FESystem::Mesh::Node* node_p;
            for (FESystemUInt ix=0; ix<n_nodes; ix++)
            {
                node_p = (*nodes)[id];
                node_p->setVal(0,points[ix]);
                node_p->setExternalID(id);
                id++;
            }
            
            
            // set the nodes for each element
            id=0;
            FESystem::Mesh::ElemBase* elem_p;
            for (FESystemUInt ix=0; ix<(n_nodes-1); ix++)
            {
                elem_p = (*elems)[id];
                dynamic_cast<FESystem::Mesh::EdgeElemBase*>(elem_p)->setVectorForXYPlane(vec);
                elem_p->setNode(0, *(*nodes)[ix]);
                elem_p->setNode(1, *(*nodes)[ix+1]);
                elem_p->setExternalID(id);
                id++;
            }
        }
            break;
            
        case FESystem::Mesh::EDGE3:
        {
            n_nodes=points.size();
            n_elems=n_nodes/2-1;
            n_elem_nodes = 3;

            nodes.reset(mesh.createNodes(n_nodes, global_cs).release());
            elems.reset(mesh.createElements(n_elems, elem_type, local_cs_same_as_global).release());
            
            FESystemUInt id=0;
            FESystem::Mesh::Node* node_p;
            for (FESystemUInt ix=0; ix<n_nodes; ix++)
            {
                node_p = (*nodes)[id];
                node_p->setVal(0,points[ix]);
                node_p->setExternalID(id);
                id++;
            }
            
            
            // set the location of individual nodes
            id=0;
            FESystem::Mesh::ElemBase* elem_p;
            for (FESystemUInt ix=0; ix<n_elems; ix++)
            {
                elem_p = (*elems)[id];
                dynamic_cast<FESystem::Mesh::EdgeElemBase*>(elem_p)->setVectorForXYPlane(vec);
                elem_p->setNode(0, *(*nodes)[2*ix]);
                elem_p->setNode(1, *(*nodes)[2*(ix+1)]);
                elem_p->setNode(2, *(*nodes)[2*ix+1]);
                elem_p->setExternalID(id);
                id++;
            }
        }
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
    }
    mesh.reinit();
    
}




void createPlaneMesh(FESystem::Mesh::ElementType elem_type, FESystem::Mesh::MeshBase& mesh, FESystem::Geometry::Point& origin,
                     const std::vector<FESystemDouble>& x_points, const std::vector<FESystemDouble>& y_points, FESystemUInt& n_elem_nodes,
                     MeshType m_type, FESystemBoolean local_cs_same_as_global)
{
    // create a nx x ny grid of nodes, and connect them by quad4 elements
    FESystemUInt nx, ny, n_nodes, n_elems;
    
    const FESystem::Geometry::CoordinateSystemBase& global_cs = origin.getCoordinateSystem();
    
    std::auto_ptr<std::vector<FESystem::Mesh::Node*> > nodes;
    std::auto_ptr<std::vector<FESystem::Mesh::ElemBase*> > elems;
    
    // create the nodes
    nx = x_points.size();
    ny = y_points.size();
    nodes.reset(mesh.createNodes(nx*ny, global_cs).release());
    FESystemUInt id=0;
    FESystem::Mesh::Node* node_p;
    for (FESystemUInt iy=0; iy<ny; iy++)
        for (FESystemUInt ix=0; ix<nx; ix++)
        {
            node_p = (*nodes)[id];
            node_p->setVal(0,x_points[ix]);
            node_p->setVal(1,y_points[iy]);
            node_p->setExternalID(id);
            id++;
        }
    
    
    // set the location of individual nodes
    switch (elem_type)
    {
        case FESystem::Mesh::QUAD4:
        {
            n_elems=(nx-1)*(ny-1);
            n_elem_nodes = 4;
            
            elems.reset(mesh.createElements(n_elems, elem_type, local_cs_same_as_global).release());
            
            // set the location of individual nodes
            id=0;
            FESystem::Mesh::ElemBase* elem_p;
            for (FESystemUInt iy=0; iy<(ny-1); iy++)
                for (FESystemUInt ix=0; ix<(nx-1); ix++)
                {
                    elem_p = (*elems)[id];
                    elem_p->setNode(0, *(*nodes)[(iy*nx+ix)]);
                    elem_p->setNode(1, *(*nodes)[(iy*nx+ix+1)]);
                    elem_p->setNode(2, *(*nodes)[((iy+1)*nx+ix+1)]);
                    elem_p->setNode(3, *(*nodes)[((iy+1)*nx+ix)]);
                    elem_p->setExternalID(id);
                    id++;
                }
        }
            break;
            
        case FESystem::Mesh::QUAD9:
        {
            FESystemAssert0(fmod(1.0*(nx-1),2.0) == 0.0, FESystem::Exception::InvalidValue);
            FESystemAssert0(fmod(1.0*(ny-1),2.0) == 0.0, FESystem::Exception::InvalidValue);
            
            FESystemUInt nx_elems = (nx-1)/2, ny_elems = (ny-1)/2;
            n_elems=nx_elems*ny_elems;
            n_elem_nodes = 9;
            
            elems.reset(mesh.createElements(n_elems, elem_type, local_cs_same_as_global).release());
                        
            // set the location of individual nodes
            id=0;
            FESystem::Mesh::ElemBase* elem_p;
            for (FESystemUInt iy=0; iy<nx_elems; iy++)
                for (FESystemUInt ix=0; ix<ny_elems; ix++)
                {
                    elem_p = (*elems)[id];
                    elem_p->setNode(0, *(*nodes)[2*iy*nx+2*ix]);
                    elem_p->setNode(1, *(*nodes)[2*iy*nx+2*(ix+1)]);
                    elem_p->setNode(2, *(*nodes)[2*(iy+1)*nx+2*(ix+1)]);
                    elem_p->setNode(3, *(*nodes)[2*(iy+1)*nx+2*ix]);
                    
                    elem_p->setNode(4, *(*nodes)[2*iy*nx+2*ix+1]);
                    elem_p->setNode(5, *(*nodes)[(2*iy+1)*nx+2*(ix+1)]);
                    elem_p->setNode(6, *(*nodes)[2*(iy+1)*nx+2*ix+1]);
                    elem_p->setNode(7, *(*nodes)[(2*iy+1)*nx+2*ix]);
                    
                    elem_p->setNode(8, *(*nodes)[(2*iy+1)*nx+(2*ix+1)]);
                    
                    elem_p->setExternalID(id);
                    id++;
                }
        }
            break;
            
            
        case FESystem::Mesh::TRI3:
        {
            switch (m_type) {
                case RIGHT_DIAGONAL:
                case LEFT_DIAGONAL:
                {
                    n_elems=(nx-1)*(ny-1)*2;
                    n_elem_nodes = 3;
                    
                    elems.reset(mesh.createElements(n_elems, elem_type, local_cs_same_as_global).release());
                    
                    // set the location of individual nodes
                    id=0;
                    FESystem::Mesh::ElemBase* elem_p;
                    for (FESystemUInt iy=0; iy<(ny-1); iy++)
                        for (FESystemUInt ix=0; ix<(nx-1); ix++)
                            if (m_type == RIGHT_DIAGONAL)
                            {
                                elem_p = (*elems)[id];
                                elem_p->setNode(0, *(*nodes)[(iy*nx+ix)]);
                                elem_p->setNode(1, *(*nodes)[(iy*nx+ix+1)]);
                                elem_p->setNode(2, *(*nodes)[((iy+1)*nx+ix+1)]);
                                elem_p->setExternalID(id);
                                id++;
                                
                                elem_p = (*elems)[id];
                                elem_p->setNode(0, *(*nodes)[((iy+1)*nx+ix+1)]);
                                elem_p->setNode(1, *(*nodes)[((iy+1)*nx+ix)]);
                                elem_p->setNode(2, *(*nodes)[(iy*nx+ix)]);
                                elem_p->setExternalID(id);
                                id++;
                            }
                            else // LEFT_DIAGONAL
                            {
                                elem_p = (*elems)[id];
                                elem_p->setNode(0, *(*nodes)[(iy*nx+ix)]);
                                elem_p->setNode(1, *(*nodes)[(iy*nx+ix+1)]);
                                elem_p->setNode(2, *(*nodes)[((iy+1)*nx+ix)]);
                                elem_p->setExternalID(id);
                                id++;
                                
                                elem_p = (*elems)[id];
                                elem_p->setNode(0, *(*nodes)[(iy*nx+ix+1)]);
                                elem_p->setNode(1, *(*nodes)[((iy+1)*nx+ix+1)]);
                                elem_p->setNode(2, *(*nodes)[((iy+1)*nx+ix)]);
                                elem_p->setExternalID(id);
                                id++;
                            }
                }
                    break;
                    
                case CROSS:
                {
                    FESystemAssert0(fmod(1.0*(nx-1),2.0) == 0.0, FESystem::Exception::InvalidValue);
                    FESystemAssert0(fmod(1.0*(ny-1),2.0) == 0.0, FESystem::Exception::InvalidValue);

                    n_elems=(nx-1)*(ny-1)*4;
                    n_elem_nodes = 3;

                    std::auto_ptr<std::vector<FESystem::Mesh::Node*> > mid_nodes;
                    mid_nodes.reset(mesh.createNodes(n_elems, global_cs).release());
                    
                    elems.reset(mesh.createElements(n_elems, elem_type, local_cs_same_as_global).release());
                    
                    // set the location of individual nodes
                    id=0;
                    FESystem::Mesh::ElemBase* elem_p;
                    for (FESystemUInt iy=0; iy<(ny-1); iy++)
                        for (FESystemUInt ix=0; ix<(nx-1); ix++)
                        {
                            (*mid_nodes)[iy*(nx-1)+ix]->setVal(0, 0.5*((*nodes)[iy*nx+ix]->getVal(0)+(*nodes)[iy*nx+ix+1]->getVal(0)));
                            (*mid_nodes)[iy*(nx-1)+ix]->setVal(1, 0.5*((*nodes)[iy*nx+ix]->getVal(1)+(*nodes)[(iy+1)*nx+ix]->getVal(1)));
                            
                            elem_p = (*elems)[id];
                            elem_p->setNode(0, *(*nodes)[iy*nx+ix]);
                            elem_p->setNode(1, *(*nodes)[iy*nx+ix+1]);
                            elem_p->setNode(2, *(*mid_nodes)[iy*(nx-1)+ix]);
                            elem_p->setExternalID(id);
                            id++;
                            
                            elem_p = (*elems)[id];
                            elem_p->setNode(0, *(*nodes)[iy*nx+ix+1]);
                            elem_p->setNode(1, *(*nodes)[(iy+1)*nx+ix+1]);
                            elem_p->setNode(2, *(*mid_nodes)[iy*(nx-1)+ix]);
                            elem_p->setExternalID(id);
                            id++;
                            
                            elem_p = (*elems)[id];
                            elem_p->setNode(0, *(*nodes)[(iy+1)*nx+ix+1]);
                            elem_p->setNode(1, *(*nodes)[(iy+1)*nx+ix]);
                            elem_p->setNode(2, *(*mid_nodes)[iy*(nx-1)+ix]);
                            elem_p->setExternalID(id);
                            id++;
                            
                            elem_p = (*elems)[id];
                            elem_p->setNode(0, *(*nodes)[(iy+1)*nx+ix]);
                            elem_p->setNode(1, *(*nodes)[iy*nx+ix]);
                            elem_p->setNode(2, *(*mid_nodes)[iy*(nx-1)+ix]);
                            elem_p->setExternalID(id);
                            id++;
                        }
                }
                    break;
                    
                default:
                    FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
                    break;
            }
            
            
        }
            break;
            
            
        case FESystem::Mesh::TRI6:
        {
            FESystemAssert0(fmod(1.0*(nx-1),2.0) == 0.0, FESystem::Exception::InvalidValue);
            FESystemAssert0(fmod(1.0*(ny-1),2.0) == 0.0, FESystem::Exception::InvalidValue);

            FESystemUInt nx_elems = (nx-1)/2, ny_elems = (ny-1)/2;
            
            n_elems=(nx_elems-1)*(ny_elems-1)*4;
            n_elem_nodes = 6;
            
            elems.reset(mesh.createElements(n_elems, elem_type, local_cs_same_as_global).release());
            
            // set the location of individual nodes
            id=0;
            FESystem::Mesh::ElemBase* elem_p;
            for (FESystemUInt iy=0; iy<(ny-1); iy++)
                for (FESystemUInt ix=0; ix<(nx-1); ix++)
                {
                    elem_p = (*elems)[id];
                    elem_p->setNode(0, *(*nodes)[((2*iy)*(2*nx_elems-1)+2*ix)]);
                    elem_p->setNode(1, *(*nodes)[((2*iy)*(2*nx_elems-1)+2*(ix+1))]);
                    elem_p->setNode(2, *(*nodes)[((2*iy+1)*(2*nx_elems-1)+(2*ix+1))]);
                    
                    elem_p->setNode(3, *(*nodes)[(2*iy*(2*nx_elems-1)+2*ix+1)]);
                    elem_p->setNode(4, *(*nodes)[((2*nx_elems-1)*(2*ny-1)+4*(nx_elems-1)*iy+4*ix+1)]);
                    elem_p->setNode(5, *(*nodes)[((2*nx_elems-1)*(2*ny-1)+4*(nx_elems-1)*iy+4*ix+0)]);
                    
                    elem_p->setExternalID(id);
                    id++;
                    
                    elem_p = (*elems)[id];
                    elem_p->setNode(0, *(*nodes)[((2*iy)*(2*nx_elems-1)+2*(ix+1))]);
                    elem_p->setNode(1, *(*nodes)[(2*(iy+1)*(2*nx_elems-1)+2*(ix+1))]);
                    elem_p->setNode(2, *(*nodes)[((2*iy+1)*(2*nx_elems-1)+(2*ix+1))]);
                    
                    elem_p->setNode(3, *(*nodes)[((2*iy+1)*(2*nx_elems-1)+2*(ix+1))]);
                    elem_p->setNode(4, *(*nodes)[((2*nx_elems-1)*(2*ny-1)+4*(nx_elems-1)*iy+4*ix+2)]);
                    elem_p->setNode(5, *(*nodes)[((2*nx_elems-1)*(2*ny-1)+4*(nx_elems-1)*iy+4*ix+1)]);
                    
                    elem_p->setExternalID(id);
                    id++;
                    
                    elem_p = (*elems)[id];
                    elem_p->setNode(0, *(*nodes)[(2*(iy+1)*(2*nx_elems-1)+2*(ix+1))]);
                    elem_p->setNode(1, *(*nodes)[(2*(iy+1)*(2*nx_elems-1)+2*ix)]);
                    elem_p->setNode(2, *(*nodes)[((2*iy+1)*(2*nx_elems-1)+(2*ix+1))]);
                    
                    elem_p->setNode(3, *(*nodes)[(2*(iy+1)*(2*nx_elems-1)+2*ix+1)]);
                    elem_p->setNode(4, *(*nodes)[((2*nx_elems-1)*(2*ny-1)+4*(nx_elems-1)*iy+4*ix+3)]);
                    elem_p->setNode(5, *(*nodes)[((2*nx_elems-1)*(2*ny-1)+4*(nx_elems-1)*iy+4*ix+2)]);
                    
                    elem_p->setExternalID(id);
                    id++;
                    
                    elem_p = (*elems)[id];
                    elem_p->setNode(0, *(*nodes)[(2*(iy+1)*(2*nx_elems-1)+2*ix)]);
                    elem_p->setNode(1, *(*nodes)[(2*iy*(2*nx_elems-1)+2*ix)]);
                    elem_p->setNode(2, *(*nodes)[((2*iy+1)*(2*nx_elems-1)+(2*ix+1))]);
                    
                    elem_p->setNode(3, *(*nodes)[((2*iy+1)*(2*nx_elems-1)+2*ix)]);
                    elem_p->setNode(4, *(*nodes)[((2*nx_elems-1)*(2*ny-1)+4*(nx_elems-1)*iy+4*ix+0)]);
                    elem_p->setNode(5, *(*nodes)[((2*nx_elems-1)*(2*ny-1)+4*(nx_elems-1)*iy+4*ix+3)]);
                    
                    elem_p->setExternalID(id);
                    id++;
                }
        }
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
    }
    mesh.reinit();
    
}




