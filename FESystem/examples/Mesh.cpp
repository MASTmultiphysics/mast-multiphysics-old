//
//  MeshTest.cpp
//  FESystemApplication
//
//  Created by Manav Bhatia on 3/19/12.
//  Copyright (c) 2012. All rights reserved.
//

// FESystem include
#include "TestingIncludes.h"






void createLineMesh(FESystem::Mesh::ElementType elem_type, FESystem::Mesh::MeshBase& mesh, FESystem::Geometry::Point& origin,
                    FESystemUInt nx, FESystemDouble x_length, FESystemUInt& n_elem_nodes, MeshType m_type, FESystemBoolean local_cs_same_as_global)
{
    // create a nx x ny grid of nodes, and connect them by quad4 elements
    FESystemDouble dx;
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
            n_nodes=nx;
            n_elems=(nx-1);
            n_elem_nodes = 2;
            dx=x_length/(nx-1);
            
            nodes.reset(mesh.createNodes(n_nodes, global_cs).release());
            elems.reset(mesh.createElements(n_elems, elem_type, local_cs_same_as_global).release());
            
            FESystemUInt id=0;
            FESystem::Mesh::Node* node_p;
            for (FESystemUInt ix=0; ix<nx; ix++)
            {
                node_p = (*nodes)[id];
                node_p->setVal(0,ix*dx);
                node_p->setExternalID(id);
                id++;
            }
            
            
            // set the location of individual nodes
            id=0;
            FESystem::Mesh::ElemBase* elem_p;
            for (FESystemUInt ix=0; ix<(nx-1); ix++)
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
            n_nodes=(2*nx-1);
            n_elems=(nx-1);
            n_elem_nodes = 3;
            dx=x_length/(2*(nx-1));
            
            nodes.reset(mesh.createNodes(n_nodes, global_cs).release());
            elems.reset(mesh.createElements(n_elems, elem_type, local_cs_same_as_global).release());
            
            FESystemUInt id=0;
            FESystem::Mesh::Node* node_p;
            for (FESystemUInt ix=0; ix<(2*nx-1); ix++)
            {
                node_p = (*nodes)[id];
                node_p->setVal(0,ix*dx);
                node_p->setExternalID(id);
                id++;
            }
            
            
            // set the location of individual nodes
            id=0;
            FESystem::Mesh::ElemBase* elem_p;
            for (FESystemUInt ix=0; ix<(nx-1); ix++)
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
                       FESystemUInt nx, FESystemUInt ny, FESystemDouble x_length, FESystemDouble y_length, FESystemUInt& n_elem_nodes,
                       MeshType m_type, FESystemBoolean local_cs_same_as_global)
{
    // create a nx x ny grid of nodes, and connect them by quad4 elements
    FESystemDouble dx, dy;
    FESystemUInt n_nodes, n_elems;
    
    const FESystem::Geometry::CoordinateSystemBase& global_cs = origin.getCoordinateSystem();
    
    std::auto_ptr<std::vector<FESystem::Mesh::Node*> > nodes;
    std::auto_ptr<std::vector<FESystem::Mesh::ElemBase*> > elems;
    
    // set the location of individual nodes
    switch (elem_type)
    {
        case FESystem::Mesh::QUAD4:
        {
            n_nodes=nx*ny;
            n_elems=(nx-1)*(ny-1);
            n_elem_nodes = 4;
            dx=x_length/(nx-1);
            dy=y_length/(ny-1);
            
            nodes.reset(mesh.createNodes(n_nodes, global_cs).release());
            elems.reset(mesh.createElements(n_elems, elem_type, local_cs_same_as_global).release());
            
            FESystemUInt id=0;
            FESystem::Mesh::Node* node_p;
            for (FESystemUInt iy=0; iy<ny; iy++)
                for (FESystemUInt ix=0; ix<nx; ix++)
                {
                    node_p = (*nodes)[id];
                    node_p->setVal(0,ix*dx);
                    node_p->setVal(1,iy*dy);
                    node_p->setExternalID(id);
                    id++;
                }
            
            
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
            n_nodes=(2*nx-1)*(2*ny-1);
            n_elems=(nx-1)*(ny-1);
            n_elem_nodes = 9;
            dx=x_length/(2*(nx-1));
            dy=y_length/(2*(ny-1));
            
            nodes.reset(mesh.createNodes(n_nodes, global_cs).release());
            elems.reset(mesh.createElements(n_elems, elem_type, local_cs_same_as_global).release());
            
            FESystemUInt id=0;
            FESystem::Mesh::Node* node_p;
            for (FESystemUInt iy=0; iy<(2*ny-1); iy++)
                for (FESystemUInt ix=0; ix<(2*nx-1); ix++)
                {
                    node_p = (*nodes)[id];
                    node_p->setVal(0,ix*dx);
                    node_p->setVal(1,iy*dy);
                    node_p->setExternalID(id);
                    id++;
                }
            
            
            // set the location of individual nodes
            id=0;
            FESystem::Mesh::ElemBase* elem_p;
            for (FESystemUInt iy=0; iy<(ny-1); iy++)
                for (FESystemUInt ix=0; ix<(nx-1); ix++)
                {
                    elem_p = (*elems)[id];
                    elem_p->setNode(0, *(*nodes)[((2*iy)*(2*nx-1)+2*ix)]);
                    elem_p->setNode(1, *(*nodes)[((2*iy)*(2*nx-1)+2*(ix+1))]);
                    elem_p->setNode(2, *(*nodes)[(2*(iy+1)*(2*nx-1)+2*(ix+1))]);
                    elem_p->setNode(3, *(*nodes)[(2*(iy+1)*(2*nx-1)+2*ix)]);
                    
                    elem_p->setNode(4, *(*nodes)[(2*iy*(2*nx-1)+2*ix+1)]);
                    elem_p->setNode(5, *(*nodes)[((2*iy+1)*(2*nx-1)+2*(ix+1))]);
                    elem_p->setNode(6, *(*nodes)[(2*(iy+1)*(2*nx-1)+2*ix+1)]);
                    elem_p->setNode(7, *(*nodes)[((2*iy+1)*(2*nx-1)+2*ix)]);
                    
                    elem_p->setNode(8, *(*nodes)[((2*iy+1)*(2*nx-1)+(2*ix+1))]);
                    
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
                    n_nodes=nx*ny;
                    n_elems=(nx-1)*(ny-1)*2;
                    n_elem_nodes = 3;
                    dx=x_length/(nx-1);
                    dy=y_length/(ny-1);
                    
                    nodes.reset(mesh.createNodes(n_nodes, global_cs).release());
                    elems.reset(mesh.createElements(n_elems, elem_type, local_cs_same_as_global).release());
                    
                    FESystemUInt id=0;
                    FESystem::Mesh::Node* node_p;
                    for (FESystemUInt iy=0; iy<ny; iy++)
                        for (FESystemUInt ix=0; ix<nx; ix++)
                        {
                            node_p = (*nodes)[id];
                            node_p->setVal(0,ix*dx);
                            node_p->setVal(1,iy*dy);
                            node_p->setExternalID(id);
                            id++;
                        }
                    
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
                    n_nodes=nx*ny + (nx-1)*(ny-1);
                    n_elems=(nx-1)*(ny-1)*4;
                    n_elem_nodes = 3;
                    dx=x_length/(nx-1);
                    dy=y_length/(ny-1);
                    
                    nodes.reset(mesh.createNodes(n_nodes, global_cs).release());
                    elems.reset(mesh.createElements(n_elems, elem_type, local_cs_same_as_global).release());
                    
                    FESystemUInt id=0;
                    FESystem::Mesh::Node* node_p;
                    for (FESystemUInt iy=0; iy<ny; iy++)
                        for (FESystemUInt ix=0; ix<nx; ix++)
                        {
                            node_p = (*nodes)[id];
                            node_p->setVal(0,ix*dx);
                            node_p->setVal(1,iy*dy);
                            node_p->setExternalID(id);
                            id++;
                        }
                    for (FESystemUInt iy=0; iy<ny-1; iy++)
                        for (FESystemUInt ix=0; ix<nx-1; ix++)
                        {
                            node_p = (*nodes)[id];
                            node_p->setVal(0,ix*dx+dx/2);
                            node_p->setVal(1,iy*dy+dy/2);
                            node_p->setExternalID(id);
                            id++;
                        }
                    
                    // set the location of individual nodes
                    id=0;
                    FESystem::Mesh::ElemBase* elem_p;
                    for (FESystemUInt iy=0; iy<(ny-1); iy++)
                        for (FESystemUInt ix=0; ix<(nx-1); ix++)
                        {
                            elem_p = (*elems)[id];
                            elem_p->setNode(0, *(*nodes)[(iy*nx+ix)]);
                            elem_p->setNode(1, *(*nodes)[(iy*nx+ix+1)]);
                            elem_p->setNode(2, *(*nodes)[nx*ny+iy*(nx-1)+ix]);
                            elem_p->setExternalID(id);
                            id++;
                            
                            elem_p = (*elems)[id];
                            elem_p->setNode(0, *(*nodes)[(iy*nx+ix+1)]);
                            elem_p->setNode(1, *(*nodes)[((iy+1)*nx+ix+1)]);
                            elem_p->setNode(2, *(*nodes)[nx*ny+iy*(nx-1)+ix]);
                            elem_p->setExternalID(id);
                            id++;
                            
                            elem_p = (*elems)[id];
                            elem_p->setNode(0, *(*nodes)[((iy+1)*nx+ix+1)]);
                            elem_p->setNode(1, *(*nodes)[((iy+1)*nx+ix)]);
                            elem_p->setNode(2, *(*nodes)[nx*ny+iy*(nx-1)+ix]);
                            elem_p->setExternalID(id);
                            id++;
                            
                            elem_p = (*elems)[id];
                            elem_p->setNode(0, *(*nodes)[((iy+1)*nx+ix)]);
                            elem_p->setNode(1, *(*nodes)[(iy*nx+ix)]);
                            elem_p->setNode(2, *(*nodes)[nx*ny+iy*(nx-1)+ix]);
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
            n_nodes=(2*nx-1)*(2*ny-1) + 4*(nx-1)*(ny-1);
            n_elems=(nx-1)*(ny-1)*4;
            n_elem_nodes = 6;
            dx=x_length/(2*(nx-1));
            dy=y_length/(2*(ny-1));
            
            nodes.reset(mesh.createNodes(n_nodes, global_cs).release());
            elems.reset(mesh.createElements(n_elems, elem_type, local_cs_same_as_global).release());
            
            FESystemUInt id=0;
            FESystem::Mesh::Node* node_p;
            for (FESystemUInt iy=0; iy<(2*ny-1); iy++)
                for (FESystemUInt ix=0; ix<(2*nx-1); ix++)
                {
                    node_p = (*nodes)[id];
                    node_p->setVal(0,ix*dx);
                    node_p->setVal(1,iy*dy);
                    node_p->setExternalID(id);
                    id++;
                }
            for (FESystemUInt iy=0; iy<(ny-1); iy++)
                for (FESystemUInt ix=0; ix<(nx-1); ix++)
                {
                    node_p = (*nodes)[id];
                    node_p->setVal(0,ix*dx*2+dx/2.0);
                    node_p->setVal(1,iy*dy*2+dy/2.0);
                    node_p->setExternalID(id);
                    id++;
                    
                    node_p = (*nodes)[id];
                    node_p->setVal(0,ix*dx*2+dx/2.0+dx);
                    node_p->setVal(1,iy*dy*2+dy/2.0);
                    node_p->setExternalID(id);
                    id++;
                    
                    node_p = (*nodes)[id];
                    node_p->setVal(0,ix*dx*2+dx/2.0+dx);
                    node_p->setVal(1,iy*dy*2+dy/2.0+dy);
                    node_p->setExternalID(id);
                    id++;
                    
                    node_p = (*nodes)[id];
                    node_p->setVal(0,ix*dx*2+dx/2.0);
                    node_p->setVal(1,iy*dy*2+dy/2.0+dy);
                    node_p->setExternalID(id);
                    id++;
                    
                }
            
            
            // set the location of individual nodes
            id=0;
            FESystem::Mesh::ElemBase* elem_p;
            for (FESystemUInt iy=0; iy<(ny-1); iy++)
                for (FESystemUInt ix=0; ix<(nx-1); ix++)
                {
                    elem_p = (*elems)[id];
                    elem_p->setNode(0, *(*nodes)[((2*iy)*(2*nx-1)+2*ix)]);
                    elem_p->setNode(1, *(*nodes)[((2*iy)*(2*nx-1)+2*(ix+1))]);
                    elem_p->setNode(2, *(*nodes)[((2*iy+1)*(2*nx-1)+(2*ix+1))]);
                    
                    elem_p->setNode(3, *(*nodes)[(2*iy*(2*nx-1)+2*ix+1)]);
                    elem_p->setNode(4, *(*nodes)[((2*nx-1)*(2*ny-1)+4*(nx-1)*iy+4*ix+1)]);
                    elem_p->setNode(5, *(*nodes)[((2*nx-1)*(2*ny-1)+4*(nx-1)*iy+4*ix+0)]);
                    
                    elem_p->setExternalID(id);
                    id++;
                    
                    elem_p = (*elems)[id];
                    elem_p->setNode(0, *(*nodes)[((2*iy)*(2*nx-1)+2*(ix+1))]);
                    elem_p->setNode(1, *(*nodes)[(2*(iy+1)*(2*nx-1)+2*(ix+1))]);
                    elem_p->setNode(2, *(*nodes)[((2*iy+1)*(2*nx-1)+(2*ix+1))]);
                    
                    elem_p->setNode(3, *(*nodes)[((2*iy+1)*(2*nx-1)+2*(ix+1))]);
                    elem_p->setNode(4, *(*nodes)[((2*nx-1)*(2*ny-1)+4*(nx-1)*iy+4*ix+2)]);
                    elem_p->setNode(5, *(*nodes)[((2*nx-1)*(2*ny-1)+4*(nx-1)*iy+4*ix+1)]);
                    
                    elem_p->setExternalID(id);
                    id++;
                    
                    elem_p = (*elems)[id];
                    elem_p->setNode(0, *(*nodes)[(2*(iy+1)*(2*nx-1)+2*(ix+1))]);
                    elem_p->setNode(1, *(*nodes)[(2*(iy+1)*(2*nx-1)+2*ix)]);
                    elem_p->setNode(2, *(*nodes)[((2*iy+1)*(2*nx-1)+(2*ix+1))]);
                    
                    elem_p->setNode(3, *(*nodes)[(2*(iy+1)*(2*nx-1)+2*ix+1)]);
                    elem_p->setNode(4, *(*nodes)[((2*nx-1)*(2*ny-1)+4*(nx-1)*iy+4*ix+3)]);
                    elem_p->setNode(5, *(*nodes)[((2*nx-1)*(2*ny-1)+4*(nx-1)*iy+4*ix+2)]);
                    
                    elem_p->setExternalID(id);
                    id++;
                    
                    elem_p = (*elems)[id];
                    elem_p->setNode(0, *(*nodes)[(2*(iy+1)*(2*nx-1)+2*ix)]);
                    elem_p->setNode(1, *(*nodes)[(2*iy*(2*nx-1)+2*ix)]);
                    elem_p->setNode(2, *(*nodes)[((2*iy+1)*(2*nx-1)+(2*ix+1))]);
                    
                    elem_p->setNode(3, *(*nodes)[((2*iy+1)*(2*nx-1)+2*ix)]);
                    elem_p->setNode(4, *(*nodes)[((2*nx-1)*(2*ny-1)+4*(nx-1)*iy+4*ix+0)]);
                    elem_p->setNode(5, *(*nodes)[((2*nx-1)*(2*ny-1)+4*(nx-1)*iy+4*ix+3)]);
                    
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




