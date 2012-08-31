//
//  MeshTest.cpp
//  FESystemApplication
//
//  Created by Manav Bhatia on 3/19/12.
//  Copyright (c) 2012. All rights reserved.
//

// C++ includes
#include <iostream>
#include <fstream>
#include <iomanip>
#include <memory.h>
#include <cassert>

// FESystem includes
#include "Base/FESystemBase.h"
#include "Utils/AutoPtrTable.h"
// Plotting
#include "Plotting/PLPlot.h"
// geometry
#include "Geom/Point.h"
#include "Geom/RectangularCoordinateSystem.h"
// mesh
#include "Mesh/MeshBase.h"
#include "Mesh/ElementType.h"
#include "Mesh/Node.h"
#include "Mesh/EdgeElemBase.h"
#include "Mesh/Edge2.h"
#include "Mesh/Edge3.h"
#include "Mesh/Quad4.h"
#include "Mesh/Quad9.h"
#include "Mesh/Quad8.h"
#include "Mesh/Tri3.h"
#include "Mesh/Tri6.h"
#include "Mesh/Tri7.h"
// finite elements
#include "FiniteElems/FELagrange.h"
#include "Functions/FunctionMappingBase.h"
// quadrature
#include "Quadrature/TrapezoidQuadrature.h"
// degrees of freedom
#include "Base/DegreeOfFreedomMap.h"
#include "Base/DegreeOfFreedomUnit.h"
// numerics
#include "Numerics/DenseMatrix.h"
#include "Numerics/SparseMatrix.h"
#include "Numerics/SparsityPattern.h"
// solver
#include "Solvers/LinearSolvers/LapackLinearSolver.h"
#include "Solvers/LinearSolvers/LUFactorizationLinearSolver.h"
#include "Solvers/LinearSolvers/PseudoTimeSteppingLinearSolver.h"
#include "Solvers/LinearSolvers/QRFactorizationLinearSolver.h"
#include "Solvers/EigenSolvers/QRMethodLinearEigensolver.h"
#include "Solvers/EigenSolvers/LapackLinearEigenSolver.h"
#include "Solvers/EigenSolvers/ArpackLinearEigenSolver.h"
#include "Solvers/TransientSolvers/NewmarkTransientSolver.h"
#include "Solvers/TransientSolvers/ExplicitRungeKuttaTransientSolver.h"
#include "Solvers/NonlinearSolvers/NewtonIterationNonlinearSolver.h"
// output processors
#include "OutputProcessors/GmshOutputProcessor.h"
#include "OutputProcessors/VtkOutputProcessor.h"
#include "OutputProcessors/TecplotOutputProcessor.h"
// Structural elements
#include "Disciplines/Structure/ReissnerMindlinPlate.h"
#include "Disciplines/Structure/DKTPlate.h"
#include "Disciplines/Structure/EulerBernoulliBeam.h"
#include "Disciplines/Structure/TimoshenkoBeam.h"
#include "Disciplines/Structure/ExtensionBar.h"
#include "Disciplines/Structure/VonKarmanStrain1D.h"


enum MeshType{RIGHT_DIAGONAL, LEFT_DIAGONAL, CROSS, INVALID_MESH};




int shape_function_main1D(int argc, char * const argv[])
{
    FESystemUInt dim = 1, n_points=15;
    
    //FESystem::Plotting::PLPlot<FESystemDouble> plot(FESystem::Plotting::REAL_AXIS, FESystem::Plotting::REAL_AXIS);
    
    FESystem::Geometry::Point origin(3);
    const FESystem::Geometry::CoordinateSystemBase global_cs = origin.getCoordinateSystem();
    FESystem::Numerics::LocalVector<FESystemDouble> y_unit;
    y_unit.resize(3);
    y_unit.setVal(1, 1.0);
    
    // create the finite element and initialize the shape functions
    FESystem::FiniteElement::FELagrange fe;
    FESystem::Geometry::Point p(dim);
    
    // now start to integrate the stiffness matrix
    std::auto_ptr<FESystem::Mesh::ElemBase> elem(new FESystem::Mesh::Edge3);
    std::auto_ptr<FESystem::Mesh::Node> node0(new FESystem::Mesh::Node(global_cs)), node1(new FESystem::Mesh::Node(global_cs)), node2(new FESystem::Mesh::Node(global_cs));
    dynamic_cast<FESystem::Mesh::EdgeElemBase*>(elem.get())->setVectorForXYPlane(y_unit);
    node0->setVal(0, -1); elem->setNode(0, *node0);
    node1->setVal(0,  1); elem->setNode(1, *node1);
    node2->setVal(0,  0); elem->setNode(2, *node2);
    
    std::vector<FESystemUInt> derivatives(dim);
    derivatives[0] = 1;
    
    FESystem::Numerics::DenseMatrix<FESystemDouble> jac_mat;
    FESystem::Numerics::LocalVector<FESystemDouble> N,dNdx, x_vals, y_vals1, y_vals2, y_vals3, dydx_vals1, dydx_vals2, dydx_vals3;
    jac_mat.resize(1,1);
    N.resize(elem->getNNodes()); dNdx.resize(elem->getNNodes());
    x_vals.resize(n_points), y_vals1.resize(n_points); y_vals2.resize(n_points); y_vals3.resize(n_points), dydx_vals1.resize(n_points); dydx_vals2.resize(n_points); dydx_vals3.resize(n_points);
    
    
    fe.reinit(*elem);
    for (FESystemUInt i=0; i<n_points; i++)
    {
        p.setVal(0, -1.0+2.0/(n_points-1)*i);
        fe.getShapeFunction(p,N);
        fe.getShapeFunctionDerivativeForLocalCoordinates(derivatives, p, dNdx);
        
        x_vals.setVal(i, p.getVal(0));
        y_vals1.setVal(i, N.getVal(0));
        y_vals2.setVal(i, N.getVal(1));
        y_vals3.setVal(i, N.getVal(2));
        dydx_vals1.setVal(i, dNdx.getVal(0));
        dydx_vals2.setVal(i, dNdx.getVal(1));
        dydx_vals3.setVal(i, dNdx.getVal(2));
    }
    
    
    //    plot.plotData2D(x_vals, y_vals1);
    //    plot.plotData2D(x_vals, y_vals2);
    //    plot.plotData2D(x_vals, y_vals3);
    //    plot.plotData2D(x_vals, dydx_vals1);
    //    plot.plotData2D(x_vals, dydx_vals2);
    //    plot.plotData2D(x_vals, dydx_vals3);
    
    //    std::cout << "Node Locations: " << std::endl;
    //    for (FESystemInt j=0; j<elem->getNNodes(); j++)
    //        elem->getNode(j).write(std::cout);
    //    std::cout << "Qpt Location: " << std::endl;
    //    p.write(std::cout);
    //    std::cout << "Shape Funcs: " << std::endl; N.write(std::cout);
    //    std::cout << "Shape Derivs: " << std::endl; dNdx.write(std::cout);
    
    return 0;
}



int shape_function_main2DQuad(int argc, char * const argv[])
{
    FESystemUInt dim = 2, n_point=30;
    
    //FESystem::Plotting::PLPlot<FESystemDouble> plot(FESystem::Plotting::REAL_AXIS, FESystem::Plotting::REAL_AXIS);
    
    FESystem::Geometry::Point origin(3);
    const FESystem::Geometry::CoordinateSystemBase global_cs = origin.getCoordinateSystem();
    
    // create the finite element and initialize the shape functions
    FESystem::FiniteElement::FELagrange fe;
    FESystem::Geometry::Point p(dim);
    
    // now start to integrate the stiffness matrix
    std::auto_ptr<FESystem::Mesh::ElemBase> elem(new FESystem::Mesh::Quad4);
    std::auto_ptr<FESystem::Mesh::Node>
    node0(new FESystem::Mesh::Node(global_cs)), node1(new FESystem::Mesh::Node(global_cs)), node2(new FESystem::Mesh::Node(global_cs)),
    node3(new FESystem::Mesh::Node(global_cs)), node4(new FESystem::Mesh::Node(global_cs)), node5(new FESystem::Mesh::Node(global_cs)),
    node6(new FESystem::Mesh::Node(global_cs)), node7(new FESystem::Mesh::Node(global_cs)), node8(new FESystem::Mesh::Node(global_cs));
    
    node0->setVal(0, -1); node0->setVal(1, -1); elem->setNode(0, *node0);
    node1->setVal(0,  1); node1->setVal(1, -1); elem->setNode(1, *node1);
    node2->setVal(0,  1); node2->setVal(1,  1); elem->setNode(2, *node2);
    node3->setVal(0, -1); node3->setVal(1,  1); elem->setNode(3, *node3);
    
    //    node4->setVal(0,  0); node4->setVal(1, -1); elem->setNode(4, *node4);
    //    node5->setVal(0,  1); node5->setVal(1,  0); elem->setNode(5, *node5);
    //    node6->setVal(0,  0); node6->setVal(1,  1); elem->setNode(6, *node6);
    //    node7->setVal(0, -1); node7->setVal(1,  0); elem->setNode(7, *node7);
    //
    //    node0->setVal(0, -1); node0->setVal(1, -1); elem->setNode(0, *node0);
    //    node1->setVal(0,  1); node1->setVal(1, -1); elem->setNode(1, *node1);
    //    node2->setVal(0, -1); node2->setVal(1,  1); elem->setNode(2, *node2);
    //    node3->setVal(0, -1); node3->setVal(1,  1); elem->setNode(3, *node3);
    //
    //    node4->setVal(0,  0); node4->setVal(1, -1); elem->setNode(4, *node4);
    //    node5->setVal(0,  0); node5->setVal(1,  0); elem->setNode(5, *node5);
    //    node6->setVal(0, -1); node6->setVal(1,  1); elem->setNode(6, *node6);
    //    node7->setVal(0, -1); node7->setVal(1,  0); elem->setNode(7, *node7);
    //
    //    node8->setVal(0,-.5); node8->setVal(1,  0); elem->setNode(8, *node8);
    
    std::vector<FESystemUInt> derivatives(dim);
    derivatives[0] = 1; derivatives[1] = 0;
    
    FESystem::Numerics::DenseMatrix<FESystemDouble> jac_mat, z_vals;
    FESystem::Numerics::LocalVector<FESystemDouble> N,dNdx, x_vals, y_vals;
    jac_mat.resize(2,2); z_vals.resize(n_point, n_point);
    N.resize(elem->getNNodes()); dNdx.resize(elem->getNNodes());
    x_vals.resize(n_point); y_vals.resize(n_point);
    
    fe.reinit(*elem);
    for (FESystemUInt i=0; i<n_point; i++)
        for (FESystemUInt j=0; j<n_point; j++)
        {
            p.setVal(0,  -1.0+2.0/(n_point-1)*i); p.setVal(1,  -1.0+2.0/(n_point-1)*j);
            fe.getShapeFunction(p,N);
            fe.getJacobianMatrix(p, jac_mat);
            fe.getShapeFunctionDerivativeForLocalCoordinates(derivatives, p, dNdx);
            //fe.getShapeFunctionDerivativeForPhysicalCoordinates(derivatives, p, dNdx);
            x_vals.setVal(i, p.getVal(0)); y_vals.setVal(j, p.getVal(1));
            z_vals.setVal(i, j, dNdx.getVal(1));
        }
    
    //plot.plotData3DSurf(x_vals, y_vals, z_vals);
    
    return 0;
}




int shape_function_main2DTri(int argc, char * const argv[])
{
    FESystemUInt dim = 2, n_point=15;
    
    //FESystem::Plotting::PLPlot<FESystemDouble> plot(FESystem::Plotting::REAL_AXIS, FESystem::Plotting::REAL_AXIS);
    
    FESystem::Geometry::Point origin(3);
    const FESystem::Geometry::CoordinateSystemBase global_cs = origin.getCoordinateSystem();
    
    // create the finite element and initialize the shape functions
    FESystem::FiniteElement::FELagrange fe;
    FESystem::Geometry::Point p(dim);
    FESystem::Numerics::DenseMatrix<FESystemDouble> T_mat; T_mat.resize(3, 3);
    
    // now start to integrate the stiffness matrix
    std::auto_ptr<FESystem::Mesh::ElemBase> elem(new FESystem::Mesh::Tri6);
    std::auto_ptr<FESystem::Mesh::Node>
    node0(new FESystem::Mesh::Node(global_cs)), node1(new FESystem::Mesh::Node(global_cs)), node2(new FESystem::Mesh::Node(global_cs)),
    node3(new FESystem::Mesh::Node(global_cs)), node4(new FESystem::Mesh::Node(global_cs)), node5(new FESystem::Mesh::Node(global_cs)),
    node6(new FESystem::Mesh::Node(global_cs));
    
    node0->setVal(0, -1); node0->setVal(1, -1); elem->setNode(0, *node0);
    node1->setVal(0,  1); node1->setVal(1, -1); elem->setNode(1, *node1);
    node2->setVal(0, -1); node2->setVal(1,  1); elem->setNode(2, *node2);
    
    node3->setVal(0,  0); node3->setVal(1, -1); elem->setNode(3, *node3);
    node4->setVal(0,  0); node4->setVal(1,  0); elem->setNode(4, *node4);
    node5->setVal(0, -1); node5->setVal(1,  0); elem->setNode(5, *node5);
    
    //    node6->setVal(0,-.5); node6->setVal(1,  0); elem->setNode(6, *node6);
    
    //    elem->getLocalPhysicalCoordinateSystem().getFunctionMappingObject().getMappingJacobian(origin, T_mat);
    //    T_mat.write(std::cout);
    
    std::vector<FESystemUInt> derivatives(dim);
    derivatives[0] = 1; derivatives[1] = 0;
    
    FESystem::Numerics::DenseMatrix<FESystemDouble> jac_mat, z_vals;
    FESystem::Numerics::LocalVector<FESystemDouble> N, dNdx, x_vals, y_vals;
    jac_mat.resize(2,2); z_vals.resize(n_point, n_point);
    N.resize(elem->getNNodes()); dNdx.resize(elem->getNNodes());
    x_vals.resize(n_point); y_vals.resize(n_point);
    
    fe.reinit(*elem);
    for (FESystemUInt i=0; i<n_point; i++)
        for (FESystemUInt j=0; j<n_point; j++)
        {
            p.setVal(0,  -1.0+2.0/(n_point-1)*i); p.setVal(1,  -1.0+2.0/(n_point-1)*j);
            fe.getShapeFunction(p,N);
            fe.getJacobianMatrix(p, jac_mat);
            //fe.getShapeFunctionDerivativeForLocalCoordinates(derivatives, p, dNdx);
            fe.getShapeFunctionDerivativeForPhysicalCoordinates(derivatives, p, dNdx);
            x_vals.setVal(i, p.getVal(0)); y_vals.setVal(j, p.getVal(1));
            z_vals.setVal(i, j, dNdx.getVal(5));
        }
    
    //plot.plotData3DSurf(x_vals, y_vals, z_vals);
    
    return 0;
}







void createLineMesh(FESystem::Mesh::ElementType elem_type, FESystem::Mesh::MeshBase& mesh, FESystem::Geometry::Point& origin,
                    FESystemUInt nx, FESystemDouble x_length, FESystemUInt& n_elem_nodes, MeshType m_type=INVALID_MESH)
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
            elems.reset(mesh.createElements(n_elems, elem_type).release());
            
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
            elems.reset(mesh.createElements(n_elems, elem_type).release());
            
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




void create_plane_mesh(FESystem::Mesh::ElementType elem_type, FESystem::Mesh::MeshBase& mesh, FESystem::Geometry::Point& origin,
                       FESystemUInt nx, FESystemUInt ny, FESystemDouble x_length, FESystemDouble y_length, FESystemUInt& n_elem_nodes,
                       MeshType m_type=INVALID_MESH)
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
            elems.reset(mesh.createElements(n_elems, elem_type).release());
            
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
            elems.reset(mesh.createElements(n_elems, elem_type).release());
            
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
                    elems.reset(mesh.createElements(n_elems, elem_type).release());
                    
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
                    elems.reset(mesh.createElements(n_elems, elem_type).release());
                    
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
            elems.reset(mesh.createElements(n_elems, elem_type).release());
            
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


int linear_potential_flow_nonconservative(int argc, char * const argv[])
{
    FESystem::Mesh::ElementType elem_type = FESystem::Mesh::QUAD4;
    
    // create the mesh object
    FESystem::Mesh::MeshBase mesh;
    
    // create a nx x ny grid of nodes
    FESystemUInt nx=50, ny=40;
    FESystemDouble x_length = 20, y_length = 20, mach=1.85, sound_speed=330.0, chord=0.5, airfoil_t = chord*.1, x_airfoil_begin=x_length/2.0-chord/2.0;
    FESystemUInt dim = 2, n_elem_nodes, n_elem_dofs;
    
    FESystem::Geometry::Point origin(3);
    
    create_plane_mesh(elem_type, mesh, origin, nx, ny, x_length, y_length, n_elem_nodes);
    n_elem_dofs = n_elem_nodes;
    
    // set the location of individual nodes
    const std::vector<FESystem::Mesh::ElemBase*>& elems = mesh.getElements();
    
    // now add the degrees of freedom
    FESystem::Base::DegreeOfFreedomMap dof_map(mesh);
    std::string name;
    name = "phi"; dof_map.addVariable(name, 0);
    dof_map.reinit(); // distribute the dofs
    
    // create the finite element and initialize the shape functions
    FESystem::FiniteElement::FELagrange fe;
    FESystem::Quadrature::TrapezoidQuadrature q_rule, q_boundary_rule;
    q_rule.init(dim, 5);
    q_boundary_rule.init(dim-1, 5);
    
    // now start to integrate the stiffness matrix
    std::vector<FESystem::Mesh::ElemBase*>::const_iterator el_it=elems.begin(), el_end=elems.end();
    std::vector<FESystemUInt> derivatives(dim);
    
    FESystem::Numerics::DenseMatrix<FESystemDouble> el_mat, el_mat2, B_mat1, B_mat2, el_mat_combined, global_stiffness_mat, global_damping_mat, global_mass_mat, global_mass_mat_inv, tmp_mat;
    FESystem::Numerics::LocalVector<FESystemDouble> rhs, sol, tmp_vec,  dNdx, dNdy, Nfunc, Nvec, el_vec, point, elem_dof_vals;
    
    el_mat.resize(n_elem_dofs,n_elem_dofs); el_mat2.resize(n_elem_dofs, n_elem_dofs); el_mat_combined.resize(n_elem_dofs,n_elem_dofs);  el_vec.resize(n_elem_dofs); elem_dof_vals.resize(n_elem_dofs);
    tmp_mat.resize(dof_map.getNDofs(), dof_map.getNDofs());   global_mass_mat.resize(dof_map.getNDofs(), dof_map.getNDofs());
    global_stiffness_mat.resize(dof_map.getNDofs(), dof_map.getNDofs()); global_damping_mat.resize(dof_map.getNDofs(), dof_map.getNDofs()); global_mass_mat_inv.resize(dof_map.getNDofs(), dof_map.getNDofs());
    B_mat1.resize(1,n_elem_dofs); B_mat2.resize(1,n_elem_dofs); Nvec.resize(n_elem_nodes);
    rhs.resize(dof_map.getNDofs()); sol.resize(dof_map.getNDofs());    tmp_vec.resize(dof_map.getNDofs());
    dNdx.resize(n_elem_nodes);     dNdy.resize(n_elem_nodes);   Nfunc.resize(n_elem_nodes);
    point.resize(3);
    
    
    // get the values of shape functions and derivatives
    const std::vector<FESystem::Geometry::Point*>& q_pts = q_rule.getQuadraturePoints();
    const std::vector<FESystemDouble>& q_weight = q_rule.getQuadraturePointWeights();
    const std::vector<FESystem::Geometry::Point*>& q_pts_b = q_boundary_rule.getQuadraturePoints();
    const std::vector<FESystemDouble>& q_weight_b = q_boundary_rule.getQuadraturePointWeights();
    FESystemDouble x_val = 0.0, jac=0.0;
    
    for ( ; el_it != el_end; el_it++)
    {
        // initialize the finite element for this element
        fe.clear();
        fe.reinit(**el_it); // first order function with first order derivative
        
        // stiffness matrix
        for (FESystemUInt j=0; j<dim; j++)
        {
            el_mat_combined.zero();
            
            for (FESystemUInt k=0; k<derivatives.size(); k++) derivatives[k] = 0;
            derivatives[j] = 1; // derivative of the j^th coordinate
            
            for (FESystemUInt i=0; i<q_pts.size(); i++)
            {
                B_mat1.zero();
                el_mat.zero();
                
                fe.getShapeFunctionDerivativeForPhysicalCoordinates(derivatives, *(q_pts[i]), dNdx);
                jac = fe.getJacobianValue(*(q_pts[i]));
                
                B_mat1.setRowVals(0, 0, n_elem_nodes-1, dNdx); // epsilon-x: phix
                B_mat1.matrixTransposeRightMultiply(1.0, B_mat1, el_mat); // B^T B
                el_mat_combined.add(q_weight[i]*jac, el_mat);
            }
            if (j == 0) // for the x-axis
                el_mat_combined.scale(1.0-mach*mach);
            
            // now add this to the global matrix: stiffness matrix
            dof_map.addToGlobalMatrix(**el_it, el_mat_combined, global_stiffness_mat);
        }
        
        
        // nonreflecting boundary conditions
        // inflow on left edge
        if ((*el_it)->getNode(0).getVal(0) == 0.0) // left edge
        {
            (*el_it)->calculateBoundaryNormal(3, point);
            
            for (FESystemUInt j=0; j<dim; j++) // each boundary is scaled differently with the surface normal component
            {
                el_mat_combined.zero();
                
                for (FESystemUInt k=0; k<derivatives.size(); k++) derivatives[k] = 0;
                derivatives[j] = 1; // derivative of the x-coordinate
                
                for (FESystemUInt i=0; i<q_pts_b.size(); i++)
                {
                    B_mat1.zero();
                    B_mat2.zero();
                    el_mat.zero();
                    
                    fe.getShapeFunctionForBoundary(3, *(q_pts_b[i]), Nvec); // bottom edge is 0 boundary
                    fe.getShapeFunctionDerivativeForPhysicalCoordinatesForBoundary(3, derivatives, *(q_pts_b[i]), dNdx);
                    jac = fe.getJacobianValueForBoundary(3, *(q_pts_b[i]));
                    
                    B_mat1.setRowVals(0, 0, n_elem_dofs-1, Nvec);
                    B_mat2.setRowVals(0, 0, n_elem_dofs-1, dNdx);
                    B_mat1.matrixTransposeRightMultiply(1.0, B_mat2, el_mat);
                    
                    el_mat_combined.add(q_weight_b[i]*jac, el_mat);
                }
                if (j == 0)
                    el_mat_combined.scale(1.0-mach*mach);
                
                el_mat_combined.scale(-point.getVal(j));
                
                // now add this to the global matrix: stiffness matrix
                dof_map.addToGlobalMatrix(**el_it, el_mat_combined, global_stiffness_mat);
            }
        }
        // non-reflecting boundary conditions at the right edge
        else if ((*el_it)->getNode(1).getVal(0) == x_length) // right edge
        {
            (*el_it)->calculateBoundaryNormal(1, point);
            for (FESystemUInt j=0; j<dim; j++) // each boundary is scaled differently with the surface normal component
            {
                el_mat_combined.zero();
                
                for (FESystemUInt k=0; k<derivatives.size(); k++) derivatives[k] = 0;
                derivatives[j] = 1; // derivative of the x-coordinate
                
                for (FESystemUInt i=0; i<q_pts_b.size(); i++)
                {
                    B_mat1.zero();
                    B_mat2.zero();
                    el_mat.zero();
                    
                    fe.getShapeFunctionForBoundary(1, *(q_pts_b[i]), Nvec); // bottom edge is 0 boundary
                    fe.getShapeFunctionDerivativeForPhysicalCoordinatesForBoundary(1, derivatives, *(q_pts_b[i]), dNdx);
                    jac = fe.getJacobianValueForBoundary(1, *(q_pts_b[i]));
                    
                    B_mat1.setRowVals(0, 0, n_elem_dofs-1, Nvec);
                    B_mat2.setRowVals(0, 0, n_elem_dofs-1, dNdx);
                    B_mat1.matrixTransposeRightMultiply(1.0, B_mat2, el_mat);
                    
                    el_mat_combined.add(q_weight_b[i]*jac, el_mat);
                }
                if (j == 0)
                    el_mat_combined.scale(1.0-mach*mach);
                
                el_mat_combined.scale(-point.getVal(j));
                
                // now add this to the global matrix: stiffness matrix
                dof_map.addToGlobalMatrix(**el_it, el_mat_combined, global_stiffness_mat);
            }
        }
        else if ((*el_it)->getNode(2).getVal(1) == y_length) // top edge
        {
            (*el_it)->calculateBoundaryNormal(2, point);
            for (FESystemUInt j=0; j<dim; j++) // each boundary is scaled differently with the surface normal component
            {
                el_mat_combined.zero();
                
                for (FESystemUInt k=0; k<derivatives.size(); k++) derivatives[k] = 0;
                derivatives[j] = 1; // derivative of the x-coordinate
                
                for (FESystemUInt i=0; i<q_pts_b.size(); i++)
                {
                    B_mat1.zero();
                    B_mat2.zero();
                    el_mat.zero();
                    
                    fe.getShapeFunctionForBoundary(2, *(q_pts_b[i]), Nvec); // bottom edge is 0 boundary
                    fe.getShapeFunctionDerivativeForPhysicalCoordinatesForBoundary(2, derivatives, *(q_pts_b[i]), dNdx);
                    jac = fe.getJacobianValueForBoundary(2, *(q_pts_b[i]));
                    
                    B_mat1.setRowVals(0, 0, n_elem_dofs-1, Nvec);
                    B_mat2.setRowVals(0, 0, n_elem_dofs-1, dNdx);
                    B_mat1.matrixTransposeRightMultiply(1.0, B_mat2, el_mat);
                    
                    el_mat_combined.add(q_weight_b[i]*jac, el_mat);
                }
                if (j == 0)
                    el_mat_combined.scale(1.0-mach*mach);
                
                el_mat_combined.scale(-point.getVal(j));
                
                // now add this to the global matrix: stiffness matrix
                dof_map.addToGlobalMatrix(**el_it, el_mat_combined, global_stiffness_mat);
            }
        }
        
        
        
        // damping matrix
        el_mat_combined.zero();
        
        for (FESystemUInt k=0; k<derivatives.size(); k++) derivatives[k] = 0;
        derivatives[0] = 1; // derivative of the x-coordinate
        
        for (FESystemUInt i=0; i<q_pts.size(); i++)
        {
            B_mat1.zero();
            B_mat2.zero();
            el_mat.zero();
            
            fe.getShapeFunction(*(q_pts[i]), Nvec);
            fe.getShapeFunctionDerivativeForPhysicalCoordinates(derivatives, *(q_pts[i]), dNdx);
            jac = fe.getJacobianValue(*(q_pts[i]));
            
            B_mat1.setRowVals(0, 0, n_elem_nodes-1, dNdx); // epsilon-x: phix
            B_mat2.setRowVals(0, 0, n_elem_nodes-1, Nvec); // N
            B_mat2.matrixTransposeRightMultiply(1.0, B_mat1, el_mat); // N^T dN/dx
            el_mat_combined.add(q_weight[i]*jac, el_mat);
        }
        el_mat_combined.scale(2.0*mach/sound_speed);
        dof_map.addToGlobalMatrix(**el_it, el_mat_combined, global_damping_mat);
        
        
        
        // mass matrix
        el_mat_combined.zero();
        
        for (FESystemUInt i=0; i<q_pts.size(); i++)
        {
            B_mat1.zero();
            el_mat.zero();
            
            fe.getShapeFunction(*(q_pts[i]), Nvec);
            jac = fe.getJacobianValue(*(q_pts[i]));
            
            B_mat1.setRowVals(0, 0, n_elem_nodes-1, Nvec); // N
            B_mat1.matrixTransposeRightMultiply(1.0, B_mat1, el_mat); // N^T N
            el_mat_combined.add(q_weight[i]*jac, el_mat);
        }
        el_mat_combined.scale(1.0/sound_speed/sound_speed);
        dof_map.addToGlobalMatrix(**el_it, el_mat_combined, global_mass_mat);
    }
    
    
    FESystem::LinearSolvers::LapackLinearSolver<FESystemDouble> linear_solver;
    std::fstream output_file;
    FESystem::OutputProcessor::VtkOutputProcessor output;
    
    std::vector<FESystemUInt> vars(1);
    
    FESystemDouble final_t=.01, time_step=1.0e-5;
    
    // initialize the solver
    FESystem::TransientSolvers::LinearNewmarkTransientSolver<FESystemDouble> transient_solver;
    std::vector<FESystemDouble> int_constants(2); int_constants[0]=0.5; int_constants[1]=0.5;
    transient_solver.initialize(2, dof_map.getNDofs(), int_constants);
    
    transient_solver.initializeStateVector(rhs); //  initialize the vector and apply the initial condition
    // the dofs at the left edge will have a unit value specified for them
    sol.setAllVals(0.0);
    //    for (FESystemUInt i=0; i<nodes.size(); i++)
    //        if (nodes[i]->getVal(0) == 0.0)
    //            sol.setVal(nodes[i]->getDegreeOfFreedomUnit(0).global_dof_id[0], 1.0);
    transient_solver.updateVectorValuesForDerivativeOrder(0, sol, rhs);
    transient_solver.setInitialTimeData(0, time_step, rhs);
    rhs.resize(dof_map.getNDofs());
    
    transient_solver.setLinearSolver(linear_solver, true);
    global_stiffness_mat.scale(-sound_speed*sound_speed); // multiplied with sound_speed^2 since the mass matrix is assumed to be diagonal, and it is scaled by 1/sound_speed^2
    transient_solver.updateJacobianValuesForDerivativeOrder(2, 0, global_stiffness_mat, transient_solver.getCurrentJacobianMatrix());
    global_damping_mat.scale(-sound_speed*sound_speed);
    transient_solver.updateJacobianValuesForDerivativeOrder(2, 1, global_damping_mat, transient_solver.getCurrentJacobianMatrix());
    
    FESystemUInt n_skip=5, n_count=0, n_write=0;
    FESystem::TransientSolvers::TransientSolverCallBack call_back;
    while (transient_solver.getCurrentTime()<final_t)
    {
        call_back = transient_solver.incrementTimeStep();
        switch (call_back)
        {
            case FESystem::TransientSolvers::TIME_STEP_CONVERGED:
            {
                if (n_count == n_skip)
                {
                    std::stringstream oss;
                    oss << "sol_" << n_write << ".vtk";
                    output_file.open(oss.str().c_str(),std::fstream::out);
                    output.writeMesh(output_file, mesh, dof_map);
                    transient_solver.extractVectorValuesForDerivativeOrder(0, transient_solver.getCurrentStateVector(), sol); // get the current X
                    output.writeSolution(output_file, "Sol", mesh, dof_map, vars, sol);
                    output_file.close();
                    
                    n_write++;
                    n_count=0;
                }
                else
                    n_count++;
            }
                break;
                
            case FESystem::TransientSolvers::EVALUATE_X_DOT:
            case FESystem::TransientSolvers::EVALUATE_X_DOT_AND_X_DOT_JACOBIAN:
                // the Jacobian is not updated since it is constant with respect to time
            {
                transient_solver.extractVectorValuesForDerivativeOrder(0, transient_solver.getCurrentStateVector(), sol);
                global_stiffness_mat.rightVectorMultiply(sol, rhs); // -K x
                
                tmp_vec.zero();
                // create the boundary conditions on all elements
                el_it=elems.begin(), el_end=elems.end();
                for ( ; el_it != el_end; el_it++)
                {
                    dof_map.getFromGlobalVector(**el_it, sol, elem_dof_vals);
                    (*el_it)->calculateBoundaryNormal(1, point);
                    
                    // inflow on left edge
                    if ((*el_it)->getNode(0).getVal(0) == 0.0) // left edge
                    {
                        //
                        //                        for (FESystemUInt j=0; j<dim; j++) // each boundary is scaled differently with the surface normal component
                        //                        {
                        //                            el_vec.zero();
                        //
                        //                            for (FESystemUInt k=0; k<derivatives.size(); k++) derivatives[k] = 0;
                        //                            derivatives[j] = 1; // derivative of the x-coordinate
                        //
                        //                            for (FESystemUInt i=0; i<q_pts_b.size(); i++)
                        //                            {
                        //                                fe.getShapeFunctionForBoundary(3, *(q_pts_b[i]), Nvec); // bottom edge is 0 boundary
                        //                                fe.getShapeFunctionDerivativeForPhysicalCoordinatesForBoundary(3, derivatives, *(q_pts_b[i]), dNdx);
                        //                                jac = fe.getJacobianValueForBoundary(3, *(q_pts_b[i]));
                        //
                        //                                el_vec.add(mach/sound_speed * q_weight_b[i]*jac, Nvec);
                        //                            }
                        //                            el_vec.scale(point.getVal(j));
                        //
                        //                            // now add this to the global matrix: stiffness matrix
                        //                            dof_map.addToGlobalVector(**el_it, el_vec, tmp_vec);
                        //                        }
                    }
                    // non-reflecting boundary conditions at the right edge
                    else if ((*el_it)->getNode(1).getVal(0) == x_length) // right edge
                    {
                        //
                        //                        for (FESystemUInt j=0; j<dim; j++) // each boundary is scaled differently with the surface normal component
                        //                        {
                        //                            el_vec.zero();
                        //
                        //                            for (FESystemUInt k=0; k<derivatives.size(); k++) derivatives[k] = 0;
                        //                            derivatives[j] = 1; // derivative of the x-coordinate
                        //
                        //                            for (FESystemUInt i=0; i<q_pts_b.size(); i++)
                        //                            {
                        //                                fe.getShapeFunctionForBoundary(1, *(q_pts_b[i]), Nvec); // bottom edge is 0 boundary
                        //                                fe.getShapeFunctionDerivativeForPhysicalCoordinatesForBoundary(1, derivatives, *(q_pts_b[i]), dNdx);
                        //                                jac = fe.getJacobianValueForBoundary(1, *(q_pts_b[i]));
                        //
                        //                                el_vec.add(mach/sound_speed * q_weight_b[i]*jac, Nvec);
                        //                            }
                        //                            el_vec.scale(-point.getVal(j));
                        //
                        //                            // now add this to the global matrix: stiffness matrix
                        //                            dof_map.addToGlobalVector(**el_it, el_vec, tmp_vec);
                        //                        }
                    }
                    else if (((*el_it)->getNode(0).getVal(1) == 0.0) &&
                             ((*el_it)->getNode(0).getVal(0) > x_airfoil_begin) && ((*el_it)->getNode(0).getVal(0) < x_airfoil_begin+chord)) // bottom edge
                    {
                        x_val = (*el_it)->getNode(0).getVal(0)-x_airfoil_begin-chord/2.0;
                        point.setVal(0, 2.0*airfoil_t*x_val/chord/chord*4);
                        point.setVal(1, 1.0);
                        point.scaleToUnitLength();
                        
                        for (FESystemUInt j=0; j<dim; j++) // each boundary is scaled differently with the surface normal component
                        {
                            el_vec.zero();
                            
                            for (FESystemUInt k=0; k<derivatives.size(); k++) derivatives[k] = 0;
                            derivatives[j] = 1; // derivative of the x-coordinate
                            
                            for (FESystemUInt i=0; i<q_pts_b.size(); i++)
                            {
                                fe.getShapeFunctionForBoundary(0, *(q_pts_b[i]), Nvec); // bottom edge is 0 boundary
                                jac = fe.getJacobianValueForBoundary(0, *(q_pts_b[i]));
                                
                                el_vec.add(q_weight_b[i]*jac, Nvec);
                            }
                            el_vec.scale(mach/sound_speed*point.getVal(j));
                            // now add this to the global matrix: stiffness matrix
                            dof_map.addToGlobalVector(**el_it, el_vec, tmp_vec);
                        }
                    }
                    else if ((*el_it)->getNode(2).getVal(1) == y_length) // top edge
                    {
                        
                    }
                }
                
                rhs.add(1.0, tmp_vec);
                transient_solver.extractVectorValuesForDerivativeOrder(1, transient_solver.getCurrentStateVector(), sol);
                global_damping_mat.rightVectorMultiply(sol, tmp_vec); // -C x_dot
                rhs.add(1.0, tmp_vec); // -K x - C x_dot
                
                transient_solver.updateVectorValuesForDerivativeOrder(1, rhs, transient_solver.getCurrentStateVelocityVector());
                transient_solver.copyDerivativeValuesFromStateToVelocityVector(transient_solver.getCurrentStateVector(), transient_solver.getCurrentStateVelocityVector());
            }
                break;
                
            default:
                break;
        }
    }
    
    return 0;
}





int test_ode_integration(int argc, char * const argv[])
{
    
    //FESystem::Plotting::PLPlot<FESystemDouble> plot(FESystem::Plotting::REAL_AXIS, FESystem::Plotting::REAL_AXIS);
    FESystemDouble omega2=250.0, final_t=1.0/(sqrt(omega2)/2.0/3.141)*10, time_step=final_t*1.0e-4;
    
    // initialize the solver
    FESystem::TransientSolvers::LinearNewmarkTransientSolver<FESystemDouble> transient_solver;
    std::vector<FESystemDouble> int_constants(2); int_constants[0]=1.0/4.0; int_constants[1]=1.0/1.5;
    transient_solver.initialize(2, 1, int_constants);
    transient_solver.setMassMatrix(true);
    std::vector<FESystemBoolean> ode_order_include(2); ode_order_include[0] = true; ode_order_include[1]=false;
    transient_solver.setActiveJacobianTerm(ode_order_include);
    
    
    //    FESystem::Solvers::ExplicitRungeKuttaTransientSolver<FESystemDouble> transient_solver;
    //    transient_solver.initialize(2, 1, 4);
    
    FESystem::Numerics::LocalVector<FESystemDouble> vec, x_vals, y_vals;
    FESystem::Numerics::DenseMatrix<FESystemDouble> jac, ode_jac;
    vec.resize(2); jac.resize(1, 1);
    vec.setVal(0, 1.0);
    x_vals.resize(final_t/time_step+1); y_vals.resize(final_t/time_step+1);
    
    transient_solver.setInitialTimeData(0, time_step, vec);
    
    FESystem::LinearSolvers::LapackLinearSolver<FESystemDouble> linear_solver;
    transient_solver.setLinearSolver(linear_solver, true);
    transient_solver.resizeMatrixToJacobianTemplate(ode_jac);
    transient_solver.setJacobianMatrix(ode_jac);
    
    jac.setVal(0, 0, -omega2);
    transient_solver.updateJacobianValuesForDerivativeOrder(2, 0, jac, transient_solver.getCurrentJacobianMatrix());
    
    FESystem::TransientSolvers::TransientSolverCallBack call_back;
    while (transient_solver.getCurrentTime()<final_t)
    {
        call_back = transient_solver.incrementTimeStep();
        switch (call_back)
        {
            case FESystem::TransientSolvers::TIME_STEP_CONVERGED:
            {
                x_vals.setVal(transient_solver.getCurrentIterationNumber()-1, transient_solver.getCurrentTime());
                y_vals.setVal(transient_solver.getCurrentIterationNumber()-1, transient_solver.getCurrentStateVector().getVal(0));
            }
                break;
                
            case FESystem::TransientSolvers::EVALUATE_X_DOT:
            case FESystem::TransientSolvers::EVALUATE_X_DOT_AND_X_DOT_JACOBIAN:
                // the Jacobian is not updated since it is constant with respect to time
            {
                transient_solver.getCurrentStateVelocityVector().setVal(1, -omega2*transient_solver.getCurrentStateVector().getVal(0));
                transient_solver.copyDerivativeValuesFromStateToVelocityVector(transient_solver.getCurrentStateVector(), transient_solver.getCurrentStateVelocityVector());
            }
                break;
                
            default:
                break;
        }
    }
    
    //plot.plotData2D(x_vals, y_vals);
    
    
    return 0;
}




void staticAnalysis(FESystemUInt dim, const FESystem::Mesh::MeshBase& mesh, const FESystem::Base::DegreeOfFreedomMap& dof_map, const FESystem::Numerics::SparsityPattern& nonbc_sparsity_pattern,
                    const std::map<FESystemUInt, FESystemUInt>& old_to_new_id_map, const std::vector<FESystemUInt>& nonbc_dofs, const FESystem::Numerics::MatrixBase<FESystemDouble>& stiff_mat,
                    const FESystem::Numerics::VectorBase<FESystemDouble>& rhs, FESystem::Numerics::VectorBase<FESystemDouble>& sol)
{
    const std::vector<FESystem::Mesh::Node*>& nodes = mesh.getNodes();
    
    FESystem::Numerics::SparseMatrix<FESystemDouble> reduced_stiff_mat;
    FESystem::Numerics::LocalVector<FESystemDouble> reduced_load_vec, reduced_sol_vec;
    
    reduced_stiff_mat.resize(nonbc_sparsity_pattern);
    reduced_load_vec.resize(nonbc_sparsity_pattern.getNDOFs()); reduced_sol_vec.resize(nonbc_sparsity_pattern.getNDOFs());
    rhs.getSubVectorValsFromIndices(nonbc_dofs, reduced_load_vec);
    
    stiff_mat.getSubMatrixValsFromRowAndColumnIndices(nonbc_dofs, nonbc_dofs, old_to_new_id_map, reduced_stiff_mat);
    
    FESystem::LinearSolvers::LUFactorizationLinearSolver<FESystemDouble> linear_solver;
    linear_solver.setSystemMatrix(reduced_stiff_mat);
    linear_solver.solve(reduced_load_vec, reduced_sol_vec);
    sol.setSubVectorValsFromIndices(nonbc_dofs, reduced_sol_vec);
    // write the solution for each node
    for (FESystemUInt i=0; i<nodes.size(); i++)
    {
        std::cout << "Node: " << std::setw(8) << i;
        // write location
        for (FESystemUInt j=0; j<dim; j++)
            std::cout << std::setw(15) << nodes[i]->getVal(j);
        // write the forces
        for (FESystemUInt j=0; j<3; j++)
            std::cout << std::setw(15) << rhs.getVal(nodes[i]->getDegreeOfFreedomUnit(j).global_dof_id[0]);
        for (FESystemUInt j=0; j<3; j++)
            std::cout << std::setw(15) << sol.getVal(nodes[i]->getDegreeOfFreedomUnit(j).global_dof_id[0]);
        std::cout << std::endl;
    }
    
    // write a gmsh format file
    std::vector<FESystemUInt> vars(3); vars[0]=0; vars[1]=1; vars[2]=2; // write all solutions
    FESystem::OutputProcessor::VtkOutputProcessor output;
    FESystem::OutputProcessor::GmshOutputProcessor gmsh_output;
    
    std::fstream output_file;
    output_file.open("gmsh_output.gmsh", std::fstream::out);
    gmsh_output.writeMesh(output_file, mesh, dof_map);
    output_file.close();
    
    
    output_file.open("vtk_output.vtk", std::fstream::out);
    output.writeMesh(output_file, mesh, dof_map);
    output.writeSolution(output_file, "Solution", mesh, dof_map, vars, sol);
    output_file.close();
}
 


void nonlinearSolution(FESystemUInt dim, FESystem::Mesh::ElementType elem_type, FESystemUInt n_elem_nodes, const
                       FESystem::Mesh::MeshBase& mesh, const FESystem::Base::DegreeOfFreedomMap& dof_map,
                       const FESystem::Numerics::SparsityPattern& nonbc_sparsity_pattern,
                       const std::map<FESystemUInt, FESystemUInt>& old_to_new_id_map, const std::vector<FESystemUInt>& nonbc_dofs,
                       FESystem::Numerics::MatrixBase<FESystemDouble>& stiff_mat,
                       FESystem::Numerics::VectorBase<FESystemDouble>& rhs, FESystem::Numerics::VectorBase<FESystemDouble>& sol,
                       void (*calculateBeamMatrices)(FESystemBoolean if_nonlinear, FESystem::Mesh::ElementType elem_type, FESystemUInt n_elem_nodes,
                                                     const FESystem::Base::DegreeOfFreedomMap& dof_map,
                                                     const FESystem::Mesh::MeshBase& mesh,
                                                     FESystem::Numerics::VectorBase<FESystemDouble>& global_sol,
                                                     FESystem::Numerics::VectorBase<FESystemDouble>& internal_force,
                                                     FESystem::Numerics::VectorBase<FESystemDouble>& external_force,
                                                     FESystem::Numerics::MatrixBase<FESystemDouble>& global_stiffness_mat,
                                                     FESystem::Numerics::VectorBase<FESystemDouble>& global_mass_vec))
{
    const std::vector<FESystem::Mesh::Node*>& nodes = mesh.getNodes();
    
    FESystem::Numerics::SparseMatrix<FESystemDouble> reduced_stiff_mat;
    FESystem::Numerics::LocalVector<FESystemDouble> reduced_load_vec, reduced_sol_vec, internal_force, dummy;

    internal_force.resize(dof_map.getNDofs());
    reduced_stiff_mat.resize(nonbc_sparsity_pattern);
    reduced_load_vec.resize(nonbc_sparsity_pattern.getNDOFs()); reduced_sol_vec.resize(nonbc_sparsity_pattern.getNDOFs());
    stiff_mat.getSubMatrixValsFromRowAndColumnIndices(nonbc_dofs, nonbc_dofs, old_to_new_id_map, reduced_stiff_mat);
    
    FESystem::LinearSolvers::LUFactorizationLinearSolver<FESystemDouble> linear_solver;
    FESystem::NonlinearSolvers::NewtonIterationNonlinearSolver<FESystemDouble> nonlinear_solver;
    
    nonlinear_solver.initialize(reduced_stiff_mat, linear_solver);
    
    FESystem::NonlinearSolvers::NonlinearSolverCallBack call_back = nonlinear_solver.getCurrentCallBack();
    
    while (call_back != FESystem::NonlinearSolvers::SOLUTION_CONVERGED)
    {
        switch (call_back)
        {
            case FESystem::NonlinearSolvers::WAITING_TO_START:
                // nothing to be done
                break;

            case FESystem::NonlinearSolvers::SET_INITIAL_GUESS:
                nonlinear_solver.getCurrentSolution().zero();
                break;

            case FESystem::NonlinearSolvers::EVALUATE_RESIDUAL:
            {
                // zero the solution vectors
                rhs.zero(); sol.zero();
                // get the latest solution vector
                sol.addVal(nonbc_dofs, nonlinear_solver.getCurrentSolution());
                calculateBeamMatrices(true, elem_type, n_elem_nodes, dof_map, mesh, sol, internal_force, rhs, stiff_mat, dummy);
                internal_force.getSubVectorValsFromIndices(nonbc_dofs, nonlinear_solver.getResidualVector());
                rhs.getSubVectorValsFromIndices(nonbc_dofs, reduced_load_vec);
                nonlinear_solver.getResidualVector().add(-1.0, reduced_load_vec);
            }
                break;

            case FESystem::NonlinearSolvers::EVALUATE_JACOBIAN:
            {
                // zero the solution vectors
                rhs.zero(); sol.zero();
                // get the latest solution vector
                sol.addVal(nonbc_dofs, nonlinear_solver.getCurrentSolution());
                calculateBeamMatrices(true, elem_type, n_elem_nodes, dof_map, mesh, sol, internal_force, rhs, stiff_mat, dummy);
                stiff_mat.getSubMatrixValsFromRowAndColumnIndices(nonbc_dofs, nonbc_dofs, old_to_new_id_map, reduced_stiff_mat);
            }
                break;
                
            case FESystem::NonlinearSolvers::EVALUATE_RESIDUAL_AND_JACOBIAN:
            {
                // zero the solution vectors
                rhs.zero(); sol.zero();
                // get the latest solution vector
                sol.addVal(nonbc_dofs, nonlinear_solver.getCurrentSolution());
                calculateBeamMatrices(true, elem_type, n_elem_nodes, dof_map, mesh, sol, internal_force, rhs, stiff_mat, dummy);
                internal_force.getSubVectorValsFromIndices(nonbc_dofs, nonlinear_solver.getResidualVector());
                rhs.getSubVectorValsFromIndices(nonbc_dofs, reduced_load_vec);
                nonlinear_solver.getResidualVector().add(-1.0, reduced_load_vec); // residual vector
                stiff_mat.getSubMatrixValsFromRowAndColumnIndices(nonbc_dofs, nonbc_dofs, old_to_new_id_map, reduced_stiff_mat); // tangent stiffness matrix
            }
                break;

            default:
                break;
        }
        call_back = nonlinear_solver.incrementSolution();
    }
    
    sol.setSubVectorValsFromIndices(nonbc_dofs, reduced_sol_vec);
    // write the solution for each node
    for (FESystemUInt i=0; i<nodes.size(); i++)
    {
        std::cout << "Node: " << std::setw(8) << i;
        // write location
        for (FESystemUInt j=0; j<dim; j++)
            std::cout << std::setw(15) << nodes[i]->getVal(j);
        // write the forces
        for (FESystemUInt j=0; j<3; j++)
            std::cout << std::setw(15) << rhs.getVal(nodes[i]->getDegreeOfFreedomUnit(j).global_dof_id[0]);
        for (FESystemUInt j=0; j<3; j++)
            std::cout << std::setw(15) << sol.getVal(nodes[i]->getDegreeOfFreedomUnit(j).global_dof_id[0]);
        std::cout << std::endl;
    }
    
    // write a gmsh format file
    std::vector<FESystemUInt> vars(3); vars[0]=0; vars[1]=1; vars[2]=2; // write all solutions
    FESystem::OutputProcessor::VtkOutputProcessor output;
    FESystem::OutputProcessor::GmshOutputProcessor gmsh_output;
    
    std::fstream output_file;
    output_file.open("gmsh_output.gmsh", std::fstream::out);
    gmsh_output.writeMesh(output_file, mesh, dof_map);
    output_file.close();
    
    
    output_file.open("vtk_output.vtk", std::fstream::out);
    output.writeMesh(output_file, mesh, dof_map);
    output.writeSolution(output_file, "Solution", mesh, dof_map, vars, sol);
    output_file.close();
}



void modalAnalysis(FESystemUInt dim, const FESystem::Mesh::MeshBase& mesh, const FESystem::Base::DegreeOfFreedomMap& dof_map, const FESystem::Numerics::SparsityPattern& nonbc_sparsity_pattern,
                   const std::map<FESystemUInt, FESystemUInt>& old_to_new_id_map, const std::vector<FESystemUInt>& nonbc_dofs, const FESystem::Numerics::MatrixBase<FESystemDouble>& stiff_mat,
                   const FESystem::Numerics::VectorBase<FESystemDouble>& mass_vec, const FESystemUInt n_modes, std::vector<FESystemUInt>& sorted_ids, FESystem::Numerics::VectorBase<FESystemDouble>& eig_vals,
                   FESystem::Numerics::MatrixBase<FESystemDouble>& eig_vec)
{
    FESystem::LinearSolvers::LUFactorizationLinearSolver<FESystemDouble> linear_solver;
    FESystem::Numerics::SparseMatrix<FESystemDouble> reduced_stiff_mat, mass_mat;
    FESystem::Numerics::LocalVector<FESystemDouble> reduced_load_vec, reduced_sol_vec, sol, reduced_mass_vec;

    reduced_stiff_mat.resize(nonbc_sparsity_pattern); mass_mat.resize(nonbc_sparsity_pattern);
    reduced_sol_vec.resize(nonbc_sparsity_pattern.getNDOFs()); reduced_mass_vec.resize(nonbc_sparsity_pattern.getNDOFs());
    sol.resize(dof_map.getNDofs());

    stiff_mat.getSubMatrixValsFromRowAndColumnIndices(nonbc_dofs, nonbc_dofs, old_to_new_id_map, reduced_stiff_mat);
    
    mass_vec.getSubVectorValsFromIndices(nonbc_dofs, reduced_mass_vec);
    mass_mat.setDiagonal(reduced_mass_vec);
    
    // modal eigensolution
    FESystem::EigenSolvers::ArpackLinearEigenSolver<FESystemDouble> eigen_solver;
    //FESystem::EigenSolvers::LapackLinearEigenSolver<FESystemDouble> eigen_solver;
    eigen_solver.setEigenProblemType(FESystem::EigenSolvers::GENERALIZED_HERMITIAN);
    eigen_solver.setMatrix(&reduced_stiff_mat, &mass_mat);
    eigen_solver.setEigenShiftType(FESystem::EigenSolvers::SHIFT_AND_INVERT); eigen_solver.setEigenShiftValue(0.0);
    linear_solver.clear(); eigen_solver.setLinearSolver(linear_solver);
    eigen_solver.setEigenSpectrumType(FESystem::EigenSolvers::LARGEST_MAGNITUDE);
    eigen_solver.init(n_modes, true);
    eigen_solver.solve();
    eig_vals.copyVector(eigen_solver.getEigenValues());
    eig_vec.copyMatrix(eigen_solver.getEigenVectorMatrix());
    eigen_solver.prepareSortingVector(FESystem::EigenSolvers::VALUE, sorted_ids);
    
    FESystemUInt id = 0;
    FESystem::OutputProcessor::VtkOutputProcessor output;
    std::fstream output_file;
    std::vector<FESystemUInt> vars(3); vars[0]=0; vars[1]=1; vars[2]=2; // write all solutions

    std::cout << "EigenValues:" << std::endl;
    output_file.open("vtk_modes.vtk", std::fstream::out);
    output.writeMesh(output_file, mesh, dof_map);
    std::pair<FESystemUInt, FESystemUInt> s_global = eig_vec.getSize();
    for (FESystemUInt i=0; i<n_modes; i++)
    {
        id = sorted_ids[i];
        std::cout << id << "  " << std::setw(15) << eig_vals.getVal(id) << std::endl;
        std::stringstream oss;
        oss << "Sol_" << i;
        eig_vec.getColumnVals(id, 0, s_global.first-1, reduced_sol_vec);
        sol.setSubVectorValsFromIndices(nonbc_dofs, reduced_sol_vec);
        output.writeSolution(output_file, oss.str(), mesh, dof_map, vars, sol);
    }
    output_file.close();
}





void transientAnalysis(FESystemUInt dim, const FESystem::Mesh::MeshBase& mesh, const FESystem::Base::DegreeOfFreedomMap& dof_map, const FESystem::Numerics::SparsityPattern& nonbc_sparsity_pattern,
                       const std::map<FESystemUInt, FESystemUInt>& old_to_new_id_map, const std::vector<FESystemUInt>& nonbc_dofs, const FESystem::Numerics::MatrixBase<FESystemDouble>& stiff_mat,
                       const FESystem::Numerics::VectorBase<FESystemDouble>& mass_vec, const std::vector<FESystemUInt>& sorted_ids, const FESystem::Numerics::VectorBase<FESystemDouble>& eig_vals,
                       const FESystem::Numerics::MatrixBase<FESystemDouble>& eig_vec)
{
    FESystem::LinearSolvers::LUFactorizationLinearSolver<FESystemDouble> linear_solver;
    FESystem::Numerics::SparseMatrix<FESystemDouble> reduced_stiff_mat, mass_mat;
    FESystem::Numerics::LocalVector<FESystemDouble> reduced_load_vec, reduced_sol_vec, rhs, sol, reduced_mass_vec;
    
    reduced_stiff_mat.resize(nonbc_sparsity_pattern); mass_mat.resize(nonbc_sparsity_pattern);
    reduced_sol_vec.resize(nonbc_sparsity_pattern.getNDOFs()); reduced_load_vec.resize(nonbc_sparsity_pattern.getNDOFs());
    reduced_mass_vec.resize(nonbc_sparsity_pattern.getNDOFs());
    sol.resize(dof_map.getNDofs());
    
    stiff_mat.getSubMatrixValsFromRowAndColumnIndices(nonbc_dofs, nonbc_dofs, old_to_new_id_map, reduced_stiff_mat);
    mass_vec.getSubVectorValsFromIndices(nonbc_dofs, reduced_mass_vec);

    FESystemUInt id= 0;
    
    // transient solution
    // initial condition is the first mode
    // scale the stiffness matrix with the mass matrix
    for (FESystemUInt i=0; i<reduced_mass_vec.getSize(); i++)
        reduced_stiff_mat.scaleRow(i, -1.0/reduced_mass_vec.getVal(i));
    
    id = sorted_ids[0];
    eig_vec.getColumnVals(id, 0, nonbc_dofs.size()-1, reduced_sol_vec);
    FESystemDouble final_t=1.0/(sqrt(eig_vals.getVal(id))/2.0/3.141)*4, time_step=final_t*1.0e-3;
    
    // initialize the solver
    FESystem::TransientSolvers::LinearNewmarkTransientSolver<FESystemDouble> transient_solver;
    FESystem::Numerics::SparsityPattern ode_sparsity;
    FESystem::Numerics::SparseMatrix<FESystemDouble> ode_jac;
    std::vector<FESystemBoolean> ode_order_include(2); ode_order_include[0] = true; ode_order_include[1]=false;
    std::vector<FESystemDouble> int_constants(2); int_constants[0]=1.0/4.0; int_constants[1]=1.0/2.0;
    transient_solver.initialize(2, nonbc_dofs.size(), int_constants);
    transient_solver.setActiveJacobianTerm(ode_order_include);
    transient_solver.setMassMatrix(true);
    
    //    FESystem::Solvers::ExplicitRungeKuttaTransientSolver<FESystemDouble> transient_solver;
    //    transient_solver.initialize(2, nonbc_dofs.size(), 4);
    
    transient_solver.initializeStateVector(rhs); //  initialize the vector and apply the initial condition
    transient_solver.updateVectorValuesForDerivativeOrder(0, reduced_sol_vec, rhs);
    transient_solver.setInitialTimeData(0, time_step, rhs);
    
    transient_solver.initializeMatrixSparsityPatterForSystem(nonbc_sparsity_pattern, ode_sparsity);
    ode_jac.resize(ode_sparsity);
    transient_solver.setJacobianMatrix(ode_jac);
    transient_solver.setLinearSolver(linear_solver, true);
    transient_solver.updateJacobianValuesForDerivativeOrder(2, 0, reduced_stiff_mat, transient_solver.getCurrentJacobianMatrix());
    
    FESystem::OutputProcessor::VtkOutputProcessor output;
    std::fstream output_file;
    std::vector<FESystemUInt> vars(3); vars[0]=0; vars[1]=1; vars[2]=2; // write all solutions
    
    FESystemUInt n_skip=20, n_count=0, n_write=0;
    FESystem::TransientSolvers::TransientSolverCallBack call_back;
    while (transient_solver.getCurrentTime()<final_t)
    {
        call_back = transient_solver.incrementTimeStep();
        switch (call_back)
        {
            case FESystem::TransientSolvers::TIME_STEP_CONVERGED:
            {
                if (n_count == n_skip)
                {
                    std::stringstream oss;
                    oss << "sol_" << n_write << ".vtk";
                    output_file.open(oss.str().c_str(),std::fstream::out);
                    output.writeMesh(output_file, mesh, dof_map);
                    transient_solver.extractVectorValuesForDerivativeOrder(0, transient_solver.getCurrentStateVector(), reduced_sol_vec); // get the current X
                    sol.setSubVectorValsFromIndices(nonbc_dofs, reduced_sol_vec);
                    output.writeSolution(output_file, "Sol", mesh, dof_map, vars, sol);
                    output_file.close();
                    
                    n_write++;
                    n_count=0;
                }
                else
                    n_count++;
            }
                break;
                
            case FESystem::TransientSolvers::EVALUATE_X_DOT:
            case FESystem::TransientSolvers::EVALUATE_X_DOT_AND_X_DOT_JACOBIAN:
                // the Jacobian is not updated since it is constant with respect to time
            {
                transient_solver.extractVectorValuesForDerivativeOrder(0, transient_solver.getCurrentStateVector(), reduced_sol_vec); // get the current X
                reduced_stiff_mat.rightVectorMultiply(reduced_sol_vec, reduced_load_vec);
                
                transient_solver.updateVectorValuesForDerivativeOrder(1, reduced_load_vec, transient_solver.getCurrentStateVelocityVector()); // set the acceleration
                transient_solver.copyDerivativeValuesFromStateToVelocityVector(transient_solver.getCurrentStateVector(), transient_solver.getCurrentStateVelocityVector());
            }
                break;
                
            default:
                FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
                break;
        }
    }
}





int plate_analysis_driver(int argc, char * const argv[])
{
    FESystem::Mesh::ElementType elem_type;
    
    // create the mesh object
    FESystem::Mesh::MeshBase mesh;
    
    // create a nx x ny grid of nodes
    FESystemUInt nx, ny, dim, n_elem_nodes, n_elem_dofs, n_modes;
    FESystemDouble x_length, y_length, p_val = 1.0e6;
    FESystem::Geometry::Point origin(3);
    
    FESystemBoolean if_mindlin = false;
    
    nx=7; ny=7;
    x_length = 2; y_length = 2;
    dim = 2; n_modes = 20;
    elem_type = FESystem::Mesh::TRI3;
    create_plane_mesh(elem_type, mesh, origin, nx, ny, x_length, y_length, n_elem_nodes, CROSS);
    
    n_elem_dofs = 6*n_elem_nodes;
    
    // set the location of individual nodes
    const std::vector<FESystem::Mesh::Node*>& nodes = mesh.getNodes();
    
    // now add the degrees of freedom
    FESystem::Base::DegreeOfFreedomMap dof_map(mesh);
    std::string name;
    name = "u"; dof_map.addVariable(name, 0);
    name = "v"; dof_map.addVariable(name, 0);
    name = "w"; dof_map.addVariable(name, 0);
    name = "tx"; dof_map.addVariable(name, 0);
    name = "ty"; dof_map.addVariable(name, 0);
    name = "tz"; dof_map.addVariable(name, 0);
    dof_map.reinit(); // distribute the dofs
    
    // create the finite element and initialize the shape functions
    FESystem::Numerics::SparseMatrix<FESystemDouble> global_stiffness_mat;
    FESystem::Numerics::LocalVector<FESystemDouble> rhs, sol, global_mass_vec, elem_vec, plate_elem_vec;
    FESystem::Numerics::DenseMatrix<FESystemDouble> elem_mat, plate_elem_mat;
    std::vector<FESystemUInt> elem_dof_indices;
    
    global_mass_vec.resize(dof_map.getNDofs());
    rhs.resize(dof_map.getNDofs()); sol.resize(dof_map.getNDofs());
    global_stiffness_mat.resize(dof_map.getSparsityPattern());
    elem_mat.resize(n_elem_dofs, n_elem_dofs);
    elem_vec.resize(n_elem_dofs);
    
    plate_elem_mat.resize(n_elem_nodes*3, n_elem_nodes*3);  plate_elem_vec.resize(n_elem_nodes*3);
    
    FESystemDouble E=72.0e9, nu=0.33, thick=0.025, rho=2700.0;
    
    
    
    // prepare the quadrature rule and FE for the
    FESystem::Quadrature::TrapezoidQuadrature q_rule_shear, q_rule_bending;
    FESystem::FiniteElement::FELagrange fe, fe_tri6;
    FESystem::Structures::ReissnerMindlinPlate mindlin_plate;
    FESystem::Structures::DKTPlate dkt_plate;
    
    if (if_mindlin)
    {
        switch (elem_type)
        {
            case FESystem::Mesh::QUAD4:
            case FESystem::Mesh::TRI3:
                q_rule_bending.init(2, 3);  // bending quadrature is higher than the shear quadrature for reduced integrations
                q_rule_shear.init(2, 0);
                break;
                
            case FESystem::Mesh::QUAD9:
            case FESystem::Mesh::TRI6:
                q_rule_bending.init(2, 5);
                q_rule_shear.init(2, 5);
                break;
                
            default:
                FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
        }
    }
    else
        q_rule_bending.init(2, 9);
    

    std::vector<FESystem::Mesh::ElemBase*>& elems = mesh.getElements();
    for (FESystemUInt i=0; i<elems.size(); i++)
    {
        fe.clear();
        fe.reinit(*(elems[i]));
        elem_mat.zero();
        elem_vec.zero();
        
        if (if_mindlin)
        {
            mindlin_plate.clear();
            mindlin_plate.initialize(*(elems[i]), fe, q_rule_bending, q_rule_shear, E, nu, rho, thick);
            
            // stiffness matrix
            mindlin_plate.calculateStiffnessMatrix(plate_elem_mat);
            mindlin_plate.transformMatrixToGlobalSystem(plate_elem_mat, elem_mat);
            
            // mass
            mindlin_plate.calculateDiagonalMassMatrix(plate_elem_vec);
            mindlin_plate.getActiveElementMatrixIndices(elem_dof_indices);
        }
        else
        {
            dkt_plate.clear();
            dkt_plate.initialize(*(elems[i]), fe, fe_tri6, q_rule_bending, q_rule_bending, E, nu, rho, thick);
            
            // stiffness matrix
            dkt_plate.calculateStiffnessMatrix(plate_elem_mat);
            dkt_plate.transformMatrixToGlobalSystem(plate_elem_mat, elem_mat);
            
            // mass
            dkt_plate.calculateDiagonalMassMatrix(plate_elem_vec);
            dkt_plate.getActiveElementMatrixIndices(elem_dof_indices);
        }
        elem_vec.zero();
        elem_vec.addVal(elem_dof_indices, plate_elem_vec);
        
        dof_map.addToGlobalMatrix(*(elems[i]), elem_mat, global_stiffness_mat);
        dof_map.addToGlobalVector(*(elems[i]), elem_vec, global_mass_vec);
        
    }
    
    // apply boundary condition and place a load on the last dof
    std::set<FESystemUInt> bc_dofs;
    for (FESystemUInt i=0; i<nodes.size(); i++)
    {
        bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(0).global_dof_id[0]); // u-disp
        bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(1).global_dof_id[0]); // v-disp
        bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(5).global_dof_id[0]); // theta-z
        if ((nodes[i]->getVal(0) == 0.0) || (nodes[i]->getVal(1) == 0.0) || (nodes[i]->getVal(0) == x_length) || (nodes[i]->getVal(1) == y_length)) // boundary nodes
            bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(2).global_dof_id[0]); // w-displacement on the bottom edge
    }
    
    // now create the vector of ids that do not have bcs
    std::vector<FESystemUInt> nonbc_dofs;
    for (FESystemUInt i=0; i<dof_map.getNDofs(); i++)
        if (!bc_dofs.count(i))
            nonbc_dofs.push_back(i);
    
    // prepare a map of the old to new ID
    std::vector<FESystemUInt>::const_iterator dof_it=nonbc_dofs.begin(), dof_end=nonbc_dofs.end();
    std::map<FESystemUInt, FESystemUInt> old_to_new_id_map;
    FESystemUInt n=0;
    for ( ; dof_it!=dof_end; dof_it++)
        old_to_new_id_map.insert(std::map<FESystemUInt, FESystemUInt>::value_type(*dof_it, n++));
    
    FESystem::Numerics::SparsityPattern nonbc_sparsity_pattern;
    dof_map.getSparsityPattern().initSparsityPatternForNonConstrainedDOFs(nonbc_dofs, old_to_new_id_map, nonbc_sparsity_pattern);
    
    std::vector<FESystemUInt> sorted_ids;
    FESystem::Numerics::LocalVector<FESystemDouble> eig_vals;
    FESystem::Numerics::DenseMatrix<FESystemDouble> eig_vecs;
    
    staticAnalysis(dim, mesh, dof_map, nonbc_sparsity_pattern, old_to_new_id_map, nonbc_dofs, global_stiffness_mat, rhs, sol);
    
    modalAnalysis(dim, mesh, dof_map, nonbc_sparsity_pattern, old_to_new_id_map, nonbc_dofs, global_stiffness_mat, global_mass_vec, n_modes, sorted_ids, eig_vals, eig_vecs);
    
    transientAnalysis(dim, mesh, dof_map, nonbc_sparsity_pattern, old_to_new_id_map, nonbc_dofs, global_stiffness_mat, global_mass_vec, sorted_ids, eig_vals, eig_vecs);
 
    return 0;
}




void calculateBeamMatrices(FESystemBoolean if_nonlinear, FESystem::Mesh::ElementType elem_type, FESystemUInt n_elem_nodes, const FESystem::Base::DegreeOfFreedomMap& dof_map,
                           const FESystem::Mesh::MeshBase& mesh, FESystem::Numerics::VectorBase<FESystemDouble>& global_sol,
                           FESystem::Numerics::VectorBase<FESystemDouble>& internal_force,
                           FESystem::Numerics::VectorBase<FESystemDouble>& external_force,
                           FESystem::Numerics::MatrixBase<FESystemDouble>& global_stiffness_mat,
                           FESystem::Numerics::VectorBase<FESystemDouble>& global_mass_vec)
{
    
    FESystemDouble E=72.0e9, nu=0.33, rho=2700.0, p_val = 1.0e6, I_tr = 6.667e-9, I_ch = 1.6667e-9, area = 2.0e-4;
    FESystemBoolean if_timoshenko_beam = false;
    FESystemUInt n_beam_dofs, n_elem_dofs;
    
    n_beam_dofs = 4*n_elem_nodes;
    n_elem_dofs = 6*n_elem_nodes;

    
    FESystem::Numerics::DenseMatrix<FESystemDouble> elem_mat, beam_elem_mat;
    FESystem::Numerics::LocalVector<FESystemDouble> elem_vec, beam_elem_vec, elem_sol, elem_local_sol;
    std::vector<FESystemUInt> elem_dof_indices;
    elem_mat.resize(n_elem_dofs, n_elem_dofs);
    elem_vec.resize(n_elem_dofs); elem_sol.resize(n_elem_dofs); elem_local_sol.resize(4*n_elem_nodes);
    
    beam_elem_mat.resize(n_beam_dofs, n_beam_dofs); beam_elem_vec.resize(n_beam_dofs);

    // prepare the quadrature rule and FE for the
    FESystem::Quadrature::TrapezoidQuadrature q_rule_shear, q_rule_bending;
    FESystem::FiniteElement::FELagrange fe;
    FESystem::Structures::TimoshenkoBeam timoshenko_beam;
    FESystem::Structures::EulerBernoulliBeam euler_beam;
    FESystem::Structures::VonKarmanStrain1D vk_beam;
    FESystem::Structures::ExtensionBar bar;
    
    if (if_timoshenko_beam)
    {
        switch (elem_type)
        {
            case FESystem::Mesh::EDGE2:
                q_rule_bending.init(1, 3);
                q_rule_shear.init(1, 0);
                break;
                
            case FESystem::Mesh::EDGE3:
            {
                q_rule_bending.init(1, 5);
                q_rule_shear.init(1, 5);
            }
                break;
                
            default:
                FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
        }
    }
    else
        q_rule_bending.init(1, 9);
    
    
    const std::vector<FESystem::Mesh::ElemBase*>& elems = mesh.getElements();
    
    for (FESystemUInt i=0; i<elems.size(); i++)
    {
        fe.clear();
        fe.reinit(*(elems[i]));
        elem_mat.zero();
        elem_vec.zero();
        
        if (if_nonlinear)
        {
            dof_map.getFromGlobalVector(*(elems[i]), global_sol, elem_sol);
            if (if_timoshenko_beam)
            {
                timoshenko_beam.clear();
                timoshenko_beam.initialize(*(elems[i]), fe, q_rule_bending, q_rule_shear, E, nu, rho, I_tr, I_ch, area);
                vk_beam.clear();
                vk_beam.initialize(*(elems[i]), fe, q_rule_bending, 0.0, timoshenko_beam);
            }
            else
            {
                euler_beam.clear();
                euler_beam.initialize(*(elems[i]), fe, q_rule_bending, E, nu, rho, I_tr, I_ch, area);
                vk_beam.clear();
                vk_beam.initialize(*(elems[i]), fe, q_rule_bending, 0.0, euler_beam);
            }
            vk_beam.getActiveElementMatrixIndices(elem_dof_indices);
            elem_sol.getSubVectorValsFromIndices(elem_dof_indices, elem_local_sol);
            vk_beam.calculateInternalForceVector(elem_local_sol, beam_elem_vec);
            vk_beam.transformVectorToGlobalSystem(beam_elem_vec, elem_vec);
            vk_beam.calculateTangentStiffnessMatrix(elem_local_sol, beam_elem_mat);
            vk_beam.transformMatrixToGlobalSystem(beam_elem_mat, elem_mat);
        }
        else
        {
            if (if_timoshenko_beam)
            {
                timoshenko_beam.clear();
                timoshenko_beam.initialize(*(elems[i]), fe, q_rule_bending, q_rule_shear, E, nu, rho, I_tr, I_ch, area);
                
                // stiffness matrix
                timoshenko_beam.calculateStiffnessMatrix(beam_elem_mat);
                timoshenko_beam.transformMatrixToGlobalSystem(beam_elem_mat, elem_mat);
                
                // mass
                timoshenko_beam.calculateDiagonalMassMatrix(beam_elem_vec);
                timoshenko_beam.getActiveElementMatrixIndices(elem_dof_indices);
            }
            else
            {
                euler_beam.clear();
                euler_beam.initialize(*(elems[i]), fe, q_rule_bending, E, nu, rho, I_tr, I_ch, area);
                
                // stiffness matrix
                euler_beam.calculateStiffnessMatrix(beam_elem_mat);
                euler_beam.transformMatrixToGlobalSystem(beam_elem_mat, elem_mat);
                
                // mass
                euler_beam.calculateDiagonalMassMatrix(beam_elem_vec);
                euler_beam.getActiveElementMatrixIndices(elem_dof_indices);
            }
        }
        elem_vec.zero();
        elem_vec.addVal(elem_dof_indices, beam_elem_vec);
        
        dof_map.addToGlobalMatrix(*(elems[i]), elem_mat, global_stiffness_mat);
        dof_map.addToGlobalVector(*(elems[i]), elem_vec, global_mass_vec);
        
    }
}




int beam_analysis_driver(int argc, char * const argv[])
{
    FESystem::Mesh::ElementType elem_type;
    
    // create the mesh object
    FESystem::Mesh::MeshBase mesh;
    
    // create a nx x ny grid of nodes
    FESystemUInt nx, dim, n_elem_nodes, n_modes;
    FESystemDouble x_length;
    FESystem::Geometry::Point origin(3);
        
    nx=15; x_length = 2; dim = 1; n_modes = 8;
    elem_type = FESystem::Mesh::EDGE2;
    createLineMesh(elem_type, mesh, origin, nx, x_length, n_elem_nodes);
        
    // set the location of individual nodes
    const std::vector<FESystem::Mesh::Node*>& nodes = mesh.getNodes();
    
    // now add the degrees of freedom
    FESystem::Base::DegreeOfFreedomMap dof_map(mesh);
    std::string name;
    name = "u"; dof_map.addVariable(name, 0);
    name = "v"; dof_map.addVariable(name, 0);
    name = "w"; dof_map.addVariable(name, 0);
    name = "tx"; dof_map.addVariable(name, 0);
    name = "ty"; dof_map.addVariable(name, 0);
    name = "tz"; dof_map.addVariable(name, 0);
    dof_map.reinit(); // distribute the dofs
    
    // create the finite element and initialize the shape functions
    FESystem::Numerics::SparseMatrix<FESystemDouble> global_stiffness_mat;
    FESystem::Numerics::LocalVector<FESystemDouble> rhs, sol, global_mass_vec;
    std::vector<FESystemUInt> elem_dof_indices;
    
    global_mass_vec.resize(dof_map.getNDofs());
    rhs.resize(dof_map.getNDofs()); sol.resize(dof_map.getNDofs());
    global_stiffness_mat.resize(dof_map.getSparsityPattern());
        
    
    // apply boundary condition and place a load on the last dof
    std::set<FESystemUInt> bc_dofs;
    for (FESystemUInt i=0; i<nodes.size(); i++)
    {
        bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(0).global_dof_id[0]); // u-disp
        //bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(1).global_dof_id[0]); // v-disp
        bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(3).global_dof_id[0]); // theta-x
        //bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(5).global_dof_id[0]); // theta-z
        if ((nodes[i]->getVal(0) == 0.0) || (nodes[i]->getVal(0) == x_length)) // boundary nodes
        {
            bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(1).global_dof_id[0]); // v-displacement on the bottom edge
            bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(2).global_dof_id[0]); // w-displacement on the bottom edge
        }
    }
    
    // now create the vector of ids that do not have bcs
    std::vector<FESystemUInt> nonbc_dofs;
    for (FESystemUInt i=0; i<dof_map.getNDofs(); i++)
        if (!bc_dofs.count(i))
            nonbc_dofs.push_back(i);
    
    // prepare a map of the old to new ID
    std::vector<FESystemUInt>::const_iterator dof_it=nonbc_dofs.begin(), dof_end=nonbc_dofs.end();
    std::map<FESystemUInt, FESystemUInt> old_to_new_id_map;
    FESystemUInt n=0;
    for ( ; dof_it!=dof_end; dof_it++)
        old_to_new_id_map.insert(std::map<FESystemUInt, FESystemUInt>::value_type(*dof_it, n++));
    
    FESystem::Numerics::SparsityPattern nonbc_sparsity_pattern;
    dof_map.getSparsityPattern().initSparsityPatternForNonConstrainedDOFs(nonbc_dofs, old_to_new_id_map, nonbc_sparsity_pattern);

    
//    nonlinearSolution(dim, elem_type, n_elem_nodes, mesh, dof_map, nonbc_sparsity_pattern, old_to_new_id_map, nonbc_dofs, global_stiffness_mat, rhs, sol, calculateBeamMatrices);
//    exit(1);
    
    // calculate the stiffness quantities for the linearized analyses
    calculateBeamMatrices(false, elem_type, n_elem_nodes, dof_map, mesh, sol, rhs, rhs, global_stiffness_mat, global_mass_vec);
    
    std::vector<FESystemUInt> sorted_ids;
    FESystem::Numerics::LocalVector<FESystemDouble> eig_vals;
    FESystem::Numerics::DenseMatrix<FESystemDouble> eig_vecs;
    
    
    staticAnalysis(dim, mesh, dof_map, nonbc_sparsity_pattern, old_to_new_id_map, nonbc_dofs, global_stiffness_mat, rhs, sol);
    
    modalAnalysis(dim, mesh, dof_map, nonbc_sparsity_pattern, old_to_new_id_map, nonbc_dofs, global_stiffness_mat, global_mass_vec, n_modes, sorted_ids, eig_vals, eig_vecs);
    
    transientAnalysis(dim, mesh, dof_map, nonbc_sparsity_pattern, old_to_new_id_map, nonbc_dofs, global_stiffness_mat, global_mass_vec, sorted_ids, eig_vals, eig_vecs);
    
    return 0;
}




int main(int argc, char * const argv[])
{
    return beam_analysis_driver(argc, argv);
}

