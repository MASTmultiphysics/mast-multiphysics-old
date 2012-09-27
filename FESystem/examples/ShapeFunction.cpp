//
//  ShapeFunctionTest.cpp
//  FESystem
//
//  Created by Manav Bhatia on 9/26/12.
//
//

// FESystem include
#include "TestingIncludes.h"


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



