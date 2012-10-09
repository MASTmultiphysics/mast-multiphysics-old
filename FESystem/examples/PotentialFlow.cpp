//
//  PotentialFlowTest.cpp
//  FESystem
//
//  Created by Manav Bhatia on 9/26/12.
//
//

// FESystem include
#include "TestingIncludes.h"


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
    
    createPlaneMesh(elem_type, mesh, origin, nx, ny, x_length, y_length, n_elem_nodes, CROSS, true);
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
    FESystem::TransientSolvers::NewmarkTransientSolver<FESystemDouble> transient_solver;
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
                
                transient_solver.updateVectorValuesForDerivativeOrder(1, rhs, transient_solver.getVelocityFunction());
                transient_solver.copyDerivativeValuesFromStateToVelocityVector(transient_solver.getCurrentStateVector(), transient_solver.getVelocityFunction());
            }
                break;
                
            default:
                break;
        }
    }
    
    return 0;
}

