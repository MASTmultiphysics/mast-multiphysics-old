////
////  FluidStructureInteraction.cpp
////  FESystem
////
////  Created by Manav Bhatia on 1/10/13.
////
////
//
//
//// FESystem include
//#include "TestingIncludes.h"
//
//const FESystemDouble  rho=1.05, u1=277.832, temp = 300.0, cp= 1.003e3, cv = 0.716e3, R=cp-cv, q0 = 0.5*rho*u1*u1, p = R*rho*temp, time_step=1.0e-4, final_t=1.0e6;
//
//
//void setBoundaryConditionTag(FESystem::Mesh::MeshBase& fluid_mesh, FESystem::Mesh::MeshBase& structural_mesh, std::set<FESystemUInt>& structural_bc_dofs)
//{
//    // first the fluid system
//    {
//        std::vector<FESystem::Mesh::Node*>& nodes = fluid_mesh.getNodes();
//        
//        for (FESystemUInt i=0; i<nodes.size(); i++)
//        {
//            if ((nodes[i]->getVal(0) == 0.0)) // left boundary nodes
//            {
//                const std::set<FESystem::Mesh::ElemBase*>& e_set = nodes[i]->getElementConnectivitySet();
//                std::set<FESystem::Mesh::ElemBase*>::const_iterator it = e_set.begin(), end = e_set.end();
//                for ( ; it != end; it++)
//                    if (!(*it)->checkForTag(0))
//                        (*it)->setTag(0);
//            }
//            
//            if (fabs(nodes[i]->getVal(0) - x_length) <= 100*FESystem::Base::getMachineEpsilon<FESystemDouble>()) // right boundary nodes
//            {
//                const std::set<FESystem::Mesh::ElemBase*>& e_set = nodes[i]->getElementConnectivitySet();
//                std::set<FESystem::Mesh::ElemBase*>::const_iterator it = e_set.begin(), end = e_set.end();
//                for ( ; it != end; it++)
//                    if (!(*it)->checkForTag(1))
//                        (*it)->setTag(1);
//            }
//            
//            if ((nodes[i]->getVal(1) == 0.0)) // lower boundary nodes
//            {
//                const std::set<FESystem::Mesh::ElemBase*>& e_set = nodes[i]->getElementConnectivitySet();
//                std::set<FESystem::Mesh::ElemBase*>::const_iterator it = e_set.begin(), end = e_set.end();
//                for ( ; it != end; it++)
//                    if (!(*it)->checkForTag(2))
//                        (*it)->setTag(2);
//                
//                // set tag 4 for airfoil surface
//                it = e_set.begin(), end = e_set.end();
//                if ((nodes[i]->getVal(0)>x0) && (nodes[i]->getVal(0)<x1))
//                    for ( ; it != end; it++)
//                        if (!(*it)->checkForTag(4))
//                            (*it)->setTag(4);
//            }
//            
//            
//            if (fabs(nodes[i]->getVal(1) - y_length) <= 100*FESystem::Base::getMachineEpsilon<FESystemDouble>()) // upper boundary nodes
//            {
//                const std::set<FESystem::Mesh::ElemBase*>& e_set = nodes[i]->getElementConnectivitySet();
//                std::set<FESystem::Mesh::ElemBase*>::const_iterator it = e_set.begin(), end = e_set.end();
//                for ( ; it != end; it++)
//                    if (!(*it)->checkForTag(3))
//                        (*it)->setTag(3);
//            }
//        }
//    }
//    
//    // next, the structural system
//    {
//        // since beam structure, constraint all w, tx and ty dofs
//        std::vector<FESystem::Mesh::Node*>& nodes = structural_mesh.getNodes();
//        for (FESystemUInt i=0; i<nodes.size(); i++)
//        {
//            bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(0).global_dof_id[0]); // u-disp
//            bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(1).global_dof_id[0]); // v-disp
//            bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(3).global_dof_id[0]); // theta-x
//            bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(5).global_dof_id[0]); // theta-z
//            
//        }
//        
//        // simply supported on the left and right dofs
//        const FESystem::Mesh::Node
//        &n1 = structural_mesh.getNodeFromInternalID(0), // left boundary node on panel
//        &n2 = structural_mesh.getNodeFromInternalID(structural_mesh.getNNodes()); // right boundary node on panel
//        
//        if (if_clamped)
//            for (FESystem)
//    }
//
//}
//
//
//
//void modifyFluidMeshForStructuralSolution()
//{
//    {
//        std::vector<FESystem::Mesh::Node*>& nodes = mesh.getNodes();
//        std::vector<FESystem::Mesh::Node*>::iterator it=nodes.begin(), end=nodes.end();
//        
//        FESystemDouble x_val, y_val;
//        
//        for ( ; it!=end; it++)
//        {
//            if (((*it)->getVal(0) >= x0) && ((*it)->getVal(0) <= x1))
//            {
//                x_val = (*it)->getVal(0);
//                y_val = (*it)->getVal(1);
//                
//                y_val += thickness*(1.0-y_val/y_length)*sin(PI_VAL*(x_val-x0)/chord);
//                
//                (*it)->setVal(1, y_val);
//            }
//        }
//        
//    }
//    
//    
//    // now update the elements after mesh modification
//    {
//        std::vector<FESystem::Mesh::ElemBase*> elems = mesh.getElements();
//        std::vector<FESystem::Mesh::ElemBase*>::iterator it=elems.begin(), end=elems.end();
//        
//        for ( ; it!=end; it++)
//            (*it)->updateAfterMeshDeformation();
//    }
//}
//
//
//void transientFluidStructureInteractions(const FESystemUInt dim, const FESystem::Mesh::ElementType fluid_elem_type, const FESystemUInt n_fluid_elem_nodes, const FESystem::Mesh::ElementType struct_elem_type, const FESystemUInt n_struct_elem_nodes,
//                                         FESystem::Mesh::MeshBase& fluid_mesh, FESystem::DegreeOfFreedom::DegreeOfFreedomMap& fluid_dof_map, FESystem::Mesh::MeshBase& structural_mesh, FESystem::DegreeOfFreedom::DegreeOfFreedomMap& structural_dof_map,
//                                         const std::set<FESystemUInt>& structural_bc_dofs)
//{
//    FESystem::LinearSolvers::LUFactorizationLinearSolver<FESystemDouble> structural_linear_solver, fluid_linear_solver;
//    const std::vector<FESystem::Mesh::Node*>& fluid_nodes = fluid_mesh.getNodes();
//    
//    FESystem::Numerics::SparseMatrix<FESystemDouble> fluid_stiff_mat, fluid_mass, structural_stiff_mat, structural_mass, structural_ode_jac, structural_ode_mass;
//    FESystem::Numerics::LocalVector<FESystemDouble>  fluid_sol, fluid_vel, structural_sol, structural_vel, primitive_sol, additional_sol;
//    
//    fluid_sol.resize(fluid_dof_map.getNDofs()); fluid_vel.resize(fluid_dof_map.getNDofs());
//    structural_sol.resize(structural_dof_map.getNDofs()); structural_vel.resize(structural_dof_map.getNDofs());
//    fluid_stiff_mat.resize(fluid_dof_map.getSparsityPattern());  fluid_mass.resize(fluid_dof_map.getSparsityPattern());
//    structural_stiff_mat.resize(structural_dof_map.getSparsityPattern());  structural_mass.resize(structural_dof_map.getSparsityPattern());
//    primitive_sol.resize(fluid_dof_map.getNDofs()); additional_sol.resize(fluid_dof_map.getNDofs());
//    
//    // set the initial condition
//    fluid_sol.zero();
//    structural_sol.zero();
//    for (FESystemUInt i=0; i<nodes.size(); i++)
//    {
//        fluid_sol.setVal(fluid_nodes[i]->getDegreeOfFreedomUnit(0).global_dof_id[0], rho);
//        fluid_sol.setVal(fluid_nodes[i]->getDegreeOfFreedomUnit(1).global_dof_id[0], rho * u1);
//        fluid_sol.setVal(fluid_nodes[i]->getDegreeOfFreedomUnit(3).global_dof_id[0], rho * (temp*cv + u1*u1*0.5));
//    }
//    
//    // initialize the solver
//    FESystem::TransientSolvers::NewmarkTransientSolver<FESystemDouble> fluid_transient_solver, structural_transient_solver;
//    FESystem::Numerics::SparsityPattern fluid_ode_sparsity, structural_ode_sparsity;
//    std::vector<FESystemBoolean> fluid_ode_order_include(1), structural_ode_order_include(2);
//    fluid_ode_order_include[0] = true; structural_ode_order_include[0] = true; structural_ode_order_include[1] = false;
//    
//    std::vector<FESystemDouble> fluid_int_constants(1), structural_int_constants(2);
//    fluid_int_constants[0]=0.5; structural_int_constants[0]=1.0/4.0; structural_int_constants[1]=1.0/2.0;
//
//    fluid_transient_solver.initialize(1, fluid_dof_map.getNDofs(), fluid_int_constants);
//    //fluid_transient_solver.enableAdaptiveTimeStepping(4, 0.4, 1.0e3);
//    fluid_transient_solver.setConvergenceTolerance(nonlin_tol, max_nonlin_iters);
//    fluid_transient_solver.setActiveJacobianTerm(fluid_ode_order_include);
//    fluid_transient_solver.setMassMatrix(false, &fluid_mass);
//    fluid_transient_solver.setInitialTimeData(0, time_step, fluid_sol);
//    fluid_transient_solver.setJacobianMatrix(fluid_stiff_mat);
//    fluid_transient_solver.setLinearSolver(fluid_linear_solver, false);
//    fluid_transient_solver.setLinearSolverDataStructureReuse(true);
//
//    // structural
//    structural_transient_solver.initialize(1, structural_dof_map.getNDofs(), structural_int_constants);
//    //structural_transient_solver.enableAdaptiveTimeStepping(4, 0.4, 1.0e3);
//    structural_transient_solver.setConvergenceTolerance(nonlin_tol, max_nonlin_iters);
//    structural_transient_solver.setActiveJacobianTerm(structural_ode_order_include);
//    structural_transient_solver.initializeMatrixSparsityPatterForSystem(nonbc_structural_sparsity_pattern, structural_ode_sparsity)
//    structural_ode_jac.resize(structural_ode_sparsity); structural_ode_mass(structural_ode_sparsity);
//    structural_transient_solver.setMassMatrix(false, &structural_ode_mass);
//    structural_transient_solver.setInitialTimeData(0, time_step, structural_sol);
//    structural_transient_solver.setJacobianMatrix(structural_ode_jac);
//    structural_transient_solver.setLinearSolver(structural_linear_solver, false);
//    structural_transient_solver.setLinearSolverDataStructureReuse(true);
//    
//    
//    FESystem::OutputProcessor::VtkOutputProcessor output;
//    std::fstream output_file;
//    std::vector<FESystemUInt> vars(4); vars[0]=0; vars[1]=1; vars[2]=2; vars[3] = 3; //vars[4] = 4; // write all solutions
//    
//    FESystemUInt n_skip=0, n_count=0, n_write=0;
//    FESystem::TransientSolvers::TransientSolverCallBack fluid_call_back, structural_call_back;
//    FESystemBoolean if_converged = false, if_time_step_converged=false;
//    
//    while (!if_converged)
//    {
//        // if both have converged
//        if ((structural_call_back == FESystem::TransientSolvers::TIME_STEP_CONVERGED) && (fluid_call_back == FESystem::TransientSolvers::TIME_STEP_CONVERGED))
//        {
//            // write the solution
//            if (n_count == n_skip)
//            {
//                std::stringstream oss;
//                oss << "sol_" << n_write << ".vtk";
//                output_file.open(oss.str().c_str(),std::fstream::out);
//                output.writeMesh(output_file, mesh, dof_map);
//                output.writeSolution(output_file, "Sol", mesh, dof_map, vars, transient_solver.getCurrentStateVector());
//                output.writeSolution(output_file, "Vel", mesh, dof_map, vars, transient_solver.getCurrentStateVelocityVector());
//                output.writeSolution(output_file, "Primitive", mesh, dof_map, vars, primitive_sol);
//                output.writeSolution(output_file, "Additional", mesh, dof_map, vars, additional_sol);
//                additional_sol.setAllVals(transient_solver.getCurrentTime());
//                output.writeSolution(output_file, "TimeValue", mesh, dof_map, vars, additional_sol);
//                
//                output_file.close();
//                
//                n_write++;
//                n_count=0;
//            }
//            else
//                n_count++;
//        }
//        else // handle the individual callbacks independentl
//        {
//            switch (fluid_call_back)
//            {
//                case FESystem::TransientSolvers::EVALUATE_X_DOT:
//                {
//                    calculateEulerQuantities(elem_type, n_elem_nodes, dof_map, mesh, transient_solver.getCurrentStepSize(), transient_solver.getCurrentIterationNumber(),
//                                             transient_solver.getCurrentStateVector(), transient_solver.getCurrentStateVelocityVector(), false,
//                                             transient_solver.getCurrentVelocityFunctionVector(), transient_solver.getCurrentJacobianMatrix(),
//                                             mass, primitive_sol, additional_sol);
//                }
//                    break;
//                    
//                case FESystem::TransientSolvers::EVALUATE_X_DOT_AND_X_DOT_JACOBIAN:
//                {
//                    calculateEulerQuantities(elem_type, n_elem_nodes, dof_map, mesh, transient_solver.getCurrentStepSize(), transient_solver.getCurrentIterationNumber(),
//                                             transient_solver.getCurrentStateVector(), transient_solver.getCurrentStateVelocityVector(), true,
//                                             transient_solver.getCurrentVelocityFunctionVector(), transient_solver.getCurrentJacobianMatrix(),
//                                             mass, primitive_sol, additional_sol);
//                }
//                    break;
//                    
//                default:
//                    break;
//            }
//            
//            switch (structural_call_back)
//            {
//                case FESystem::TransientSolvers::EVALUATE_X_DOT:
//                {
//                    calculateEulerQuantities(elem_type, n_elem_nodes, dof_map, mesh, transient_solver.getCurrentStepSize(), transient_solver.getCurrentIterationNumber(),
//                                             transient_solver.getCurrentStateVector(), transient_solver.getCurrentStateVelocityVector(), false,
//                                             transient_solver.getCurrentVelocityFunctionVector(), transient_solver.getCurrentJacobianMatrix(),
//                                             mass, primitive_sol, additional_sol);
//                }
//                    break;
//                    
//                case FESystem::TransientSolvers::EVALUATE_X_DOT_AND_X_DOT_JACOBIAN:
//                {
//                    calculateEulerQuantities(elem_type, n_elem_nodes, dof_map, mesh, transient_solver.getCurrentStepSize(), transient_solver.getCurrentIterationNumber(),
//                                             transient_solver.getCurrentStateVector(), transient_solver.getCurrentStateVelocityVector(), true,
//                                             transient_solver.getCurrentVelocityFunctionVector(), transient_solver.getCurrentJacobianMatrix(),
//                                             mass, primitive_sol, additional_sol);
//                }
//                    break;
//                    
//                default:
//                    break;
//            }
//
//            // now increment the call backs
//            if (fluid_call_back != FESystem::TransientSolvers::TIME_STEP_CONVERGED)
//                fluid_call_back = fluid_transient_solver.incrementTimeStep();
//            if (structural_call_back != FESystem::TransientSolvers::TIME_STEP_CONVERGED)
//                structural_call_back = structural_transient_solver.incrementTimeStep();
//            
//        }
//        if (fluid_transient_solver.getCurrentTime()>final_t)
//            if_converged = true;
//    }
//    
//}
//
//
//int fluidStructutreInteractionDriver(int argc, char * const argv[])
//{
//    // create the mesh object
//    FESystem::Mesh::MeshBase fluid_mesh, structural_mesh;
//    
//    // create a nx x ny grid of nodes
//    FESystemUInt n_fluid_elem_nodes, n_fluid_elem_dofs, n_struct_elem_nodes, n_struct_elem_dofs;
//    FESystem::Geometry::Point origin(3);
//    
//    FESystem::Mesh::ElementType fluid_elem_type = FESystem::Mesh::QUAD4, struct_elem_type = FESystem::Mesh::EDGE2;
//    
//    
//    // fluid mesh parameters and mesh initialization
//    {
//        FESystemUInt x_n_divs, y_n_divs;
//        std::vector<FESystemDouble> x_div_locations, x_relative_mesh_size_in_div, x_points, y_div_locations, y_relative_mesh_size_in_div, y_points;
//        std::vector<FESystemUInt> x_n_subdivs_in_div, y_n_subdivs_in_div;
//
//        // create fluid mesh and set BC tags
//        x_n_divs = 3;
//        x_div_locations.resize(x_n_divs+1);
//        x_relative_mesh_size_in_div.resize(x_n_divs+1);
//        x_n_subdivs_in_div.resize(x_n_divs);
//        
//        x_div_locations[0] = 0.0;
//        x_div_locations[1] = (x_length-chord)/2.0;
//        x_div_locations[2] = (x_length+chord)/2.0;
//        x_div_locations[3] = x_length;
//        x_relative_mesh_size_in_div[0] = 100.0;
//        x_relative_mesh_size_in_div[1] = 1.0;
//        x_relative_mesh_size_in_div[2] = 1.0;
//        x_relative_mesh_size_in_div[3] = 100.0;
//        
//        x_n_subdivs_in_div[0] = 30;
//        x_n_subdivs_in_div[1] = 40;
//        x_n_subdivs_in_div[2] = 30;
//        
//        y_n_divs = 1;
//        y_div_locations.resize(y_n_divs+1);
//        y_relative_mesh_size_in_div.resize(y_n_divs+1);
//        y_n_subdivs_in_div.resize(y_n_divs);
//        
//        y_div_locations[0] = 0.0;
//        y_div_locations[1] = y_length;
//        
//        y_relative_mesh_size_in_div[0] = 1.0;
//        y_relative_mesh_size_in_div[1] = 50.0;
//        
//        y_n_subdivs_in_div[0] = 40;
//
//        
//        distributePoints(x_n_divs, x_div_locations, x_n_subdivs_in_div, x_relative_mesh_size_in_div, x_points);
//        distributePoints(y_n_divs, y_div_locations, y_n_subdivs_in_div, y_relative_mesh_size_in_div, y_points);
//        
//        createPlaneMesh(fluid_elem_type, fluid_mesh, origin, x_points, y_points, n_fluid_elem_nodes, CROSS, true);
//        
//        n_fluid_elem_dofs = 4*n_fluid_elem_nodes;
//        
//        // structural mesh is based on the fluid panel mesh to create a one-to-one correspondence between nodes
//        y_div_locations[0] = x_div_locations[1];
//        y_div_locations[1] = x_div_locations[2];
//        
//        y_relative_mesh_size_in_div[0] = 1.0;
//        y_relative_mesh_size_in_div[1] = 1.0;
//        
//        y_n_subdivs_in_div[0] = x_n_subdivs_in_div[1];
//        
//        y_points.clear();
//        
//        distributePoints(y_n_divs, y_div_locations, y_n_subdivs_in_div, y_relative_mesh_size_in_div, y_points);
//        createLineMesh(struct_elem_type, structural_mesh, origin, y_points, n_struct_elem_nodes, CROSS, true);
//     
//        n_struct_elem_dofs = 6*n_struct_elem_nodes;
//    }
//
//    
//
//    
//    // now add the degrees of freedom
//    FESystem::DegreeOfFreedom::DegreeOfFreedomMap fluid_dof_map(fluid_mesh), structural_dof_map(structural_mesh);
//    std::string name;
//    name = "rho"; fluid_dof_map.addVariable(name, 0);
//    name = "rhou"; fluid_dof_map.addVariable(name, 0);
//    name = "rhov"; fluid_dof_map.addVariable(name, 0);
//    name = "rhoe"; fluid_dof_map.addVariable(name, 0);
//    fluid_dof_map.reinit(); // distribute the dofs
//
//    
//    name = "u"; structural_mesh.addVariable(name, 0);
//    name = "v"; structural_mesh.addVariable(name, 0);
//    name = "w"; structural_mesh.addVariable(name, 0);
//    name = "tx"; structural_mesh.addVariable(name, 0);
//    name = "ty"; structural_mesh.addVariable(name, 0);
//    name = "tz"; structural_mesh.addVariable(name, 0);
//    structural_dof_map.reinit(); // distribute the dofs
//
//    
//    // apply boundary condition and place a load on the last dof
//    std::set<FESystemUInt> structural_bc_dofs;
//    setBoundaryConditionTag(fluid_mesh, structural_mesh, structural_bc_dofs);
//
//    // prepare a map of the old to new ID
//    std::vector<FESystemUInt> structural_nonbc_dofs;
//    for (FESystemUInt i=0; i<structural_dof_map.getNDofs(); i++)
//        if (!structural_bc_dofs.count(i))
//            structural_nonbc_dofs.push_back(i);
//    
//    std::vector<FESystemUInt>::const_iterator dof_it=structural_nonbc_dofs.begin(), dof_end=nonbc_dofs.end();
//    std::map<FESystemUInt, FESystemUInt> structural_old_to_new_id_map;
//    FESystemUInt n=0;
//    for ( ; dof_it!=dof_end; dof_it++)
//        structural_old_to_new_id_map.insert(std::map<FESystemUInt, FESystemUInt>::value_type(*dof_it, n++));
//    
//    FESystem::Numerics::SparsityPattern nonbc_structural_sparsity_pattern;
//    structural_dof_map.getSparsityPattern().initSparsityPatternForNonConstrainedDOFs(structural_nonbc_dofs, structural_old_to_new_id_map, nonbc_structural_sparsity_pattern);
//
//    // transient solution
//    transientFluidStructureInteraction(dim, fluid_elem_type, n_fluid_elem_nodes, struct_elem_type, n_struct_elem_nodes, fluid_mesh, fluid_dof_map, structural_mesh, structural_dof_map, nonbc_structural_sparsity_pattern);
//
//    return 1;
//}
//
