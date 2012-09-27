//
//  EulerTest.cpp
//  FESystem
//
//  Created by Manav Bhatia on 9/26/12.
//
//

// FESystem include
#include "TestingIncludes.h"








void calculateEulerQuantities(FESystem::Mesh::ElementType elem_type, FESystemUInt n_elem_nodes, const FESystem::Base::DegreeOfFreedomMap& dof_map,
                              const FESystem::Mesh::MeshBase& mesh, 
                              const FESystem::Numerics::VectorBase<FESystemDouble>& sol,
                              const FESystem::Numerics::VectorBase<FESystemDouble>& vel,
                              FESystem::Numerics::VectorBase<FESystemDouble>& residual,
                              FESystem::Numerics::MatrixBase<FESystemDouble>& global_stiffness_mat,
                              FESystem::Numerics::VectorBase<FESystemDouble>& global_mass_vec)
{
    //FESystemDouble E=30.0e6, nu=0.3, rho=2700.0, p_val = 67.5, thick = 0.1; // (Reddy's parameters)
    FESystemDouble cp= 1.003e3, cv = 0.716e3, R=cp-cv, gamma=cp/cv; // (R should be 287.04)
    FESystemUInt n_elem_dofs;
    
    n_elem_dofs = 4*n_elem_nodes;
    
    
    FESystem::Numerics::DenseMatrix<FESystemDouble> elem_mat1, elem_mat2;
    FESystem::Numerics::LocalVector<FESystemDouble> elem_vec, elem_sol, elem_vel;
    std::vector<FESystemUInt> elem_dof_indices;
    elem_mat1.resize(n_elem_dofs, n_elem_dofs); elem_mat2.resize(n_elem_dofs, n_elem_dofs); elem_vec.resize(n_elem_dofs);  elem_sol.resize(n_elem_dofs); elem_vel.resize(n_elem_dofs);
    
    // prepare the quadrature rule and FE for the
    FESystem::Quadrature::TrapezoidQuadrature q_rule;
    FESystem::FiniteElement::FELagrange fe;
    FESystem::Fluid::FluidElementBase fluid_elem;
    
    switch (elem_type)
    {
        case FESystem::Mesh::QUAD4:
        case FESystem::Mesh::TRI3:
            q_rule.init(2, 3);  // bending quadrature is higher than the shear quadrature for reduced integrations
            break;
            
        case FESystem::Mesh::QUAD9:
        case FESystem::Mesh::TRI6:
            q_rule.init(2, 5);
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
    }
    
    residual.zero(); global_stiffness_mat.zero(); global_mass_vec.zero();
    
    const std::vector<FESystem::Mesh::ElemBase*>& elems = mesh.getElements();
    
    
    for (FESystemUInt i=0; i<elems.size(); i++)
    {
        fe.clear();
        fe.reinit(*(elems[i]));
        elem_mat1.zero();
        elem_mat2.zero();
        elem_vec.zero();
        
        dof_map.getFromGlobalVector(*(elems[i]), sol, elem_sol);
        dof_map.getFromGlobalVector(*(elems[i]), vel, elem_vel);
        
        fluid_elem.clear();
        fluid_elem.initialize(*(elems[i]), fe, q_rule, cp, cv);

        fluid_elem.calculateResidualVector(elem_sol, elem_vel, elem_vec);
        fluid_elem.calculateTangentMatrix(elem_sol, elem_vel, elem_mat1, elem_mat2);
        
        dof_map.addToGlobalVector(*(elems[i]), elem_vec, residual);
        dof_map.addToGlobalMatrix(*(elems[i]), elem_mat1, global_stiffness_mat);
    }
}



void nonlinearEulerSolution(FESystemUInt dim, FESystem::Mesh::ElementType elem_type, FESystemUInt n_elem_nodes, const
                            FESystem::Mesh::MeshBase& mesh, const FESystem::Base::DegreeOfFreedomMap& dof_map,
                            const FESystem::Numerics::SparsityPattern& nonbc_sparsity_pattern,
                            const std::map<FESystemUInt, FESystemUInt>& old_to_new_id_map, const std::vector<FESystemUInt>& nonbc_dofs)
{
    const std::vector<FESystem::Mesh::Node*>& nodes = mesh.getNodes();
    
    FESystem::Numerics::SparseMatrix<FESystemDouble> stiff_mat;
    FESystem::Numerics::LocalVector<FESystemDouble> residual, sol, vel, mass;
    
    residual.resize(dof_map.getNDofs()); sol.resize(dof_map.getNDofs()); vel.resize(dof_map.getNDofs()); mass.resize(dof_map.getNDofs());
    stiff_mat.resize(dof_map.getSparsityPattern());
    
    FESystem::LinearSolvers::LUFactorizationLinearSolver<FESystemDouble> linear_solver;
    FESystem::NonlinearSolvers::NewtonIterationNonlinearSolver<FESystemDouble> nonlinear_solver;
    
    nonlinear_solver.initialize(stiff_mat, linear_solver);
    nonlinear_solver.setConvergenceLimits(20, 1.0e-10);
    
    // for output
    FESystem::OutputProcessor::VtkOutputProcessor output;
    std::fstream output_file;
    std::vector<FESystemUInt> vars(3); vars[0]=0; vars[1]=1; vars[2]=2; // write all solutions
    
    
    FESystem::NonlinearSolvers::NonlinearSolverCallBack call_back = nonlinear_solver.getCurrentCallBack();
    
    while (call_back != FESystem::NonlinearSolvers::SOLUTION_CONVERGED)
    {
        switch (call_back)
        {
            case FESystem::NonlinearSolvers::WAITING_TO_START:
                // nothing to be done
                break;
                
            case FESystem::NonlinearSolvers::SET_INITIAL_GUESS:
            {
                // initialize solution to a uniform density and velocity value
                FESystem::Numerics::VectorBase<FESystemDouble>& vec = nonlinear_solver.getCurrentSolution();
                vec.zero();
                for (FESystemUInt i=0; i<nodes.size(); i++)
                {
                    vec.setVal(nodes[i]->getDegreeOfFreedomUnit(0).global_dof_id[0], 1.05); // rho
                    vec.setVal(nodes[i]->getDegreeOfFreedomUnit(1).global_dof_id[0], 1.05 * 300.0); // rho u
                    vec.setVal(nodes[i]->getDegreeOfFreedomUnit(3).global_dof_id[0], 1.05 * (300.0*716.0 + 300.0*300.0*0.5)); // rho * (cv * T + u^2/2)
                }
            }
                break;
                
            case FESystem::NonlinearSolvers::EVALUATE_RESIDUAL:
            {
                calculateEulerQuantities(elem_type, n_elem_nodes, dof_map, mesh, nonlinear_solver.getCurrentSolution(), nonlinear_solver.getCurrentSolution(), // the last arg is actually supposed to be velocity
                                         nonlinear_solver.getResidualVector(), nonlinear_solver.getJacobianMatrix(), mass);
                
                // write the solution
                std::stringstream oss;
                oss << "sol_" << nonlinear_solver.getCurrentIterationNumber() << ".vtk";
                output_file.open(oss.str().c_str(),std::fstream::out);
                output.writeMesh(output_file, mesh, dof_map);
                output.writeSolution(output_file, "Sol", mesh, dof_map, vars, sol);
                output_file.close();
            }
                break;
                
            case FESystem::NonlinearSolvers::EVALUATE_JACOBIAN:
            {
                // zero the solution vectors
                calculateEulerQuantities(elem_type, n_elem_nodes, dof_map, mesh, nonlinear_solver.getCurrentSolution(), nonlinear_solver.getCurrentSolution(), // the last arg is actually supposed to be velocity
                                         nonlinear_solver.getResidualVector(), nonlinear_solver.getJacobianMatrix(), mass);
            }
                break;
                
            case FESystem::NonlinearSolvers::EVALUATE_RESIDUAL_AND_JACOBIAN:
            {
                // zero the solution vectors
                calculateEulerQuantities(elem_type, n_elem_nodes, dof_map, mesh, nonlinear_solver.getCurrentSolution(), nonlinear_solver.getCurrentSolution(), // the last arg is actually supposed to be velocity
                                         nonlinear_solver.getResidualVector(), nonlinear_solver.getJacobianMatrix(), mass);
            }
                break;
                
            default:
                break;
        }
        call_back = nonlinear_solver.incrementSolution();
        if (call_back == FESystem::NonlinearSolvers::MAXIMUM_ITERATIONS_REACHED)
        {
            std::cout << "Maximum iterations reached" << std::endl;
            break;
        }
    }
    
    // copy the last solution and output
    sol.copyVector(nonlinear_solver.getCurrentSolution());
    // write the solution for each node
    for (FESystemUInt i=0; i<nodes.size(); i++)
    {
        std::cout << "Node: " << std::setw(8) << i;
        // write location
        for (FESystemUInt j=0; j<dim; j++)
            std::cout << std::setw(15) << nodes[i]->getVal(j);
//        // write the forces
//        for (FESystemUInt j=0; j<3; j++)
//            std::cout << std::setw(15) << rhs.getVal(nodes[i]->getDegreeOfFreedomUnit(j).global_dof_id[0]);
        for (FESystemUInt j=0; j<3; j++)
            std::cout << std::setw(15) << sol.getVal(nodes[i]->getDegreeOfFreedomUnit(j).global_dof_id[0]);
        std::cout << std::endl;
    }
}




int euler_analysis_driver(int argc, char * const argv[])
{
    FESystem::Mesh::ElementType elem_type;
    
    // create the mesh object
    FESystem::Mesh::MeshBase mesh;
    
    // create a nx x ny grid of nodes
    FESystemUInt nx, ny, dim, n_elem_nodes, n_elem_dofs, n_modes;
    FESystemDouble x_length, y_length;
    FESystem::Geometry::Point origin(3);
    
    nx=11; ny=11; x_length = 10.8; y_length = 10.8; dim = 2; n_modes = 10;
    elem_type = FESystem::Mesh::QUAD4;
    createPlaneMesh(elem_type, mesh, origin, nx, ny, x_length, y_length, n_elem_nodes, CROSS);
    
    n_elem_dofs = 6*n_elem_nodes;
    
    // set the location of individual nodes
    const std::vector<FESystem::Mesh::Node*>& nodes = mesh.getNodes();
    
    // now add the degrees of freedom
    FESystem::Base::DegreeOfFreedomMap dof_map(mesh);
    std::string name;
    name = "rho"; dof_map.addVariable(name, 0);
    name = "rhou"; dof_map.addVariable(name, 0);
    name = "rhov"; dof_map.addVariable(name, 0);
    //name = "rhow"; dof_map.addVariable(name, 0);
    name = "rhoe"; dof_map.addVariable(name, 0);
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
//        bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(0).global_dof_id[0]); // u-displacement
//        bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(1).global_dof_id[0]); // v-displacement
//        bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(5).global_dof_id[0]); // theta-z
//        if ((nodes[i]->getVal(0) == 0.0) || (nodes[i]->getVal(1) == 0.0) || (nodes[i]->getVal(0) == x_length) || (nodes[i]->getVal(1) == y_length)) // boundary nodes
//        {
//            //bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(0).global_dof_id[0]); // u-displacement on edges
//            //bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(1).global_dof_id[0]); // v-displacement on edges
//            bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(2).global_dof_id[0]); // w-displacement on edges
//        }
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
    
    nonlinearEulerSolution(dim, elem_type, n_elem_nodes, mesh, dof_map, nonbc_sparsity_pattern, old_to_new_id_map, nonbc_dofs);
    
    exit(1);
    
    return 0;
}


