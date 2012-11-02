//
//  EulerTest.cpp
//  FESystem
//
//  Created by Manav Bhatia on 9/26/12.
//
//

// FESystem include
#include "TestingIncludes.h"




const FESystemDouble  rho=1.05, u1=1300.0, temp = 300.0, cp= 1.003e3, cv = 0.716e3, R=cp-cv, p = R*rho*temp, time_step=1.0e-2, final_t=time_step*1.0e5;
const FESystemDouble x_length = 2, y_length = .5, t_by_c = 0.02, chord = 0.5, thickness = 0.5*t_by_c*chord, x0=x_length/2-chord/2, x1=x0+chord, nonlin_tol = 1.0e-10;
const FESystemUInt nx=200, ny=70, dim = 2, max_nonlin_iters = 10;



void calculateEulerQuantities(FESystem::Mesh::ElementType elem_type, FESystemUInt n_elem_nodes, const FESystem::Base::DegreeOfFreedomMap& dof_map,
                              const FESystem::Mesh::MeshBase& mesh,
                              const FESystem::Numerics::VectorBase<FESystemDouble>& sol,
                              const FESystem::Numerics::VectorBase<FESystemDouble>& vel,
                              FESystem::Numerics::VectorBase<FESystemDouble>& residual,
                              FESystem::Numerics::MatrixBase<FESystemDouble>& global_stiffness_mat,
                              FESystem::Numerics::MatrixBase<FESystemDouble>& global_mass_mat)
{
    FESystemUInt n_elem_dofs;
    
    n_elem_dofs = 4*n_elem_nodes;
    
    
    FESystem::Numerics::DenseMatrix<FESystemDouble> elem_mat1, elem_mat2, momentum_flux_tensor;
    FESystem::Numerics::LocalVector<FESystemDouble> elem_vec, elem_sol, elem_vel, bc_vec, mass_flux, energy_flux;
    std::vector<FESystemUInt> elem_dof_indices;
    elem_mat1.resize(n_elem_dofs, n_elem_dofs); elem_mat2.resize(n_elem_dofs, n_elem_dofs); elem_vec.resize(n_elem_dofs);
    elem_sol.resize(n_elem_dofs); elem_vel.resize(n_elem_dofs); bc_vec.resize(n_elem_dofs);
    mass_flux.resize(2); energy_flux.resize(2); momentum_flux_tensor.resize(2, 2);
    
    // prepare the quadrature rule and FE for the
    FESystem::Quadrature::TrapezoidQuadrature q_rule, q_boundary;
    FESystem::FiniteElement::FELagrange fe;
    FESystem::Fluid::FluidElementBase fluid_elem;
    
    switch (elem_type)
    {
        case FESystem::Mesh::QUAD4:
        case FESystem::Mesh::TRI3:
            q_rule.init(2, 2);  // bending quadrature is higher than the shear quadrature for reduced integrations
            q_boundary.init(1,2);
            break;
            
        case FESystem::Mesh::QUAD9:
        case FESystem::Mesh::TRI6:
            q_rule.init(2, 5);
            q_boundary.init(1, 5);
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
    }
    
    residual.zero(); global_stiffness_mat.zero(); global_mass_mat.zero();
    
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
        fluid_elem.initialize(*(elems[i]), fe, q_rule, time_step, cp, cv, elem_sol, elem_vel);

        // set the flux value for the left and right boundary
        mass_flux.zero(); mass_flux.setVal(0, rho*u1);
        momentum_flux_tensor.zero(); momentum_flux_tensor.setVal(0, 0, rho*u1*u1+p); momentum_flux_tensor.setVal(1, 1, p);
        energy_flux.zero(); energy_flux.setVal(0, rho*u1*(cv*temp+0.5*u1*u1)+p*u1);
                
        if (elems[i]->checkForTag(0)) // left edge
        {
            fluid_elem.calculateFluxBoundaryCondition(3, q_boundary, mass_flux, momentum_flux_tensor, energy_flux, bc_vec);
            dof_map.addToGlobalVector(*(elems[i]), bc_vec, residual);
        }
        if (elems[i]->checkForTag(1)) // right edge
        {
            fluid_elem.calculateFluxBoundaryConditionUsingLocalSolution(1, q_boundary, bc_vec);
            fluid_elem.calculateTangentMatrixForFluxBoundaryConditionUsingLocalSolution(1, q_boundary, elem_mat1);
            dof_map.addToGlobalVector(*(elems[i]), bc_vec, residual);
            dof_map.addToGlobalMatrix(*(elems[i]), elem_mat1, global_stiffness_mat);
        }
        // set the flux value for the lower and upper boundary
        if (elems[i]->checkForTag(2)) // lower edge
        {
            fluid_elem.calculateSolidWallFluxBoundaryCondition(0, q_boundary, bc_vec);
            fluid_elem.calculateTangentMatrixForSolidWallFluxBoundaryCondition(0, q_boundary, elem_mat1);
            dof_map.addToGlobalVector(*(elems[i]), bc_vec, residual);
            dof_map.addToGlobalMatrix(*(elems[i]), elem_mat1, global_stiffness_mat);
        }
        if (elems[i]->checkForTag(3)) // upper edge
        {
            fluid_elem.calculateSolidWallFluxBoundaryCondition(2, q_boundary, bc_vec);
            fluid_elem.calculateTangentMatrixForSolidWallFluxBoundaryCondition(2, q_boundary, elem_mat1);
            dof_map.addToGlobalVector(*(elems[i]), bc_vec, residual);
            dof_map.addToGlobalMatrix(*(elems[i]), elem_mat1, global_stiffness_mat);
        }

        
        fluid_elem.calculateResidualVector(elem_vec);
        fluid_elem.calculateTangentMatrix(elem_mat1, elem_mat2);
        
        dof_map.addToGlobalVector(*(elems[i]), elem_vec, residual);
        dof_map.addToGlobalMatrix(*(elems[i]), elem_mat1, global_stiffness_mat);
        dof_map.addToGlobalMatrix(*(elems[i]), elem_mat2, global_mass_mat);
    }
}




void testJacobian(FESystem::Mesh::ElementType elem_type, FESystemUInt n_elem_nodes, const FESystem::Base::DegreeOfFreedomMap& dof_map,
                  const FESystem::Mesh::MeshBase& mesh,
                  const FESystem::Numerics::VectorBase<FESystemDouble>& sol,
                  const FESystem::Numerics::VectorBase<FESystemDouble>& vel)
{
    FESystemUInt n_elem_dofs;
    
    n_elem_dofs = 4*n_elem_nodes;
    
    
    FESystem::Numerics::DenseMatrix<FESystemDouble> elem_mat1, elem_mat2, momentum_flux_tensor;
    FESystem::Numerics::LocalVector<FESystemDouble> elem_vec, elem_sol, elem_vel, bc_vec, mass_flux, energy_flux, elem_res, vec, elem_sol_delta;
    std::vector<FESystemUInt> elem_dof_indices;
    elem_mat1.resize(n_elem_dofs, n_elem_dofs); elem_mat2.resize(n_elem_dofs, n_elem_dofs); elem_vec.resize(n_elem_dofs);  elem_sol.resize(n_elem_dofs); elem_vel.resize(n_elem_dofs), bc_vec.resize(n_elem_dofs);
    mass_flux.resize(2); energy_flux.resize(2); momentum_flux_tensor.resize(2, 2); vec.resize(n_elem_dofs), elem_sol_delta.resize(n_elem_dofs); elem_res.resize(n_elem_dofs);
    
    // prepare the quadrature rule and FE for the
    FESystem::Quadrature::TrapezoidQuadrature q_rule, q_boundary;
    FESystem::FiniteElement::FELagrange fe;
    FESystem::Fluid::FluidElementBase fluid_elem;
    
    switch (elem_type)
    {
        case FESystem::Mesh::QUAD4:
        case FESystem::Mesh::TRI3:
            q_rule.init(2, 3);  // bending quadrature is higher than the shear quadrature for reduced integrations
            q_boundary.init(1,3);
            break;
            
        case FESystem::Mesh::QUAD9:
        case FESystem::Mesh::TRI6:
            q_rule.init(2, 5);
            q_boundary.init(1, 5);
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
    }
        
    const std::vector<FESystem::Mesh::ElemBase*>& elems = mesh.getElements();
    FESystemUInt elem_num=0;
    
    fe.clear();
    fe.reinit(*(elems[elem_num]));
    elem_mat1.zero();
    elem_mat2.zero();
    elem_vec.zero();
    
    dof_map.getFromGlobalVector(*(elems[elem_num]), sol, elem_sol);
    dof_map.getFromGlobalVector(*(elems[elem_num]), vel, elem_vel);

    // set the flux value for the left and right boundary
    mass_flux.zero(); mass_flux.setVal(0, rho*u1);
    momentum_flux_tensor.zero(); momentum_flux_tensor.setVal(0, 0, rho*u1*u1+p); momentum_flux_tensor.setVal(1, 1, p);
    energy_flux.zero(); energy_flux.setVal(0, rho*u1*(cv*temp+0.5*u1*u1)+p*u1);

    
    // calculate the baseline quantities
    elem_res.zero();
    fluid_elem.clear();
    fluid_elem.initialize(*(elems[elem_num]), fe, q_rule, time_step, cp, cv, elem_sol, elem_vel);
    
    if (elems[elem_num]->checkForTag(0)) // left edge
    {
        fluid_elem.calculateFluxBoundaryCondition(3, q_boundary, mass_flux, momentum_flux_tensor, energy_flux, bc_vec);
        elem_res.add(1.0, bc_vec);
    }
    if (elems[elem_num]->checkForTag(1)) // right edge
    {
        fluid_elem.calculateFluxBoundaryCondition(1, q_boundary, mass_flux, momentum_flux_tensor, energy_flux, bc_vec);
        elem_res.add(1.0, bc_vec);
    }
    // set the flux value for the lower and upper boundary
    if (elems[elem_num]->checkForTag(2)) // lower edge
    {
        fluid_elem.calculateFluxBoundaryCondition(0, q_boundary, mass_flux, momentum_flux_tensor, energy_flux, bc_vec);
        elem_res.add(1.0, bc_vec);
    }
    if (elems[elem_num]->checkForTag(3)) // upper edge
    {
        fluid_elem.calculateFluxBoundaryCondition(2, q_boundary, mass_flux, momentum_flux_tensor, energy_flux, bc_vec);
        elem_res.add(1.0, bc_vec);
    }
    
    fluid_elem.calculateResidualVector(elem_vec);
    elem_res.add(1.0, elem_vec);
    
    fluid_elem.calculateTangentMatrix(elem_mat1, elem_mat2); // mat1 is the jacobian to be compared
    
    FESystemDouble delta = 0.01;
    
    // now calculate the perturbed quantity
    for (FESystemUInt i=0; i<n_elem_dofs; i++)
    {
        elem_sol_delta.copyVector(elem_sol);
        elem_sol_delta.setVal(i, elem_sol.getVal(i)+delta);
        vec.zero();
        
        fluid_elem.clear();
        fluid_elem.initialize(*(elems[elem_num]), fe, q_rule, time_step, cp, cv, elem_sol_delta, elem_vel);
        
        if (elems[elem_num]->checkForTag(0)) // left edge
        {
            fluid_elem.calculateFluxBoundaryCondition(3, q_boundary, mass_flux, momentum_flux_tensor, energy_flux, bc_vec);
            vec.add(1.0, bc_vec);
        }
        if (elems[elem_num]->checkForTag(1)) // right edge
        {
            fluid_elem.calculateFluxBoundaryCondition(1, q_boundary, mass_flux, momentum_flux_tensor, energy_flux, bc_vec);
            vec.add(1.0, bc_vec);
        }
        // set the flux value for the lower and upper boundary
        if (elems[elem_num]->checkForTag(2)) // lower edge
        {
            fluid_elem.calculateFluxBoundaryCondition(0, q_boundary, mass_flux, momentum_flux_tensor, energy_flux, bc_vec);
            vec.add(1.0, bc_vec);
        }
        if (elems[elem_num]->checkForTag(3)) // upper edge
        {
            fluid_elem.calculateFluxBoundaryCondition(2, q_boundary, mass_flux, momentum_flux_tensor, energy_flux, bc_vec);
            vec.add(1.0, bc_vec);
        }
        
        fluid_elem.calculateResidualVector(elem_vec);
        vec.add(1.0, elem_vec);
        vec.add(-1.0, elem_res);
        vec.scale(1.0/delta);
        
        elem_mat2.setColumnVals(i, 0, n_elem_dofs-1, vec);
    }
    
    elem_mat1.write(std::cout);
    elem_mat2.write(std::cout);
}


void transientEulerAnalysis(FESystemUInt dim, FESystem::Mesh::ElementType elem_type, FESystemUInt n_elem_nodes, const FESystem::Mesh::MeshBase& mesh, const FESystem::Base::DegreeOfFreedomMap& dof_map,
                            const FESystem::Numerics::SparsityPattern& nonbc_sparsity_pattern, const std::map<FESystemUInt, FESystemUInt>& old_to_new_id_map, const std::vector<FESystemUInt>& nonbc_dofs)
{
    FESystem::LinearSolvers::LUFactorizationLinearSolver<FESystemDouble> linear_solver;
    const std::vector<FESystem::Mesh::Node*>& nodes = mesh.getNodes();
    
    FESystem::Numerics::SparseMatrix<FESystemDouble> stiff_mat, mass, stiff_mat_reduced, mass_mat_reduced;
    FESystem::Numerics::LocalVector<FESystemDouble>  sol, sol_reduced, vel, vel_func;
    
    sol.resize(dof_map.getNDofs()); vel_func.resize(dof_map.getNDofs()); vel.resize(dof_map.getNDofs()); sol_reduced.resize(nonbc_sparsity_pattern.getNDOFs());
    stiff_mat.resize(dof_map.getSparsityPattern());  mass.resize(dof_map.getSparsityPattern());
    stiff_mat_reduced.resize(nonbc_sparsity_pattern); mass_mat_reduced.resize(nonbc_sparsity_pattern);
    
    FESystemUInt id= 0;
    
    // initialize the solver
    FESystem::TransientSolvers::NewmarkTransientSolver<FESystemDouble> transient_solver;
    FESystem::Numerics::SparsityPattern ode_sparsity;
    std::vector<FESystemBoolean> ode_order_include(1); ode_order_include[0] = true;
    std::vector<FESystemDouble> int_constants(1); int_constants[0]=0.5;
    transient_solver.initialize(1, nonbc_dofs.size(), int_constants);
    transient_solver.setConvergenceTolerance(nonlin_tol, max_nonlin_iters);
    transient_solver.setActiveJacobianTerm(ode_order_include);
    transient_solver.setMassMatrix(false, &mass_mat_reduced);
    
    transient_solver.initializeStateVector(sol_reduced); //  initialize the vector and apply the initial condition
    // set the initial condition
    sol.zero();
    for (FESystemUInt i=0; i<nodes.size(); i++)
    {
        sol.setVal(nodes[i]->getDegreeOfFreedomUnit(0).global_dof_id[0], rho);
        sol.setVal(nodes[i]->getDegreeOfFreedomUnit(1).global_dof_id[0], rho * u1);
        sol.setVal(nodes[i]->getDegreeOfFreedomUnit(3).global_dof_id[0], rho * (temp*cv + u1*u1*0.5));
    }
    
    sol.getSubVectorValsFromIndices(nonbc_dofs, sol_reduced);
    
    transient_solver.setInitialTimeData(0, time_step, sol_reduced);
    
    transient_solver.initializeMatrixSparsityPatterForSystem(nonbc_sparsity_pattern, ode_sparsity);
    stiff_mat_reduced.resize(ode_sparsity);
    transient_solver.setJacobianMatrix(stiff_mat_reduced);
    transient_solver.setLinearSolver(linear_solver, false);
    
    FESystem::OutputProcessor::VtkOutputProcessor output;
    std::fstream output_file;
    std::vector<FESystemUInt> vars(4); vars[0]=0; vars[1]=1; vars[2]=2; vars[3] = 3; //vars[4] = 4; // write all solutions
    
    FESystemUInt n_skip=0, n_count=0, n_write=0;
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
                    sol.setSubVectorValsFromIndices(nonbc_dofs, transient_solver.getCurrentStateVector());

                    std::stringstream oss;
                    oss << "sol_" << n_write << ".vtk";
                    output_file.open(oss.str().c_str(),std::fstream::out);
                    output.writeMesh(output_file, mesh, dof_map);
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
                if (false)
                {
                    sol.copyVector(transient_solver.getCurrentStateVector());
                    // write the solution for each node
                    for (FESystemUInt i=0; i<nodes.size(); i++)
                    {
                        std::cout << "Node: " << std::setw(8) << i;
                        // write location
                        for (FESystemUInt j=0; j<3; j++)
                            std::cout << std::setw(15) << nodes[i]->getVal(j);
                        for (FESystemUInt j=0; j<4; j++)
                            std::cout << std::setw(15) << sol.getVal(nodes[i]->getDegreeOfFreedomUnit(j).global_dof_id[0]);
                        std::cout << std::endl;
                    }
                    std::cout << std::endl<<  std::endl;;

                }

                sol.setSubVectorValsFromIndices(nonbc_dofs, transient_solver.getCurrentStateVector());
                vel.setSubVectorValsFromIndices(nonbc_dofs, transient_solver.getCurrentStateVelocityVector());
                
                //testJacobian(elem_type, n_elem_nodes, dof_map, mesh, sol, vel, vel_func);
                calculateEulerQuantities(elem_type, n_elem_nodes, dof_map, mesh, sol, vel, vel_func, stiff_mat, mass);
                
                vel_func.getSubVectorValsFromIndices(nonbc_dofs, transient_solver.getCurrentVelocityFunctionVector());
                stiff_mat.getSubMatrixValsFromRowAndColumnIndices(nonbc_dofs, nonbc_dofs, old_to_new_id_map, transient_solver.getCurrentJacobianMatrix());
                mass.getSubMatrixValsFromRowAndColumnIndices(nonbc_dofs, nonbc_dofs, old_to_new_id_map, mass_mat_reduced);
                
                
                if (false)
                {
                    sol.copyVector(transient_solver.getCurrentStateVelocityVector());
                    // write the solution for each node
                    for (FESystemUInt i=0; i<nodes.size(); i++)
                    {
                        std::cout << "Node: " << std::setw(8) << i;
                        // write location
                        for (FESystemUInt j=0; j<3; j++)
                            std::cout << std::setw(15) << nodes[i]->getVal(j);
                        for (FESystemUInt j=0; j<4; j++)
                            std::cout << std::setw(15) << sol.getVal(nodes[i]->getDegreeOfFreedomUnit(j).global_dof_id[0]);
                        std::cout << std::endl;
                    }
                    std::cout << std::endl<<  std::endl;
                }
            }
                break;
                
            default:
                FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
                break;
        }
    }
}



void nonlinearEulerSolution(FESystemUInt dim, FESystem::Mesh::ElementType elem_type, FESystemUInt n_elem_nodes, const
                            FESystem::Mesh::MeshBase& mesh, const FESystem::Base::DegreeOfFreedomMap& dof_map,
                            const FESystem::Numerics::SparsityPattern& nonbc_sparsity_pattern,
                            const std::map<FESystemUInt, FESystemUInt>& old_to_new_id_map, const std::vector<FESystemUInt>& nonbc_dofs)
{
    const std::vector<FESystem::Mesh::Node*>& nodes = mesh.getNodes();
    
    FESystem::Numerics::SparseMatrix<FESystemDouble> stiff_mat, mass;
    FESystem::Numerics::LocalVector<FESystemDouble> residual, sol, vel;
    
    residual.resize(dof_map.getNDofs()); sol.resize(dof_map.getNDofs()); vel.resize(dof_map.getNDofs());
    stiff_mat.resize(dof_map.getSparsityPattern()); mass.resize(dof_map.getSparsityPattern());

    
    FESystem::LinearSolvers::LUFactorizationLinearSolver<FESystemDouble> linear_solver;
    FESystem::NonlinearSolvers::NewtonIterationNonlinearSolver<FESystemDouble> nonlinear_solver;
    
    nonlinear_solver.initialize(stiff_mat, linear_solver);
    nonlinear_solver.setConvergenceLimits(20, 1.0e-10);
    
    // for output
    FESystem::OutputProcessor::VtkOutputProcessor output;
    std::fstream output_file;
    std::vector<FESystemUInt> vars(5); vars[0]=0; vars[1]=1; vars[2]=2; vars[3]=3; vars[4]=4; // write all solutions
    
    
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
                    vec.setVal(nodes[i]->getDegreeOfFreedomUnit(0).global_dof_id[0], rho); // rho
                    vec.setVal(nodes[i]->getDegreeOfFreedomUnit(1).global_dof_id[0], rho * u1); // rho u
                    vec.setVal(nodes[i]->getDegreeOfFreedomUnit(3).global_dof_id[0], rho * (temp*cv + u1*u1*0.5)); // rho * (cv * T + u^2/2)
                }
            }
                break;
                
            case FESystem::NonlinearSolvers::SOLUTION_CONVERGED:
            {
                // write the solution
                std::stringstream oss;
                oss << "sol_" << nonlinear_solver.getCurrentIterationNumber() << ".vtk";
                output_file.open(oss.str().c_str(),std::fstream::out);
                output.writeMesh(output_file, mesh, dof_map);
                output.writeSolution(output_file, "Sol", mesh, dof_map, vars, nonlinear_solver.getCurrentSolution());
                output_file.close();
            }
                break;
                
            case FESystem::NonlinearSolvers::EVALUATE_RESIDUAL:
            case FESystem::NonlinearSolvers::EVALUATE_JACOBIAN:
            case FESystem::NonlinearSolvers::EVALUATE_RESIDUAL_AND_JACOBIAN:
            {
                calculateEulerQuantities(elem_type, n_elem_nodes, dof_map, mesh, nonlinear_solver.getCurrentSolution(), vel,
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
        for (FESystemUInt j=0; j<3; j++)
            std::cout << std::setw(15) << sol.getVal(nodes[i]->getDegreeOfFreedomUnit(j).global_dof_id[0]);
        std::cout << std::endl;
    }
}




void modifyMeshForCase(FESystem::Mesh::MeshBase& mesh)
{
    {
        std::vector<FESystem::Mesh::Node*>& nodes = mesh.getNodes();
        std::vector<FESystem::Mesh::Node*>::iterator it=nodes.begin(), end=nodes.end();
        
        FESystemDouble x_val, y_val;
        
        for ( ; it!=end; it++)
        {
            if (((*it)->getVal(0) >= x0) && ((*it)->getVal(0) <= x1))
            {
                x_val = (*it)->getVal(0);
                y_val = (*it)->getVal(1);
                
                y_val += thickness*(1.0-y_val/y_length)*sin(PI_VAL*(x_val-x0)/chord);
                
                (*it)->setVal(1, y_val);
            }
        }
    }

    {
        std::vector<FESystem::Mesh::ElemBase*> elems = mesh.getElements();
        std::vector<FESystem::Mesh::ElemBase*>::iterator it=elems.begin(), end=elems.end();
        
        for ( ; it!=end; it++)
            (*it)->updateAfterMeshDeformation();
    }
}





int euler_analysis_driver(int argc, char * const argv[])
{
    FESystem::Mesh::ElementType elem_type;
    
    // create the mesh object
    FESystem::Mesh::MeshBase mesh;
    
    // create a nx x ny grid of nodes
    FESystemUInt n_elem_nodes, n_elem_dofs;
    FESystem::Geometry::Point origin(3);
    
    elem_type = FESystem::Mesh::QUAD4;
    createPlaneMesh(elem_type, mesh, origin, nx, ny, x_length, y_length, n_elem_nodes, CROSS, true);
    
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
    FESystem::Numerics::SparseMatrix<FESystemDouble> global_stiffness_mat, global_mass_mat;
    FESystem::Numerics::LocalVector<FESystemDouble> rhs, sol;
    std::vector<FESystemUInt> elem_dof_indices;
    
    rhs.resize(dof_map.getNDofs()); sol.resize(dof_map.getNDofs());
    global_mass_mat.resize(dof_map.getSparsityPattern());
    global_stiffness_mat.resize(dof_map.getSparsityPattern());
    
    // apply boundary condition and place a load on the last dof
    std::set<FESystemUInt> bc_dofs;
    for (FESystemUInt i=0; i<nodes.size(); i++)
    {
        if ((nodes[i]->getVal(0) == 0.0)) // left boundary nodes
        {
            bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(0).global_dof_id[0]); // rho
            bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(1).global_dof_id[0]); // rho u1
            bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(2).global_dof_id[0]); // rho u2
            bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(3).global_dof_id[0]); // rho e
            const std::set<FESystem::Mesh::ElemBase*>& e_set = nodes[i]->getElementConnectivitySet();
            std::set<FESystem::Mesh::ElemBase*>::const_iterator it = e_set.begin(), end = e_set.end();
            for ( ; it != end; it++)
                if (!(*it)->checkForTag(0))
                    (*it)->setTag(0);
        }
        
        if ((nodes[i]->getVal(0) == x_length)) // right boundary nodes
        {
            const std::set<FESystem::Mesh::ElemBase*>& e_set = nodes[i]->getElementConnectivitySet();
            std::set<FESystem::Mesh::ElemBase*>::const_iterator it = e_set.begin(), end = e_set.end();
            for ( ; it != end; it++)
                if (!(*it)->checkForTag(1))
                    (*it)->setTag(1);
        }

        if ((nodes[i]->getVal(1) == 0.0)) // lower boundary nodes
        {
            const std::set<FESystem::Mesh::ElemBase*>& e_set = nodes[i]->getElementConnectivitySet();
            std::set<FESystem::Mesh::ElemBase*>::const_iterator it = e_set.begin(), end = e_set.end();
            for ( ; it != end; it++)
                if (!(*it)->checkForTag(2))
                    (*it)->setTag(2);
        }
        
        
        if ((nodes[i]->getVal(1) == y_length)) // upper boundary nodes
        {
            const std::set<FESystem::Mesh::ElemBase*>& e_set = nodes[i]->getElementConnectivitySet();
            std::set<FESystem::Mesh::ElemBase*>::const_iterator it = e_set.begin(), end = e_set.end();
            for ( ; it != end; it++)
                if (!(*it)->checkForTag(3))
                    (*it)->setTag(3);
        }
    }
    
    // now create the vector of ids that do not have bcs
    std::vector<FESystemUInt> nonbc_dofs;
    for (FESystemUInt i=0; i<dof_map.getNDofs(); i++)
        if (!bc_dofs.count(i))
            nonbc_dofs.push_back(i);

    
    // modify mesh after application of boundary condition
    modifyMeshForCase(mesh);

    
    // prepare a map of the old to new ID
    std::vector<FESystemUInt>::const_iterator dof_it=nonbc_dofs.begin(), dof_end=nonbc_dofs.end();
    std::map<FESystemUInt, FESystemUInt> old_to_new_id_map;
    FESystemUInt n=0;
    for ( ; dof_it!=dof_end; dof_it++)
        old_to_new_id_map.insert(std::map<FESystemUInt, FESystemUInt>::value_type(*dof_it, n++));
    
    FESystem::Numerics::SparsityPattern nonbc_sparsity_pattern;
    dof_map.getSparsityPattern().initSparsityPatternForNonConstrainedDOFs(nonbc_dofs, old_to_new_id_map, nonbc_sparsity_pattern);
    
    //nonlinearEulerSolution(dim, elem_type, n_elem_nodes, mesh, dof_map, nonbc_sparsity_pattern, old_to_new_id_map, nonbc_dofs);
    transientEulerAnalysis(dim, elem_type, n_elem_nodes, mesh, dof_map, nonbc_sparsity_pattern, old_to_new_id_map, nonbc_dofs);
    
    //exit(1);
    
    return 0;
}





