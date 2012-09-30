//
//  SolutionDrivers.cpp
//  FESystem
//
//  Created by Manav Bhatia on 9/26/12.
//
//


// FESystem include
#include "TestingIncludes.h"



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
                       void (*calculateStructuralMatrices)(FESystemBoolean if_nonlinear, FESystem::Mesh::ElementType elem_type, FESystemUInt n_elem_nodes,
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
                nonlinear_solver.getCurrentSolution().zero();
                break;
                
            case FESystem::NonlinearSolvers::EVALUATE_RESIDUAL:
            {
                // zero the solution vectors
                rhs.zero(); sol.zero();
                // get the latest solution vector
                sol.addVal(nonbc_dofs, nonlinear_solver.getCurrentSolution());
                calculateStructuralMatrices(true, elem_type, n_elem_nodes, dof_map, mesh, sol, internal_force, rhs, stiff_mat, dummy);
                internal_force.getSubVectorValsFromIndices(nonbc_dofs, nonlinear_solver.getResidualVector());
                rhs.getSubVectorValsFromIndices(nonbc_dofs, reduced_load_vec);
                nonlinear_solver.getResidualVector().add(-1.0, reduced_load_vec);
                
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
                rhs.zero(); sol.zero();
                // get the latest solution vector
                sol.addVal(nonbc_dofs, nonlinear_solver.getCurrentSolution());
                calculateStructuralMatrices(true, elem_type, n_elem_nodes, dof_map, mesh, sol, internal_force, rhs, stiff_mat, dummy);
                stiff_mat.getSubMatrixValsFromRowAndColumnIndices(nonbc_dofs, nonbc_dofs, old_to_new_id_map, reduced_stiff_mat);
            }
                break;
                
            case FESystem::NonlinearSolvers::EVALUATE_RESIDUAL_AND_JACOBIAN:
            {
                // zero the solution vectors
                rhs.zero(); sol.zero();
                // get the latest solution vector
                sol.addVal(nonbc_dofs, nonlinear_solver.getCurrentSolution());
                calculateStructuralMatrices(true, elem_type, n_elem_nodes, dof_map, mesh, sol, internal_force, rhs, stiff_mat, dummy);
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
        if (call_back == FESystem::NonlinearSolvers::MAXIMUM_ITERATIONS_REACHED)
        {
            std::cout << "Maximum iterations reached" << std::endl;
            break;
        }
    }
    
    // copy the last solution and output
    reduced_sol_vec.copyVector(nonlinear_solver.getCurrentSolution());
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
    eig_vec.copyMatrix(eigen_solver.getRightEigenVectorMatrix());
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
    FESystemUInt n_cycles = 10;
    FESystemDouble final_t=1.0/(sqrt(eig_vals.getVal(id))/2.0/3.141)*n_cycles, time_step=final_t/n_cycles*1.0e-3;
    
    // initialize the solver
    FESystem::TransientSolvers::NewmarkTransientSolver<FESystemDouble> transient_solver;
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





int test_ode_integration(int argc, char * const argv[])
{
    
    FESystem::Plotting::PLPlot<FESystemDouble> plot(FESystem::Plotting::REAL_AXIS, FESystem::Plotting::REAL_AXIS);
    FESystemDouble omega2=250.0, final_t=1.0/(sqrt(omega2)/2.0/3.141)*10, time_step=final_t*1.0e-4;
    
    // initialize the solver
    FESystem::TransientSolvers::NewmarkTransientSolver<FESystemDouble> transient_solver;
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
    
    plot.plotData2D(x_vals, y_vals);
    
    
    return 0;
}




