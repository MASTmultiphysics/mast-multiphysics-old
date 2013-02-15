//
//  PistonTheoryTest.cpp
//  FESystem
//
//  Created by Manav Bhatia on 9/26/12.
//
//



// FESystem include
#include "TestingIncludes.h"



void transientFlutterSolution(FESystemUInt dim, FESystem::Mesh::ElementType elem_type, FESystemUInt n_elem_nodes,
                              const FESystem::Mesh::MeshBase& mesh, const FESystem::DegreeOfFreedom::DegreeOfFreedomMap& dof_map, const FESystem::Numerics::SparsityPattern& nonbc_sparsity_pattern,
                              const std::map<FESystemUInt, FESystemUInt>& old_to_new_id_map, const std::vector<FESystemUInt>& nonbc_dofs, const FESystem::Numerics::MatrixBase<FESystemDouble>& stiff_mat,
                              const FESystem::Numerics::VectorBase<FESystemDouble>& mass_vec, const FESystemUInt n_modes, std::vector<FESystemUInt>& sorted_ids, FESystem::Numerics::VectorBase<FESystemDouble>& eig_vals,
                              FESystem::Numerics::MatrixBase<FESystemDouble>& eig_vec,
                              void (*calculatePistonTheoryMatrices)(FESystemDouble mach, FESystemDouble rho, FESystemDouble gamma, FESystemDouble a_inf, FESystemDouble u_inf,
                                                                    FESystemBoolean if_nonlinear, FESystem::Mesh::ElementType elem_type, FESystemUInt n_elem_nodes, const FESystem::DegreeOfFreedom::DegreeOfFreedomMap& dof_map,
                                                                    const std::map<FESystemUInt, FESystemUInt>& old_to_new_id_map, const std::vector<FESystemUInt>& nonbc_dofs, const FESystem::Mesh::MeshBase& mesh,
                                                                    FESystem::Numerics::VectorBase<FESystemDouble>& global_sol, FESystem::Numerics::VectorBase<FESystemDouble>& global_vel,
                                                                    FESystem::Numerics::VectorBase<FESystemDouble>& force, FESystem::Numerics::VectorBase<FESystemDouble>& generalized_force,
                                                                    FESystem::Numerics::MatrixBase<FESystemDouble>& aero_stiffness_mat, FESystem::Numerics::MatrixBase<FESystemDouble>& reduced_stiffness_mat, FESystem::Numerics::MatrixBase<FESystemDouble>& generalized_stiffness_mat,
                                                                    FESystem::Numerics::MatrixBase<FESystemDouble>& aero_damp_mat, FESystem::Numerics::MatrixBase<FESystemDouble>& reduced_damp_mat, FESystem::Numerics::MatrixBase<FESystemDouble>& generalized_damp_mat,
                                                                    const FESystemUInt n_modes, std::vector<FESystemUInt>& sorted_ids, FESystem::Numerics::VectorBase<FESystemDouble>& eig_vals, FESystem::Numerics::MatrixBase<FESystemDouble>& eig_vec))
{
    FESystem::Numerics::SparseMatrix<FESystemDouble> aero_stiff, aero_damp, reduced_aero_stiff, reduced_aero_damp;
    FESystem::Numerics::DenseMatrix<FESystemDouble> aero_force, generalized_aero_stiff, generalized_aero_damp, elem_mat, eig_mat;
    FESystem::Numerics::LocalVector<FESystemDouble> dummy_vec, real_val, imag_val;
    real_val.resize(n_modes); imag_val.resize(n_modes);
    aero_damp.resize(dof_map.getSparsityPattern()); aero_stiff.resize(dof_map.getSparsityPattern());
    reduced_aero_damp.resize(nonbc_sparsity_pattern); reduced_aero_stiff.resize(nonbc_sparsity_pattern);
    generalized_aero_damp.resize(n_modes, n_modes); generalized_aero_stiff.resize(n_modes, n_modes);
    aero_force.resize(dof_map.getNDofs(), dof_map.getNDofs()); elem_mat.resize(n_elem_nodes, n_elem_nodes); eig_mat.resize(2*n_modes, 2*n_modes);
    
    FESystem::EigenSolvers::LapackLinearEigenSolver<FESystemDouble> eigen_solver;
    
    FESystemDouble mass_prop_damp_coeff = 0.0, stiff_prop_damp_coeff = 0.0, q_dyn = 0.0, mach = 2.00, rho = 0.05, gamma = 1.4, a_inf = 0, u_inf = 29;
    
    //FESystem::Plotting::PLPlot<FESystemDouble> plot(FESystem::Plotting::REAL_AXIS, FESystem::Plotting::REAL_AXIS);
    FESystemComplexDouble val;
    
    // set the values in the first order matrix
    q_dyn = 0.5 * rho * u_inf * u_inf;
    
    // calculate the aeroelastic quantities
    calculatePistonTheoryMatrices(mach, rho, gamma, a_inf, u_inf, false, elem_type, n_elem_nodes, dof_map, old_to_new_id_map, nonbc_dofs, mesh, dummy_vec, dummy_vec, dummy_vec, dummy_vec,
                                  aero_stiff, reduced_aero_stiff, generalized_aero_stiff, aero_damp, reduced_aero_damp, generalized_aero_damp,
                                  n_modes, sorted_ids, eig_vals, eig_vec);
    
    // add to the eigenvalue matrix: scale with the dynamic pressure
    eig_mat.addSubMatrixVals(n_modes, 2*n_modes-1, 0, n_modes-1, 0, n_modes-1, 0, n_modes-1, q_dyn, generalized_aero_stiff);
    eig_mat.addSubMatrixVals(n_modes, 2*n_modes-1, n_modes, 2*n_modes-1, 0, n_modes-1, 0, n_modes-1, q_dyn, generalized_aero_damp);
    
    
    
    
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
    reduced_stiff_mat.scale(-1.0);
    reduced_stiff_mat.add(q_dyn, reduced_aero_stiff);
    reduced_aero_stiff.scale(q_dyn);
    for (FESystemUInt i=0; i<reduced_mass_vec.getSize(); i++)
    {
        reduced_stiff_mat.scaleRow(i, 1.0/reduced_mass_vec.getVal(i));
        reduced_aero_damp.scaleRow(i, 1.0/reduced_mass_vec.getVal(i));
    }
    
    id = sorted_ids[0];
    eig_vec.getColumnVals(id, 0, nonbc_dofs.size()-1, reduced_sol_vec);
    FESystemUInt n_cycles = 300;
    FESystemDouble final_t=1.0/(sqrt(eig_vals.getVal(id))/2.0/3.141)*n_cycles, time_step=final_t/n_cycles*1.0e-3;
    
    // initialize the solver
    FESystem::TransientSolvers::NewmarkTransientSolver<FESystemDouble> transient_solver;
    FESystem::Numerics::SparsityPattern ode_sparsity;
    FESystem::Numerics::SparseMatrix<FESystemDouble> ode_jac;
    std::vector<FESystemBoolean> ode_order_include(2); ode_order_include[0] = true; ode_order_include[1]=true;
    std::vector<FESystemDouble> int_constants(2); int_constants[0]=1.0/4.0; int_constants[1]=1.0/2.0;
    transient_solver.initialize(2, nonbc_dofs.size(), int_constants);
    transient_solver.setActiveJacobianTerm(ode_order_include);
    transient_solver.setMassMatrix(true);
    
    transient_solver.initializeStateVector(rhs); //  initialize the vector and apply the initial condition
    transient_solver.updateVectorValuesForDerivativeOrder(0, reduced_sol_vec, rhs);
    transient_solver.setInitialTimeData(0, time_step, rhs);
    
    transient_solver.initializeMatrixSparsityPatterForSystem(nonbc_sparsity_pattern, ode_sparsity);
    ode_jac.resize(ode_sparsity);
    transient_solver.setJacobianMatrix(ode_jac);
    transient_solver.setLinearSolver(linear_solver, true);
    transient_solver.updateJacobianValuesForDerivativeOrder(2, 0, reduced_stiff_mat, transient_solver.getCurrentJacobianMatrix());
    transient_solver.updateJacobianValuesForDerivativeOrder(2, 1, reduced_aero_damp, transient_solver.getCurrentJacobianMatrix());
    
    
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
            {
                // the Jacobian is not updated since it is constant with respect to time
                transient_solver.extractVectorValuesForDerivativeOrder(0, transient_solver.getCurrentStateVector(), reduced_sol_vec); // get the current X
                reduced_stiff_mat.rightVectorMultiply(reduced_sol_vec, reduced_load_vec);
                transient_solver.extractVectorValuesForDerivativeOrder(1, transient_solver.getCurrentStateVector(), reduced_sol_vec); // get the current X
                reduced_aero_damp.rightVectorMultiply(reduced_sol_vec, reduced_mass_vec);
                reduced_load_vec.add(1.0, reduced_mass_vec);
                
                transient_solver.updateVectorValuesForDerivativeOrder(1, reduced_load_vec, transient_solver.getCurrentVelocityFunctionVector()); // set the acceleration
                transient_solver.copyDerivativeValuesFromStateToVelocityVector(transient_solver.getCurrentStateVector(), transient_solver.getCurrentVelocityFunctionVector());
            }
                break;
                
            default:
                FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
                break;
        }
    }
}




void flutterSolution(FESystemUInt dim, FESystem::Mesh::ElementType elem_type, FESystemUInt n_elem_nodes,
                     const FESystem::Mesh::MeshBase& mesh, const FESystem::DegreeOfFreedom::DegreeOfFreedomMap& dof_map, const FESystem::Numerics::SparsityPattern& nonbc_sparsity_pattern,
                     const std::map<FESystemUInt, FESystemUInt>& old_to_new_id_map, const std::vector<FESystemUInt>& nonbc_dofs, const FESystem::Numerics::MatrixBase<FESystemDouble>& stiff_mat,
                     const FESystem::Numerics::VectorBase<FESystemDouble>& mass_vec, const FESystemUInt n_modes, std::vector<FESystemUInt>& sorted_ids, FESystem::Numerics::VectorBase<FESystemDouble>& eig_vals,
                     FESystem::Numerics::MatrixBase<FESystemDouble>& eig_vec,
                     void (*calculatePistonTheoryMatrices)(FESystemDouble mach, FESystemDouble rho, FESystemDouble gamma, FESystemDouble a_inf, FESystemDouble u_inf,
                                                           FESystemBoolean if_nonlinear, FESystem::Mesh::ElementType elem_type, FESystemUInt n_elem_nodes, const FESystem::DegreeOfFreedom::DegreeOfFreedomMap& dof_map,
                                                           const std::map<FESystemUInt, FESystemUInt>& old_to_new_id_map, const std::vector<FESystemUInt>& nonbc_dofs, const FESystem::Mesh::MeshBase& mesh,
                                                           FESystem::Numerics::VectorBase<FESystemDouble>& global_sol, FESystem::Numerics::VectorBase<FESystemDouble>& global_vel,
                                                           FESystem::Numerics::VectorBase<FESystemDouble>& force, FESystem::Numerics::VectorBase<FESystemDouble>& generalized_force,
                                                           FESystem::Numerics::MatrixBase<FESystemDouble>& aero_stiffness_mat, FESystem::Numerics::MatrixBase<FESystemDouble>& reduced_stiffness_mat, FESystem::Numerics::MatrixBase<FESystemDouble>& generalized_stiffness_mat,
                                                           FESystem::Numerics::MatrixBase<FESystemDouble>& aero_damp_mat, FESystem::Numerics::MatrixBase<FESystemDouble>& reduced_damp_mat, FESystem::Numerics::MatrixBase<FESystemDouble>& generalized_damp_mat,
                                                           const FESystemUInt n_modes, std::vector<FESystemUInt>& sorted_ids, FESystem::Numerics::VectorBase<FESystemDouble>& eig_vals, FESystem::Numerics::MatrixBase<FESystemDouble>& eig_vec))
{
    FESystem::Numerics::SparseMatrix<FESystemDouble> aero_stiff, aero_damp, reduced_stiff, reduced_damp;
    FESystem::Numerics::DenseMatrix<FESystemDouble> aero_force, generalized_aero_stiff, generalized_aero_damp, elem_mat, eig_mat;
    FESystem::Numerics::LocalVector<FESystemDouble> dummy_vec, real_val, imag_val;
    real_val.resize(n_modes); imag_val.resize(n_modes);
    aero_damp.resize(dof_map.getSparsityPattern()); aero_stiff.resize(dof_map.getSparsityPattern());
    reduced_damp.resize(nonbc_sparsity_pattern); reduced_stiff.resize(nonbc_sparsity_pattern);
    generalized_aero_damp.resize(n_modes, n_modes); generalized_aero_stiff.resize(n_modes, n_modes);
    aero_force.resize(dof_map.getNDofs(), dof_map.getNDofs()); elem_mat.resize(n_elem_nodes, n_elem_nodes); eig_mat.resize(2*n_modes, 2*n_modes);
    
    FESystem::EigenSolvers::LapackLinearEigenSolver<FESystemDouble> eigen_solver;
    
    FESystemDouble mass_prop_damp_coeff = 0.0, stiff_prop_damp_coeff = 0.0, q_dyn = 0.0, mach = 2.00, rho = 0.05, gamma = 1.4, a_inf = 0, u_inf = 1;
    
    //FESystem::Plotting::PLPlot<FESystemDouble> plot(FESystem::Plotting::REAL_AXIS, FESystem::Plotting::REAL_AXIS);
    FESystemComplexDouble val;
    
    // set the values in the first order matrix
    for (FESystemUInt i=0; i<100; i++)
    {
        q_dyn = 0.5 * rho * u_inf * u_inf;
        
        eig_mat.zero();
        for (FESystemUInt i=0; i<n_modes; i++)
        {
            eig_mat.setVal(i, n_modes+i, 1.0); // xdot
            eig_mat.setVal(n_modes+i, i, -eig_vals.getVal(sorted_ids[i])); // stiffness: (-omega^2)
            eig_mat.setVal(n_modes+i, n_modes+i, -eig_vals.getVal(sorted_ids[i])*stiff_prop_damp_coeff - mass_prop_damp_coeff); // stiffness: (-omega^2)
        }
        
        // calculate the aeroelastic quantities
        calculatePistonTheoryMatrices(mach, rho, gamma, a_inf, u_inf, false, elem_type, n_elem_nodes, dof_map, old_to_new_id_map, nonbc_dofs, mesh, dummy_vec, dummy_vec, dummy_vec, dummy_vec,
                                      aero_stiff, reduced_stiff, generalized_aero_stiff, aero_damp, reduced_damp, generalized_aero_damp,
                                      n_modes, sorted_ids, eig_vals, eig_vec);
        
        // add to the eigenvalue matrix: scale with the dynamic pressure
        eig_mat.addSubMatrixVals(n_modes, 2*n_modes-1, 0, n_modes-1, 0, n_modes-1, 0, n_modes-1, q_dyn, generalized_aero_stiff);
        eig_mat.addSubMatrixVals(n_modes, 2*n_modes-1, n_modes, 2*n_modes-1, 0, n_modes-1, 0, n_modes-1, q_dyn, generalized_aero_damp);
        
        //eig_mat.write(std::cout);
        std::cout << "Iter: " << i << ":  v = " << u_inf << ":  q_dyn = " << q_dyn << std::endl;
        
        //eig_mat.write(std::cout);
        eigen_solver.setEigenProblemType(FESystem::EigenSolvers::NONHERMITIAN);
        eigen_solver.setMatrix(&eig_mat);
        eigen_solver.solve();
        eigen_solver.getComplexEigenValues().write(std::cout);
        //eigen_solver.getRightComplexEigenVectorMatrix().write(std::cout);
        //eigen_solver.getLeftComplexEigenVectorMatrix().write(std::cout);
        
        for (FESystemUInt i=0; i<n_modes; i++)
        {
            val = eigen_solver.getComplexEigenValues().getVal(i);
            real_val.setVal(i, std::real(val));
            imag_val.setVal(i, std::imag(val));
        }
        //plot.plotData2D(real_val, imag_val);
        u_inf += 1;
    }
    
}




void calculatePlatePistonTheoryMatrices(FESystemDouble mach, FESystemDouble rho, FESystemDouble gamma, FESystemDouble a_inf, FESystemDouble u_inf,
                                        FESystemBoolean if_nonlinear, FESystem::Mesh::ElementType elem_type, FESystemUInt n_elem_nodes, const FESystem::DegreeOfFreedom::DegreeOfFreedomMap& dof_map,
                                        const std::map<FESystemUInt, FESystemUInt>& old_to_new_id_map, const std::vector<FESystemUInt>& nonbc_dofs, const FESystem::Mesh::MeshBase& mesh,
                                        FESystem::Numerics::VectorBase<FESystemDouble>& global_sol, FESystem::Numerics::VectorBase<FESystemDouble>& global_vel,
                                        FESystem::Numerics::VectorBase<FESystemDouble>& force, FESystem::Numerics::VectorBase<FESystemDouble>& generalized_force,
                                        FESystem::Numerics::MatrixBase<FESystemDouble>& aero_stiffness_mat, FESystem::Numerics::MatrixBase<FESystemDouble>& reduced_stiffness_mat, FESystem::Numerics::MatrixBase<FESystemDouble>& generalized_stiffness_mat,
                                        FESystem::Numerics::MatrixBase<FESystemDouble>& aero_damp_mat, FESystem::Numerics::MatrixBase<FESystemDouble>& reduced_damp_mat, FESystem::Numerics::MatrixBase<FESystemDouble>& generalized_damp_mat,
                                        const FESystemUInt n_modes, std::vector<FESystemUInt>& sorted_ids, FESystem::Numerics::VectorBase<FESystemDouble>& eig_vals, FESystem::Numerics::MatrixBase<FESystemDouble>& eig_vec)
{
    FESystemUInt n_plate_dofs, n_elem_dofs, order = 1;
    
    n_plate_dofs = n_elem_nodes;
    n_elem_dofs = 6*n_elem_nodes;
    
    
    FESystem::Numerics::DenseMatrix<FESystemDouble> elem_mat, dfdw, dfdwdot;
    FESystem::Numerics::LocalVector<FESystemDouble> elem_vec, beam_elem_vec, elem_sol, elem_vel, elem_local_sol, elem_local_vel;
    std::vector<FESystemUInt> elem_dof_indices;
    elem_mat.resize(n_elem_dofs, n_elem_dofs);
    elem_vec.resize(n_elem_dofs); elem_sol.resize(n_elem_dofs); elem_vel.resize(n_elem_dofs); elem_local_sol.resize(n_plate_dofs); elem_local_vel.resize(n_plate_dofs);
    dfdw.resize(n_plate_dofs, n_plate_dofs); dfdwdot.resize(n_plate_dofs, n_plate_dofs);
    
    // prepare the quadrature rule and FE for the
    FESystem::Quadrature::TrapezoidQuadrature q_rule;
    FESystem::FiniteElement::FELagrange fe;
    FESystem::Structures::PistonTheory2D piston_elem;
    
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
    
    force.zero(); aero_stiffness_mat.zero(); aero_damp_mat.zero();
    
    const std::vector<FESystem::Mesh::ElemBase*>& elems = mesh.getElements();
    
    for (FESystemUInt i=0; i<elems.size(); i++)
    {
        fe.clear();
        fe.reinit(*(elems[i]));
        elem_mat.zero();
        elem_vec.zero();
        elem_sol.zero();
        
        if (if_nonlinear)
        {
            dof_map.getFromGlobalVector(*(elems[i]), global_sol, elem_sol);
            dof_map.getFromGlobalVector(*(elems[i]), global_vel, elem_vel);
        }
        
        // initialize element
        piston_elem.clear();
        piston_elem.initialize(*(elems[i]), fe, q_rule, order, mach, a_inf, u_inf, gamma);
        
        // calculate quantities
        piston_elem.calculateTangentMatrix(elem_local_sol, elem_local_vel, dfdw, dfdwdot);
        piston_elem.transformMatrixToGlobalSystem(dfdw, elem_mat);
        dof_map.addToGlobalMatrix(*(elems[i]), elem_mat, aero_stiffness_mat);
        piston_elem.transformMatrixToGlobalSystem(dfdwdot, elem_mat);
        dof_map.addToGlobalMatrix(*(elems[i]), elem_mat, aero_damp_mat);
    }
    
    // extract the constrained matrices from this
    aero_stiffness_mat.getSubMatrixValsFromRowAndColumnIndices(nonbc_dofs, nonbc_dofs, old_to_new_id_map, reduced_stiffness_mat);
    aero_damp_mat.getSubMatrixValsFromRowAndColumnIndices(nonbc_dofs, nonbc_dofs, old_to_new_id_map, reduced_damp_mat);
    
    // now calculate the generalized matrices
    const std::pair<FESystemUInt, FESystemUInt> s = eig_vec.getSize();
    FESystem::Numerics::LocalVector<FESystemDouble> tmp_vec1, tmp_vec2, tmp_vec3;
    tmp_vec1.resize(s.first); tmp_vec2.resize(s.first); tmp_vec3.resize(s.first);
    
    for (FESystemUInt j=0; j<n_modes; j++)
    {
        tmp_vec1.zero(); tmp_vec2.zero(); tmp_vec3.zero();
        eig_vec.getColumnVals(sorted_ids[j], 0, s.first-1, tmp_vec1);
        reduced_stiffness_mat.rightVectorMultiply(tmp_vec1, tmp_vec2);
        reduced_damp_mat.rightVectorMultiply(tmp_vec1, tmp_vec3);
        for (FESystemUInt i=0; i<n_modes; i++)
        {
            tmp_vec1.zero();
            eig_vec.getColumnVals(sorted_ids[i], 0, s.first-1, tmp_vec1);
            generalized_stiffness_mat.setVal(i, j, tmp_vec1.dotProduct(tmp_vec2));
            generalized_damp_mat.setVal(i, j, tmp_vec1.dotProduct(tmp_vec3));
        }
    }
}





void calculateBeamPistonTheoryMatrices(FESystemDouble mach, FESystemDouble rho, FESystemDouble gamma, FESystemDouble a_inf, FESystemDouble u_inf,
                                       FESystemBoolean if_nonlinear, FESystem::Mesh::ElementType elem_type, FESystemUInt n_elem_nodes, const FESystem::DegreeOfFreedom::DegreeOfFreedomMap& dof_map,
                                       const std::map<FESystemUInt, FESystemUInt>& old_to_new_id_map, const std::vector<FESystemUInt>& nonbc_dofs, const FESystem::Mesh::MeshBase& mesh,
                                       FESystem::Numerics::VectorBase<FESystemDouble>& global_sol, FESystem::Numerics::VectorBase<FESystemDouble>& global_vel,
                                       FESystem::Numerics::VectorBase<FESystemDouble>& force, FESystem::Numerics::VectorBase<FESystemDouble>& generalized_force,
                                       FESystem::Numerics::MatrixBase<FESystemDouble>& aero_stiffness_mat, FESystem::Numerics::MatrixBase<FESystemDouble>& reduced_stiffness_mat, FESystem::Numerics::MatrixBase<FESystemDouble>& generalized_stiffness_mat,
                                       FESystem::Numerics::MatrixBase<FESystemDouble>& aero_damp_mat, FESystem::Numerics::MatrixBase<FESystemDouble>& reduced_damp_mat, FESystem::Numerics::MatrixBase<FESystemDouble>& generalized_damp_mat,
                                       const FESystemUInt n_modes, std::vector<FESystemUInt>& sorted_ids, FESystem::Numerics::VectorBase<FESystemDouble>& eig_vals, FESystem::Numerics::MatrixBase<FESystemDouble>& eig_vec)
{
    FESystemUInt n_beam_dofs, n_elem_dofs, order = 1;
    
    n_beam_dofs = n_elem_nodes;
    n_elem_dofs = 6*n_elem_nodes;
    
    
    FESystem::Numerics::DenseMatrix<FESystemDouble> elem_mat, dfdw, dfdwdot;
    FESystem::Numerics::LocalVector<FESystemDouble> elem_vec, beam_elem_vec, elem_sol, elem_vel, elem_local_sol, elem_local_vel;
    std::vector<FESystemUInt> elem_dof_indices;
    elem_mat.resize(n_elem_dofs, n_elem_dofs);
    elem_vec.resize(n_elem_dofs); elem_sol.resize(n_elem_dofs); elem_vel.resize(n_elem_dofs); elem_local_sol.resize(n_beam_dofs); elem_local_vel.resize(n_beam_dofs);
    dfdw.resize(n_beam_dofs, n_beam_dofs); dfdwdot.resize(n_beam_dofs, n_beam_dofs);
    
    // prepare the quadrature rule and FE for the
    FESystem::Quadrature::TrapezoidQuadrature q_rule;
    FESystem::FiniteElement::FELagrange fe;
    FESystem::Structures::PistonTheory1D piston_elem;
    
    switch (elem_type)
    {
        case FESystem::Mesh::EDGE2:
            q_rule.init(1, 3);
            break;
            
        case FESystem::Mesh::EDGE3:
            q_rule.init(1, 5);
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
    }
    
    force.zero(); aero_stiffness_mat.zero(); aero_damp_mat.zero();
    
    const std::vector<FESystem::Mesh::ElemBase*>& elems = mesh.getElements();
    
    for (FESystemUInt i=0; i<elems.size(); i++)
    {
        fe.clear();
        fe.reinit(*(elems[i]));
        elem_mat.zero();
        elem_vec.zero();
        elem_sol.zero();
        
        if (if_nonlinear)
        {
            dof_map.getFromGlobalVector(*(elems[i]), global_sol, elem_sol);
            dof_map.getFromGlobalVector(*(elems[i]), global_vel, elem_vel);
        }
        
        // initialize element
        piston_elem.clear();
        piston_elem.initialize(*(elems[i]), fe, q_rule, order, mach, a_inf, u_inf, gamma);
        
        // calculate quantities
        piston_elem.calculateTangentMatrix(elem_local_sol, elem_local_vel, dfdw, dfdwdot);
        piston_elem.transformMatrixToGlobalSystem(dfdw, elem_mat);
        dof_map.addToGlobalMatrix(*(elems[i]), elem_mat, aero_stiffness_mat);
        piston_elem.transformMatrixToGlobalSystem(dfdwdot, elem_mat);
        dof_map.addToGlobalMatrix(*(elems[i]), elem_mat, aero_damp_mat);
    }
    
    // extract the constrained matrices from this
    aero_stiffness_mat.getSubMatrixValsFromRowAndColumnIndices(nonbc_dofs, nonbc_dofs, old_to_new_id_map, reduced_stiffness_mat);
    aero_damp_mat.getSubMatrixValsFromRowAndColumnIndices(nonbc_dofs, nonbc_dofs, old_to_new_id_map, reduced_damp_mat);
    
    // now calculate the generalized matrices
    const std::pair<FESystemUInt, FESystemUInt> s = eig_vec.getSize();
    FESystem::Numerics::LocalVector<FESystemDouble> tmp_vec1, tmp_vec2, tmp_vec3;
    tmp_vec1.resize(s.first); tmp_vec2.resize(s.first); tmp_vec3.resize(s.first);
    
    for (FESystemUInt j=0; j<n_modes; j++)
    {
        tmp_vec1.zero(); tmp_vec2.zero(); tmp_vec3.zero();
        eig_vec.getColumnVals(sorted_ids[j], 0, s.first-1, tmp_vec1);
        reduced_stiffness_mat.rightVectorMultiply(tmp_vec1, tmp_vec2);
        reduced_damp_mat.rightVectorMultiply(tmp_vec1, tmp_vec3);
        for (FESystemUInt i=0; i<n_modes; i++)
        {
            tmp_vec1.zero();
            eig_vec.getColumnVals(sorted_ids[i], 0, s.first-1, tmp_vec1);
            generalized_stiffness_mat.setVal(i, j, tmp_vec1.dotProduct(tmp_vec2));
            generalized_damp_mat.setVal(i, j, tmp_vec1.dotProduct(tmp_vec3));
        }
    }
    
}

