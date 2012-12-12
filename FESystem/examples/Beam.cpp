//
//  BeamTest.cpp
//  FESystem
//
//  Created by Manav Bhatia on 9/26/12.
//
//

// FESystem include
#include "TestingIncludes.h"



void calculateBeamStructuralMatrices(FESystemBoolean if_nonlinear, FESystem::Mesh::ElementType elem_type, FESystemUInt n_elem_nodes, const FESystem::Base::DegreeOfFreedomMap& dof_map,
                                     const FESystem::Mesh::MeshBase& mesh, FESystem::Numerics::VectorBase<FESystemDouble>& global_sol,
                                     FESystem::Numerics::VectorBase<FESystemDouble>& internal_force,
                                     FESystem::Numerics::VectorBase<FESystemDouble>& external_force,
                                     FESystem::Numerics::MatrixBase<FESystemDouble>& global_stiffness_mat,
                                     FESystem::Numerics::VectorBase<FESystemDouble>& global_mass_vec)
{
    
    //FESystemDouble E=30.0e6, nu=0.33, rho=2700.0, v_p_val = 0.0, w_p_val = 9.0e0, I_tr = 1.0/12.0, I_ch = 1.0/12.0, area = 1.0; // (Reddy's parameters)
    FESystemDouble E=72.0e9, nu=0.33, rho=2700.0, v_p_val = 0.0, w_p_val = 1.0e0, I_tr = 6.667e-9, I_ch = 1.6667e-9, area = 2.0e-4;
    FESystemBoolean if_timoshenko_beam = true;
    FESystemUInt n_beam_dofs, n_elem_dofs;
    
    n_beam_dofs = 4*n_elem_nodes;
    n_elem_dofs = 6*n_elem_nodes;
    
    
    FESystem::Numerics::DenseMatrix<FESystemDouble> elem_mat, beam_elem_mat, vk_elem_mat;
    FESystem::Numerics::LocalVector<FESystemDouble> elem_vec, beam_elem_vec, elem_sol, elem_local_sol, vk_elem_vec;
    std::vector<FESystemUInt> elem_dof_indices;
    elem_mat.resize(n_elem_dofs, n_elem_dofs); vk_elem_mat.resize(5*n_elem_nodes, 5*n_elem_nodes);
    elem_vec.resize(n_elem_dofs); elem_sol.resize(n_elem_dofs); elem_local_sol.resize(5*n_elem_nodes);
    
    beam_elem_mat.resize(n_beam_dofs, n_beam_dofs); beam_elem_vec.resize(n_beam_dofs); vk_elem_vec.resize(5*n_elem_nodes);
    
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
    
    internal_force.zero(); external_force.zero(); global_stiffness_mat.zero(); global_mass_vec.zero();
    
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
            
            bar.clear();
            bar.initialize(*(elems[i]), fe, q_rule_bending, E, nu, rho, area);
            
            if (if_timoshenko_beam)
            {
                timoshenko_beam.clear();
                timoshenko_beam.initialize(*(elems[i]), fe, q_rule_bending, q_rule_shear, E, nu, rho, I_tr, I_ch, area);
                
                // force
                timoshenko_beam.calculateDistributedLoad(v_p_val, w_p_val, beam_elem_vec);
                timoshenko_beam.transformVectorToGlobalSystem(beam_elem_vec, elem_vec);
                dof_map.addToGlobalVector(*(elems[i]), elem_vec, external_force);
                
                vk_beam.clear();
                vk_beam.initialize(*(elems[i]), fe, q_rule_bending, 0.0, bar, timoshenko_beam);
            }
            else
            {
                euler_beam.clear();
                euler_beam.initialize(*(elems[i]), fe, q_rule_bending, E, nu, rho, I_tr, I_ch, area);
                
                // force
                euler_beam.calculateDistributedLoad(v_p_val, w_p_val, beam_elem_vec);
                euler_beam.transformVectorToGlobalSystem(beam_elem_vec, elem_vec);
                dof_map.addToGlobalVector(*(elems[i]), elem_vec, external_force);
                
                vk_beam.clear();
                vk_beam.initialize(*(elems[i]), fe, q_rule_bending, 0.0, bar, euler_beam);
            }
            vk_beam.getActiveElementMatrixIndices(elem_dof_indices);
            elem_sol.getSubVectorValsFromIndices(elem_dof_indices, elem_local_sol);
            
            vk_beam.calculateInternalForceVector(elem_local_sol, vk_elem_vec);
            vk_beam.transformVectorToGlobalSystem(vk_elem_vec, elem_vec);
            dof_map.addToGlobalVector(*(elems[i]), elem_vec, internal_force);
            
            vk_beam.calculateTangentStiffnessMatrix(elem_local_sol, vk_elem_mat);
            vk_beam.transformMatrixToGlobalSystem(vk_elem_mat, elem_mat);
            
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
                
                // force
                timoshenko_beam.calculateDistributedLoad(v_p_val, w_p_val, beam_elem_vec);
                timoshenko_beam.transformVectorToGlobalSystem(beam_elem_vec, elem_vec);
                dof_map.addToGlobalVector(*(elems[i]), elem_vec, external_force);
                
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
                
                // force
                euler_beam.calculateDistributedLoad(v_p_val, w_p_val, beam_elem_vec);
                euler_beam.transformVectorToGlobalSystem(beam_elem_vec, elem_vec);
                dof_map.addToGlobalVector(*(elems[i]), elem_vec, external_force);
                
                // mass
                euler_beam.calculateDiagonalMassMatrix(beam_elem_vec);
                euler_beam.getActiveElementMatrixIndices(elem_dof_indices);
            }
            
            elem_vec.zero();
            elem_vec.addVal(elem_dof_indices, beam_elem_vec);
            dof_map.addToGlobalVector(*(elems[i]), elem_vec, global_mass_vec);
        }
        
        dof_map.addToGlobalMatrix(*(elems[i]), elem_mat, global_stiffness_mat);
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
    
    nx=15; x_length = 10.8; dim = 1; n_modes = 5;
    elem_type = FESystem::Mesh::EDGE2;

    FESystemUInt x_n_divs;
    std::vector<FESystemDouble> x_div_locations, x_relative_mesh_size_in_div, x_points;
    std::vector<FESystemUInt> x_n_subdivs_in_div;
    
    x_n_divs = 1;
    x_div_locations.resize(x_n_divs+1);
    x_relative_mesh_size_in_div.resize(x_n_divs+1);
    x_n_subdivs_in_div.resize(x_n_divs);
    
    x_div_locations[0] = 0.0;
    x_div_locations[1] = x_length;
    
    x_relative_mesh_size_in_div[0] = 1.0;
    x_relative_mesh_size_in_div[1] = 1.0;
    
    x_n_subdivs_in_div[0] = 60;
    
    distributePoints(x_n_divs, x_div_locations, x_n_subdivs_in_div, x_relative_mesh_size_in_div, x_points);
    
    createLineMesh(elem_type, mesh, origin, x_points, n_elem_nodes, INVALID_MESH, true);
    
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
        bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(1).global_dof_id[0]); // v-disp
        bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(3).global_dof_id[0]); // theta-x
        bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(5).global_dof_id[0]); // theta-z
        if ((nodes[i]->getVal(0) == 0.0) || (nodes[i]->getVal(0) == x_length)) // boundary nodes
        {
            //bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(0).global_dof_id[0]); // u-displacement on the edges
            //bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(1).global_dof_id[0]); // v-displacement on the edges
            bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(2).global_dof_id[0]); // w-displacement on the edges
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
    
    //nonlinearSolution(dim, elem_type, n_elem_nodes, mesh, dof_map, nonbc_sparsity_pattern, old_to_new_id_map, nonbc_dofs, global_stiffness_mat, rhs, sol, calculateBeamStructuralMatrices);
    //exit(1);
    
    // calculate the stiffness quantities for the linearized analyses
    calculateBeamStructuralMatrices(false, elem_type, n_elem_nodes, dof_map, mesh, sol, rhs, rhs, global_stiffness_mat, global_mass_vec);
    
    std::vector<FESystemUInt> sorted_ids;
    FESystem::Numerics::LocalVector<FESystemDouble> eig_vals;
    FESystem::Numerics::DenseMatrix<FESystemDouble> eig_vecs;
    
    
    //staticAnalysis(dim, mesh, dof_map, nonbc_sparsity_pattern, old_to_new_id_map, nonbc_dofs, global_stiffness_mat, rhs, sol);
    
    modalAnalysis(dim, mesh, dof_map, nonbc_sparsity_pattern, old_to_new_id_map, nonbc_dofs, global_stiffness_mat, global_mass_vec, n_modes, sorted_ids, eig_vals, eig_vecs);
    
    //transientAnalysis(dim, mesh, dof_map, nonbc_sparsity_pattern, old_to_new_id_map, nonbc_dofs, global_stiffness_mat, global_mass_vec, sorted_ids, eig_vals, eig_vecs);
    
    //flutterSolution(dim, elem_type, n_elem_nodes, mesh, dof_map, nonbc_sparsity_pattern, old_to_new_id_map, nonbc_dofs, global_stiffness_mat, global_mass_vec, n_modes, sorted_ids, eig_vals, eig_vecs, calculateBeamPistonTheoryMatrices);
    
    //transientFlutterSolution(dim, elem_type, n_elem_nodes, mesh, dof_map, nonbc_sparsity_pattern, old_to_new_id_map, nonbc_dofs, global_stiffness_mat, global_mass_vec, n_modes, sorted_ids, eig_vals, eig_vecs, calculateBeamPistonTheoryMatrices);
    
    return 0;
}





