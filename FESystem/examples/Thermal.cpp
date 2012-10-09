//
//  ThermalTest.cpp
//  FESystem
//
//  Created by Manav Bhatia on 9/26/12.
//
//


// FESystem include
#include "TestingIncludes.h"






int thermal_1D_analysis_driver(int argc, char * const argv[])
{
    FESystem::Mesh::ElementType elem_type;
    
    // create the mesh object
    FESystem::Mesh::MeshBase mesh;
    
    // create a nx x ny grid of nodes
    FESystemUInt nx, dim, n_elem_nodes, n_modes = 10;
    FESystemDouble x_length;
    FESystem::Geometry::Point origin(3);
    
    nx=15; x_length = 2; dim = 1;
    elem_type = FESystem::Mesh::EDGE2;
    createLineMesh(elem_type, mesh, origin, nx, x_length, n_elem_nodes, INVALID_MESH, true);
    
    // set the location of individual nodes
    const std::vector<FESystem::Mesh::Node*>& nodes = mesh.getNodes();
    
    // now add the degrees of freedom
    FESystem::Base::DegreeOfFreedomMap dof_map(mesh);
    std::string name;
    name = "T"; dof_map.addVariable(name, 0);
    dof_map.reinit(); // distribute the dofs
    
    // create the finite element and initialize the shape functions
    FESystem::Numerics::SparseMatrix<FESystemDouble> global_stiffness_mat, constrained_dof_stiffness_matrix;
    FESystem::Numerics::LocalVector<FESystemDouble> rhs, sol, constrained_sol, constrained_force, global_mass_vec;
    std::vector<FESystemUInt> elem_dof_indices;
    
    global_mass_vec.resize(dof_map.getNDofs());
    rhs.resize(dof_map.getNDofs()); sol.resize(dof_map.getNDofs());
    global_stiffness_mat.resize(dof_map.getSparsityPattern());
    
    // zero solution vector and set the boundary conditions
    sol.zero();
    
    // apply boundary condition and place a load on the last dof
    std::set<FESystemUInt> bc_dofs;
    for (FESystemUInt i=0; i<nodes.size(); i++)
    {
        if (nodes[i]->getVal(0) == 0.0) // boundary nodes
        {
            bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(0).global_dof_id[0]); // T on left edge
            sol.setVal(nodes[i]->getDegreeOfFreedomUnit(0).global_dof_id[0], -10.0);
        }
        if (nodes[i]->getVal(0) == x_length) // boundary nodes
        {
            bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(0).global_dof_id[0]); // T on right edge
            sol.setVal(nodes[i]->getDegreeOfFreedomUnit(0).global_dof_id[0], 10.0);
        }
    }
    
    // now create the vector of ids that do not have bcs
    std::vector<FESystemUInt> nonbc_dofs, bc_dofs_vec;
    for (FESystemUInt i=0; i<dof_map.getNDofs(); i++)
        if (!bc_dofs.count(i))
            nonbc_dofs.push_back(i);
        else
            bc_dofs_vec.push_back(i);
    
    // prepare a map of the old to new ID
    std::vector<FESystemUInt>::const_iterator dof_it=nonbc_dofs.begin(), dof_end=nonbc_dofs.end();
    std::map<FESystemUInt, FESystemUInt> old_to_new_id_map;
    FESystemUInt n=0;
    for ( ; dof_it!=dof_end; dof_it++)
        old_to_new_id_map.insert(std::map<FESystemUInt, FESystemUInt>::value_type(*dof_it, n++));
    
    FESystem::Numerics::SparsityPattern nonbc_sparsity_pattern;
    dof_map.getSparsityPattern().initSparsityPatternForNonConstrainedDOFs(nonbc_dofs, old_to_new_id_map, nonbc_sparsity_pattern);
    
    // calculate the stiffness quantities for the linearized analyses
    calculateBeamStructuralMatrices(false, elem_type, n_elem_nodes, dof_map, mesh, sol, rhs, rhs, global_stiffness_mat, global_mass_vec);
    
    // create the force vector from dof matrices
    constrained_dof_stiffness_matrix.resize(nonbc_dofs.size(), bc_dofs.size());
    constrained_sol.resize(bc_dofs.size()); constrained_force.resize(nonbc_dofs.size());
    sol.getSubVectorValsFromIndices(bc_dofs_vec, constrained_force);
    global_stiffness_mat.getSubMatrixValsFromRowAndColumnIndices(nonbc_dofs, bc_dofs_vec, old_to_new_id_map, constrained_dof_stiffness_matrix);
    constrained_dof_stiffness_matrix.rightVectorMultiply(constrained_sol, constrained_force);
    rhs.add(-1.0, constrained_force);
    
    
    std::vector<FESystemUInt> sorted_ids;
    FESystem::Numerics::LocalVector<FESystemDouble> eig_vals;
    FESystem::Numerics::DenseMatrix<FESystemDouble> eig_vecs;
    
    
    staticAnalysis(dim, mesh, dof_map, nonbc_sparsity_pattern, old_to_new_id_map, nonbc_dofs, global_stiffness_mat, rhs, sol);
    
    modalAnalysis(dim, mesh, dof_map, nonbc_sparsity_pattern, old_to_new_id_map, nonbc_dofs, global_stiffness_mat, global_mass_vec, n_modes, sorted_ids, eig_vals, eig_vecs);
    
    transientAnalysis(dim, mesh, dof_map, nonbc_sparsity_pattern, old_to_new_id_map, nonbc_dofs, global_stiffness_mat, global_mass_vec, sorted_ids, eig_vals, eig_vecs);
    
    return 0;
}

