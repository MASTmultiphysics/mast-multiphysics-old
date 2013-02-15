//
//  PlateTest.cpp
//  FESystem
//
//  Created by Manav Bhatia on 9/26/12.
//
//

// FESystem include
#include "TestingIncludes.h"


void calculatePlateStructuralMatrices(FESystemBoolean if_nonlinear, FESystem::Mesh::ElementType elem_type, FESystemUInt n_elem_nodes, const FESystem::DegreeOfFreedom::DegreeOfFreedomMap& dof_map,
                                      const FESystem::Mesh::MeshBase& mesh, FESystem::Numerics::VectorBase<FESystemDouble>& global_sol,
                                      FESystem::Numerics::VectorBase<FESystemDouble>& internal_force,
                                      FESystem::Numerics::VectorBase<FESystemDouble>& external_force,
                                      FESystem::Numerics::MatrixBase<FESystemDouble>& global_stiffness_mat,
                                      FESystem::Numerics::VectorBase<FESystemDouble>& global_mass_vec)
{
    //FESystemDouble E=30.0e6, nu=0.3, rho=2700.0, p_val = 67.5, thick = 0.1; // (Reddy's parameters)
    FESystemDouble E=72.0e9, nu=0.33, rho=2700.0, p_val = 1.0e2, thick = 0.002;
    FESystemBoolean if_mindlin = true;
    FESystemUInt n_plate_dofs, n_elem_dofs;
    
    n_plate_dofs = 3*n_elem_nodes;
    n_elem_dofs = 6*n_elem_nodes;
    
    
    FESystem::Numerics::DenseMatrix<FESystemDouble> elem_mat, plate_elem_mat, vk_elem_mat, pre_stress;
    FESystem::Numerics::LocalVector<FESystemDouble> elem_vec, plate_elem_vec, elem_sol, elem_local_sol, vk_elem_vec;
    std::vector<FESystemUInt> elem_dof_indices;
    elem_mat.resize(n_elem_dofs, n_elem_dofs); vk_elem_mat.resize(5*n_elem_nodes, 5*n_elem_nodes);
    elem_vec.resize(n_elem_dofs); elem_sol.resize(n_elem_dofs); elem_local_sol.resize(5*n_elem_nodes);
    pre_stress.resize(2, 2);
    
    plate_elem_mat.resize(n_plate_dofs, n_plate_dofs); plate_elem_vec.resize(n_plate_dofs); vk_elem_vec.resize(5*n_elem_nodes);
    
    // prepare the quadrature rule and FE for the
    FESystem::Quadrature::TrapezoidQuadrature q_rule_shear, q_rule_bending;
    FESystem::FiniteElement::FELagrange fe, fe_tri6;
    FESystem::Structures::ReissnerMindlinPlate mindlin_plate;
    FESystem::Structures::DKTPlate dkt_plate;
    FESystem::Structures::Membrane membrane;
    FESystem::Structures::VonKarmanStrain2D vk_plate;
    
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
            
            membrane.clear();
            membrane.initialize(*(elems[i]), fe, q_rule_bending, E, nu, rho, thick);
            
            if (if_mindlin)
            {
                mindlin_plate.clear();
                mindlin_plate.initialize(*(elems[i]), fe, q_rule_bending, q_rule_shear, E, nu, rho, thick);
                
                // force
                mindlin_plate.calculateDistributedLoad(p_val, plate_elem_vec);
                mindlin_plate.transformVectorToGlobalSystem(plate_elem_vec, elem_vec);
                dof_map.addToGlobalVector(*(elems[i]), elem_vec, external_force);
                
                vk_plate.clear();
                vk_plate.initialize(*(elems[i]), fe, q_rule_bending, pre_stress, membrane, mindlin_plate);
            }
            else
            {
                dkt_plate.clear();
                dkt_plate.initialize(*(elems[i]), fe, fe_tri6, q_rule_bending, q_rule_bending, E, nu, rho, thick);
                
                // force
                dkt_plate.calculateDistributedLoad(p_val, plate_elem_vec);
                dkt_plate.transformVectorToGlobalSystem(plate_elem_vec, elem_vec);
                dof_map.addToGlobalVector(*(elems[i]), elem_vec, external_force);
                
                vk_plate.clear();
                vk_plate.initialize(*(elems[i]), fe, q_rule_bending, pre_stress, membrane, dkt_plate);
            }
            vk_plate.getActiveElementMatrixIndices(elem_dof_indices);
            elem_sol.getSubVectorValsFromIndices(elem_dof_indices, elem_local_sol);
            
            vk_plate.calculateInternalForceVector(elem_local_sol, vk_elem_vec);
            vk_plate.transformVectorToGlobalSystem(vk_elem_vec, elem_vec);
            dof_map.addToGlobalVector(*(elems[i]), elem_vec, internal_force);
            
            vk_plate.calculateTangentStiffnessMatrix(elem_local_sol, vk_elem_mat);
            vk_plate.transformMatrixToGlobalSystem(vk_elem_mat, elem_mat);
            
        }
        else
        {
            if (if_mindlin)
            {
                mindlin_plate.clear();
                mindlin_plate.initialize(*(elems[i]), fe, q_rule_bending, q_rule_shear, E, nu, rho, thick);
                
                // stiffness matrix
                mindlin_plate.calculateStiffnessMatrix(plate_elem_mat);
                mindlin_plate.transformMatrixToGlobalSystem(plate_elem_mat, elem_mat);
                
                // force
                mindlin_plate.calculateDistributedLoad(p_val, plate_elem_vec);
                mindlin_plate.transformVectorToGlobalSystem(plate_elem_vec, elem_vec);
                dof_map.addToGlobalVector(*(elems[i]), elem_vec, external_force);
                
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
                
                // force
                dkt_plate.calculateDistributedLoad(p_val, plate_elem_vec);
                dkt_plate.transformVectorToGlobalSystem(plate_elem_vec, elem_vec);
                dof_map.addToGlobalVector(*(elems[i]), elem_vec, external_force);
                
                // mass
                dkt_plate.calculateDiagonalMassMatrix(plate_elem_vec);
                dkt_plate.getActiveElementMatrixIndices(elem_dof_indices);
            }
            
            elem_vec.zero();
            elem_vec.addVal(elem_dof_indices, plate_elem_vec);
            dof_map.addToGlobalVector(*(elems[i]), elem_vec, global_mass_vec);
        }
        
        dof_map.addToGlobalMatrix(*(elems[i]), elem_mat, global_stiffness_mat);
    }
}




int plate_analysis_driver(int argc, char * const argv[])
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

    FESystemUInt x_n_divs, y_n_divs;
    std::vector<FESystemDouble> x_div_locations, x_relative_mesh_size_in_div, x_points, y_div_locations, y_relative_mesh_size_in_div, y_points;
    std::vector<FESystemUInt> x_n_subdivs_in_div, y_n_subdivs_in_div;

    x_n_divs = 1;
    x_div_locations.resize(x_n_divs+1);
    x_relative_mesh_size_in_div.resize(x_n_divs+1);
    x_n_subdivs_in_div.resize(x_n_divs);
    
    x_div_locations[0] = 0.0;
    x_div_locations[1] = x_length;
    
    x_relative_mesh_size_in_div[0] = 20.0;
    x_relative_mesh_size_in_div[1] = 1.0;
    
    x_n_subdivs_in_div[0] = 60;
    
    y_n_divs = 1;
    y_div_locations.resize(y_n_divs+1);
    y_relative_mesh_size_in_div.resize(y_n_divs+1);
    y_n_subdivs_in_div.resize(y_n_divs);
    
    y_div_locations[0] = 0.0;
    y_div_locations[1] = y_length;
    
    y_relative_mesh_size_in_div[0] = 1.0;
    y_relative_mesh_size_in_div[1] = 1.0;
    
    y_n_subdivs_in_div[0] = 80;

    distributePoints(x_n_divs, x_div_locations, x_n_subdivs_in_div, x_relative_mesh_size_in_div, x_points);
    distributePoints(y_n_divs, y_div_locations, y_n_subdivs_in_div, y_relative_mesh_size_in_div, y_points);
    
    createPlaneMesh(elem_type, mesh, origin, x_points, y_points, n_elem_nodes, CROSS, true);

    n_elem_dofs = 6*n_elem_nodes;
    
    // set the location of individual nodes
    const std::vector<FESystem::Mesh::Node*>& nodes = mesh.getNodes();
    
    // now add the degrees of freedom
    FESystem::DegreeOfFreedom::DegreeOfFreedomMap dof_map(mesh);
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
        bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(0).global_dof_id[0]); // u-displacement
        bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(1).global_dof_id[0]); // v-displacement
        bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(5).global_dof_id[0]); // theta-z
        if ((nodes[i]->getVal(0) == 0.0) || (nodes[i]->getVal(1) == 0.0) || (nodes[i]->getVal(0) == x_length) || (nodes[i]->getVal(1) == y_length)) // boundary nodes
        {
            //bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(0).global_dof_id[0]); // u-displacement on edges
            //bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(1).global_dof_id[0]); // v-displacement on edges
            bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(2).global_dof_id[0]); // w-displacement on edges
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
    
    //nonlinearSolution(dim, elem_type, n_elem_nodes, mesh, dof_map, nonbc_sparsity_pattern, old_to_new_id_map, nonbc_dofs, global_stiffness_mat, rhs, sol, calculatePlateStructuralMatrices);
    //exit(1);
    
    // calculate the stiffness quantities for the linearized analyses
    calculatePlateStructuralMatrices(false, elem_type, n_elem_nodes, dof_map, mesh, sol, rhs, rhs, global_stiffness_mat, global_mass_vec);
    
    std::vector<FESystemUInt> sorted_ids;
    FESystem::Numerics::LocalVector<FESystemDouble> eig_vals;
    FESystem::Numerics::DenseMatrix<FESystemDouble> eig_vecs;
    
    
    //staticAnalysis(dim, mesh, dof_map, nonbc_sparsity_pattern, old_to_new_id_map, nonbc_dofs, global_stiffness_mat, rhs, sol);
    
    modalAnalysis(dim, mesh, dof_map, nonbc_sparsity_pattern, old_to_new_id_map, nonbc_dofs, global_stiffness_mat, global_mass_vec, n_modes, sorted_ids, eig_vals, eig_vecs);
    
    //transientAnalysis(dim, mesh, dof_map, nonbc_sparsity_pattern, old_to_new_id_map, nonbc_dofs, global_stiffness_mat, global_mass_vec, sorted_ids, eig_vals, eig_vecs);
    
    //flutterSolution(dim, elem_type, n_elem_nodes, mesh, dof_map, nonbc_sparsity_pattern, old_to_new_id_map, nonbc_dofs, global_stiffness_mat, global_mass_vec, n_modes, sorted_ids, eig_vals, eig_vecs, calculatePlatePistonTheoryMatrices);
    
    transientFlutterSolution(dim, elem_type, n_elem_nodes, mesh, dof_map, nonbc_sparsity_pattern, old_to_new_id_map, nonbc_dofs, global_stiffness_mat, global_mass_vec, n_modes, sorted_ids, eig_vals, eig_vecs, calculatePlatePistonTheoryMatrices);
    
    return 0;
}

