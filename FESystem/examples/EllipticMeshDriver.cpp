//
//  EllipticMeshDriver.cpp
//  FESystem
//
//  Created by Manav Bhatia on 12/23/12.
//
//

// C++ includes
#include <map>
#include <algorithm>

// FESystem include
#include "TestingIncludes.h"
#include "Disciplines/Meshing/EllipticMeshGenerator.h"

// TBB includes
#include "tbb/tbb.h"
#include "tbb/mutex.h"
#include "tbb/spin_mutex.h"


const FESystem::Mesh::ElementType elem_type = FESystem::Mesh::QUAD4;

tbb::mutex assembly_mutex;

class AssembleEllipticMeshingElementMatrices
{
public:
    AssembleEllipticMeshingElementMatrices(const std::vector<FESystem::Mesh::ElemBase*>& e,
                                           const FESystemUInt ne_dofs,
                                           const FESystem::Base::DegreeOfFreedomMap& dof,
                                           const FESystem::Quadrature::QuadratureBase& q,
                                           const FESystem::Quadrature::QuadratureBase& q_b,
                                           const FESystem::Numerics::VectorBase<FESystemDouble>& s,
                                           const FESystemBoolean if_jacobian,
                                           FESystem::Numerics::VectorBase<FESystemDouble>& r,
                                           FESystem::Numerics::MatrixBase<FESystemDouble>& stiff):
    elems(e),
    n_elem_dofs(ne_dofs),
    dof_map(dof),
    q_rule(q),
    q_boundary(q_b),
    sol(s),
    if_calculate_jacobian(if_jacobian),
    residual(r),
    global_stiffness_mat(stiff)
    {
    }
    
    
    void operator() (const tbb::blocked_range<FESystemUInt>& r) const
    {
        this->assemble(r.begin(), r.end());
    }
    
    
    void assembleAll() const
    {
        this->assemble(0, elems.size());
    }
    
    
    void assemble(FESystemUInt e1, FESystemUInt e2) const
    {
        FESystem::Numerics::DenseMatrix<FESystemDouble> elem_mat;
        FESystem::Numerics::LocalVector<FESystemDouble> elem_vec, elem_sol;
        elem_mat.resize(n_elem_dofs, n_elem_dofs); elem_vec.resize(n_elem_dofs); elem_sol.resize(n_elem_dofs);
        
        FESystem::FiniteElement::FELagrange fe;
        FESystem::Meshing::EllipticMeshGenerator meshing_elem;
        
        for (FESystemUInt i=e1; i!=e2; i++)
        {
            fe.clear();
            fe.reinit(*(elems[i]));
            elem_mat.zero();
            elem_vec.zero();
            
            dof_map.getFromGlobalVector(*(elems[i]), sol, elem_sol);
            
            meshing_elem.clear();
            meshing_elem.initialize(*(elems[i]), fe, q_rule, elem_sol, 0.0, 0.0);
            
            meshing_elem.calculateResidual(elem_vec);
            
            if (if_calculate_jacobian)
                meshing_elem.calculateTangentMatrix(elem_mat);
            
            {
                tbb::mutex::scoped_lock my_lock(assembly_mutex);
                dof_map.addToGlobalVector(*(elems[i]), elem_vec, residual);
                if (if_calculate_jacobian)
                    dof_map.addToGlobalMatrix(*(elems[i]), elem_mat, global_stiffness_mat);
            }
        }
    }
    

protected:
    
    const std::vector<FESystem::Mesh::ElemBase*>& elems;
    const FESystemUInt n_elem_dofs;
    const FESystem::Base::DegreeOfFreedomMap& dof_map;
    const FESystem::Quadrature::QuadratureBase& q_rule;
    const FESystem::Quadrature::QuadratureBase& q_boundary;
    const FESystem::Numerics::VectorBase<FESystemDouble>& sol;
    const FESystemBoolean if_calculate_jacobian;
    FESystem::Numerics::VectorBase<FESystemDouble>& residual;
    FESystem::Numerics::MatrixBase<FESystemDouble>& global_stiffness_mat;
};





void calculateEllipticMeshingMatrices(const FESystem::Mesh::MeshBase& mesh, FESystemUInt n_elem_dofs, const FESystem::Base::DegreeOfFreedomMap& dof_map, const FESystem::Numerics::VectorBase<FESystemDouble>& sol, FESystem::Numerics::VectorBase<FESystemDouble>& res, FESystem::Numerics::MatrixBase<FESystemDouble>& mat, FESystemBoolean if_evaluate_jac)
{
    // prepare the quadrature rule and FE for the
    FESystem::Quadrature::TrapezoidQuadrature q_rule, q_boundary_rule;
    
    switch (elem_type)
    {
        case FESystem::Mesh::QUAD4:
        case FESystem::Mesh::TRI3:
            q_rule.init(2, 1);
            q_boundary_rule.init(1, 1);

            break;
            
        case FESystem::Mesh::QUAD9:
        case FESystem::Mesh::TRI6:
            q_rule.init(2, 4);
            q_boundary_rule.init(1, 4);
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
    }
    
    res.zero(); mat.zero();
    
    const std::vector<FESystem::Mesh::ElemBase*>& elems = mesh.getElements();
    
    tbb::parallel_for(tbb::blocked_range<FESystemUInt>(0, elems.size()),
                      AssembleEllipticMeshingElementMatrices(elems, n_elem_dofs, dof_map, q_rule, q_boundary_rule, sol, if_evaluate_jac, res, mat));
//    AssembleEllipticMeshingElementMatrices a(elems, n_elem_dofs, dof_map, q_rule, q_boundary_rule, sol, if_evaluate_jac, res, mat);
//    a.assembleAll();
}




void ellipticMeshingNonlinearSolution(FESystemUInt dim, FESystem::Mesh::ElementType elem_type, FESystemUInt n_elem_nodes, const
                                      FESystem::Mesh::MeshBase& mesh, const FESystem::Base::DegreeOfFreedomMap& dof_map,
                                      const FESystem::Numerics::SparsityPattern& nonbc_sparsity_pattern,
                                      const std::map<FESystemUInt, FESystemUInt>& old_to_new_id_map,
                                      const std::vector<FESystemUInt>& bc_dofs,
                                      const std::vector<FESystemUInt>& nonbc_dofs,
                                      FESystem::Numerics::MatrixBase<FESystemDouble>& stiff_mat,
                                      FESystem::Numerics::VectorBase<FESystemDouble>& rhs, FESystem::Numerics::VectorBase<FESystemDouble>& sol)
{
    const std::vector<FESystem::Mesh::Node*>& nodes = mesh.getNodes();
    const FESystemUInt n_elem_dofs = 2*n_elem_nodes;
    
    
    FESystem::Numerics::SparseMatrix<FESystemDouble> reduced_stiff_mat;
    FESystem::Numerics::LocalVector<FESystemDouble> reduced_load_vec, reduced_sol_vec, internal_force, dummy;
    
    internal_force.resize(dof_map.getNDofs());
    reduced_stiff_mat.resize(nonbc_sparsity_pattern);
    reduced_load_vec.resize(nonbc_sparsity_pattern.getNDOFs()); reduced_sol_vec.resize(nonbc_sparsity_pattern.getNDOFs());
    stiff_mat.getSubMatrixValsFromRowAndColumnIndices(nonbc_dofs, nonbc_dofs, old_to_new_id_map, reduced_stiff_mat);
    
    FESystem::LinearSolvers::LUFactorizationLinearSolver<FESystemDouble> linear_solver;
    FESystem::NonlinearSolvers::NewtonIterationNonlinearSolver<FESystemDouble> nonlinear_solver;
    
    nonlinear_solver.initialize(reduced_stiff_mat, linear_solver);
    nonlinear_solver.setConvergenceLimits(50, 1.0e-10);
    
    
    FESystem::NonlinearSolvers::NonlinearSolverCallBack call_back = nonlinear_solver.getCurrentCallBack();
    
    while (call_back != FESystem::NonlinearSolvers::SOLUTION_CONVERGED)
    {
        switch (call_back)
        {
            case FESystem::NonlinearSolvers::WAITING_TO_START:
                // nothing to be done
                break;
                
            case FESystem::NonlinearSolvers::SET_INITIAL_GUESS:
                sol.getSubVectorValsFromIndices(nonbc_dofs, nonlinear_solver.getCurrentSolution());
                break;
                
            case FESystem::NonlinearSolvers::EVALUATE_RESIDUAL:
            {
                // zero the solution vectors
                rhs.zero();
                // get the latest solution vector
                sol.setSubVectorValsFromIndices(nonbc_dofs, nonlinear_solver.getCurrentSolution());
                calculateEllipticMeshingMatrices(mesh, n_elem_dofs, dof_map, sol, rhs, stiff_mat, false);
                rhs.getSubVectorValsFromIndices(nonbc_dofs, nonlinear_solver.getResidualVector());
            }
                break;
                
            case FESystem::NonlinearSolvers::EVALUATE_JACOBIAN:
            case FESystem::NonlinearSolvers::EVALUATE_RESIDUAL_AND_JACOBIAN:
            {
                // zero the solution vectors
                rhs.zero(); 
                // get the latest solution vector
                sol.setSubVectorValsFromIndices(nonbc_dofs, nonlinear_solver.getCurrentSolution());
                calculateEllipticMeshingMatrices(mesh, n_elem_dofs, dof_map, sol, rhs, stiff_mat, true);
                rhs.getSubVectorValsFromIndices(nonbc_dofs, nonlinear_solver.getResidualVector());
                stiff_mat.getSubMatrixValsFromRowAndColumnIndices(nonbc_dofs, nonbc_dofs, old_to_new_id_map, reduced_stiff_mat);
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
    
    // write the solution
    // for output
    FESystem::OutputProcessor::VtkOutputProcessor output;
    std::fstream output_file;
    std::vector<FESystemUInt> vars(2); vars[0]=0; vars[1]=1;; // write all solutions
    
    output_file.open("sol.vtk",std::fstream::out);
    output.writeMesh(output_file, mesh, dof_map);
    output.writeSolution(output_file, "Sol", mesh, dof_map, vars, sol);
    output_file.close();
}




FESystemBoolean compareInt(FESystemUInt i, FESystemUInt j)
{
    return i<j;
}



int ellipticMeshingDriver(int argc, char * const argv[])
{
    
    // create the mesh object
    FESystem::Mesh::MeshBase mesh;
    
    // create a nx x ny grid of nodes
    FESystemUInt dim, n_elem_nodes;
    FESystemDouble x_length, y_length;
    FESystem::Geometry::Point origin(3);
    
    x_length = 1.0; y_length = 1.0; dim = 2;
    
    
    FESystemUInt x_n_divs, y_n_divs;
    std::vector<FESystemDouble> x_div_locations, x_relative_mesh_size_in_div, x_points, y_div_locations, y_relative_mesh_size_in_div, y_points;
    std::vector<FESystemUInt> x_n_subdivs_in_div, y_n_subdivs_in_div;
    
    x_n_divs = 1;
    x_div_locations.resize(x_n_divs+1);
    x_relative_mesh_size_in_div.resize(x_n_divs+1);
    x_n_subdivs_in_div.resize(x_n_divs);
    
    x_div_locations[0] = 0.0;
    x_div_locations[1] = x_length;
    
    x_relative_mesh_size_in_div[0] = 1.0;
    x_relative_mesh_size_in_div[1] = 1.0;
    
    x_n_subdivs_in_div[0] = 30;
    
    y_n_divs = 1;
    y_div_locations.resize(y_n_divs+1);
    y_relative_mesh_size_in_div.resize(y_n_divs+1);
    y_n_subdivs_in_div.resize(y_n_divs);
    
    y_div_locations[0] = 0.0;
    y_div_locations[1] = y_length;
    
    y_relative_mesh_size_in_div[0] = 1.0;
    y_relative_mesh_size_in_div[1] = 1.0;
    
    y_n_subdivs_in_div[0] = 20;
    
    distributePoints(x_n_divs, x_div_locations, x_n_subdivs_in_div, x_relative_mesh_size_in_div, x_points);
    distributePoints(y_n_divs, y_div_locations, y_n_subdivs_in_div, y_relative_mesh_size_in_div, y_points);
    
    createPlaneMesh(elem_type, mesh, origin, x_points, y_points, n_elem_nodes, CROSS, true);
        
    // now add the degrees of freedom
    FESystem::Base::DegreeOfFreedomMap dof_map(mesh);
    std::string name;
    name = "x"; dof_map.addVariable(name, 0);
    name = "y"; dof_map.addVariable(name, 0);
    dof_map.reinit(); // distribute the dofs
    
    // create the finite element and initialize the shape functions
    FESystem::Numerics::SparseMatrix<FESystemDouble> global_stiffness_mat;
    FESystem::Numerics::LocalVector<FESystemDouble> rhs, sol, global_mass_vec;
    std::vector<FESystemUInt> elem_dof_indices;
    
    global_mass_vec.resize(dof_map.getNDofs());
    rhs.resize(dof_map.getNDofs()); sol.resize(dof_map.getNDofs());
    global_stiffness_mat.resize(dof_map.getSparsityPattern());
    
    std::vector<FESystem::Mesh::Node*>& nodes = mesh.getNodes();
    FESystemDouble xval, yval, domain_x = 10.0, domain_y = 5.0, y_height=domain_y*0.2;
    FESystemUInt x_id, y_id;
    std::vector<FESystemUInt> nonbc_dofs, bc_dofs;
    // apply boundary condition and set the boundary dof locations
    for (FESystemUInt i=0; i<nodes.size(); i++)
    {
        xval = nodes[i]->getVal(0);
        yval = nodes[i]->getVal(1);
        x_id = nodes[i]->getDegreeOfFreedomUnit(0).global_dof_id[0];
        y_id = nodes[i]->getDegreeOfFreedomUnit(1).global_dof_id[0];
        
        if (xval == 0.0) // left boundary
        {
            sol.setVal(x_id, 0.0);
            sol.setVal(y_id, yval*domain_y); // initial value, allowed to float based on zero flux
            bc_dofs.push_back(x_id);
            bc_dofs.push_back(y_id);
        }
        else if (xval == x_length) // right boundary
        {
            sol.setVal(x_id, domain_x);
            sol.setVal(y_id, yval*domain_y);
            bc_dofs.push_back(x_id);
            bc_dofs.push_back(y_id);
        }
        else if (yval == 0.0) // bottom boundary
        {
            sol.setVal(x_id, xval*domain_x);
            sol.setVal(y_id, y_height*(.25-pow(xval-0.5,2))/.25);
            bc_dofs.push_back(x_id);
            bc_dofs.push_back(y_id);
        }
        else if (yval == y_length) // top boundary
        {
            sol.setVal(x_id, xval*domain_x);
            sol.setVal(y_id, domain_y);
            bc_dofs.push_back(x_id);
            bc_dofs.push_back(y_id);
        }
        else
        {
            sol.setVal(x_id, xval*domain_x);
            sol.setVal(y_id, yval*domain_y);
            nonbc_dofs.push_back(x_id);
            nonbc_dofs.push_back(y_id);
        }
    }
    
    std::sort(bc_dofs.begin(), bc_dofs.end(), compareInt);
    std::sort(nonbc_dofs.begin(), nonbc_dofs.end(), compareInt);
    
    
    // prepare a map of the old to new ID
    std::vector<FESystemUInt>::const_iterator dof_it=nonbc_dofs.begin(), dof_end=nonbc_dofs.end();
    std::map<FESystemUInt, FESystemUInt> old_to_new_id_map;
    FESystemUInt n=0;
    for ( ; dof_it!=dof_end; dof_it++)
        old_to_new_id_map.insert(std::map<FESystemUInt, FESystemUInt>::value_type(*dof_it, n++));
    
    FESystem::Numerics::SparsityPattern nonbc_sparsity_pattern;
    dof_map.getSparsityPattern().initSparsityPatternForNonConstrainedDOFs(nonbc_dofs, old_to_new_id_map, nonbc_sparsity_pattern);
        
    ellipticMeshingNonlinearSolution(dim, elem_type, n_elem_nodes, mesh, dof_map, nonbc_sparsity_pattern, old_to_new_id_map, bc_dofs, nonbc_dofs, global_stiffness_mat, rhs, sol);
    
    // now set the nodal location
    std::vector<FESystem::Mesh::Node*>::iterator node_it=nodes.begin(), node_end=nodes.end();
    for ( ; node_it != node_end; node_it++)
    {
        (*node_it)->setVal(0, sol.getVal((*node_it)->getDegreeOfFreedomUnit(0).global_dof_id[0]));
        (*node_it)->setVal(1, sol.getVal((*node_it)->getDegreeOfFreedomUnit(1).global_dof_id[0]));
    }
    
    std::vector<FESystem::Mesh::ElemBase*>& elems = mesh.getElements();
    std::vector<FESystem::Mesh::ElemBase*>::iterator elem_it=elems.begin(), elem_end=elems.end();
    for ( ; elem_it != elem_end; elem_it++)
        (*elem_it)->updateAfterMeshDeformation();
    
    // write the solution
    
    std::fstream output_file;
    FESystem::OutputProcessor::VtkOutputProcessor output;
    output_file.open("mesh.vtk",std::fstream::out);
    output.writeMesh(output_file, mesh, dof_map);
    output_file.close();

    return 0;
}


