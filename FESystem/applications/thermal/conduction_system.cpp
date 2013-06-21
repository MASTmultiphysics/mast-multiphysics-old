//
//  assembleEuler.cpp
//  FESystem
//
//  Created by Manav Bhatia on 2/21/13.
//
//


// FESystem includes
#include "thermal/conduction_system.h"
#include "numerics/fem_operator_matrix.h"

// Basic include files
#include "libmesh/equation_systems.h"
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/parameters.h"
#include "libmesh/quadrature.h"
#include "libmesh/dof_map.h"
#include "libmesh/const_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/mesh_function.h"
#include "libmesh/fem_function_base.h"



// Bring in everything from the libMesh namespace
using namespace libMesh;


Real temperature_solution_value(const Point& p,
                                const Parameters& parameters,
                                const std::string& sys_name,
                                const std::string& var_name)
{
    libmesh_assert_equal_to (sys_name, "ConductionSystem");
    
    // since we are initializing the solution, variable values at all points is same
    
    if (var_name == "T")
        return parameters.get<Real> ("Tinit");
    else
        libmesh_assert(false);
}



void init_temperature_variables(EquationSystems& es,
                                const std::string& system_name)
{
    // It is a good idea to make sure we are initializing
    // the proper system.
    libmesh_assert_equal_to (system_name, "ConductionSystem");
    
    // Get a reference to the Convection-Diffusion system object.
    ConductionSystem & system = es.get_system<ConductionSystem>("ConductionSystem");
    
    // Project initial conditions at time 0
    libmesh_assert_equal_to(system.time, 0.);
    
    system.project_solution(temperature_solution_value, NULL, es.parameters);
}


void ConductionSystem::init_data ()
{
    this->use_fixed_solution = true;
    
    // Check the input file for Reynolds number, application type,
    // approximation type
    
    GetPot infile("conduction.in");

    // set parameter values
    Parameters& params = this->get_equation_systems().parameters;
    
    unsigned int o = infile("fe_order", 1);
    std::string fe_family = infile("fe_family", std::string("LAGRANGE"));
    FEFamily fefamily = Utility::string_to_enum<FEFamily>(fe_family);
    
    T_var  = this->add_variable ( "T", static_cast<Order>(o), fefamily);
    this->time_evolving(T_var);
    params.set<Real> ("Tinit") = infile("Tinit", 0.);

    // Useful debugging options
    // Set verify_analytic_jacobians to 1e-6 to use
    this->verify_analytic_jacobians = infile("verify_analytic_jacobians", 0.);
    this->print_jacobians = infile("print_jacobians", false);
    this->print_element_jacobians = infile("print_element_jacobians", false);
    
    // initialize the Dirichlet boundary conditions
    std::multimap<unsigned, ConductionBoundaryConditionType>::const_iterator
    bc_it = _boundary_condition.begin(), bc_end = _boundary_condition.end();
    
    std::set<boundary_id_type> temperature_boundary;
    
    for ( ; bc_it!=bc_end; bc_it++)
    {
        switch (bc_it->second)
        {
            case TEMPERATURE: // specified temperature
                temperature_boundary.insert(bc_it->first);
                break;
                
            default:
                break;
        }
    }
    
    if (temperature_boundary.size() > 0)
    {
        // Dirichlet Boundary condition
        ConstFunction<Real> constant_function(0.);
        this->get_dof_map().add_dirichlet_boundary(DirichletBoundary(temperature_boundary,
                                                                     std::vector<unsigned int>(1,T_var),
                                                                     &constant_function));
    }

    // Do the parent's initialization after variables and boundary constraints are defined
    FEMSystem::init_data();
}



void ConductionSystem::init_context(DiffContext &context)
{
    FEMContext &c = libmesh_cast_ref<FEMContext&>(context);
    
    FEBase* elem_fe;
    
    c.get_element_fe(T_var, elem_fe);
    elem_fe->get_JxW();
    elem_fe->get_phi();
    elem_fe->get_dphi();
    elem_fe->get_xyz();
    
    
    c.get_side_fe( T_var, elem_fe);
    elem_fe->get_JxW();
    elem_fe->get_phi();
    elem_fe->get_xyz();
}



bool ConductionSystem::element_time_derivative (bool request_jacobian,
                                           DiffContext &context)
{
    FEMContext &c = libmesh_cast_ref<FEMContext&>(context);
    
    FEBase* elem_fe;

    c.get_element_fe(T_var, elem_fe);
    
    // Element Jacobian * quadrature weights for interior integration
    const std::vector<Real> &JxW = elem_fe->get_JxW();
        
    // The number of local degrees of freedom in each variable
    const unsigned int n_dofs = c.get_dof_indices( T_var ).size();

    // The subvectors and submatrices we need to fill:
    DenseMatrix<Number>& Kmat = c.get_elem_jacobian();
    DenseVector<Number>& Fvec = c.get_elem_residual();
    
    const unsigned int n_qpoints = (c.get_element_qrule())->n_points(),
    dim = this->get_mesh().mesh_dimension();

    std::vector<FEMOperatorMatrix> dB_mat(dim);
    
    DenseMatrix<Real> tmp_mat_n1n2, tmp_mat_n2n2;
    DenseVector<Real> tmp_vec_n1, tmp_vec_n2;

    tmp_mat_n1n2.resize(1, n_dofs); tmp_mat_n2n2.resize(n_dofs, n_dofs);
    tmp_vec_n1.resize(1); tmp_vec_n2.resize(n_dofs);
    

    const std::vector<std::vector<RealVectorValue> >& dphi = elem_fe->get_dphi();
    
    for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
        // initialize the FEM operator matrices
        for ( unsigned int i_dim=0; i_dim<dim; i_dim++ )
        {
            for ( unsigned int i_nd=0; i_nd<n_dofs; i_nd++ )
                tmp_vec_n2(i_nd) = dphi[i_nd][qp](i_dim);
            dB_mat[i_dim].reinit(1, tmp_vec_n2);
        }

        // now calculate the flux value
        for (unsigned int i_dim=0; i_dim<dim; i_dim++)
        {
            dB_mat[i_dim].vector_mult(tmp_vec_n1, c.elem_solution);      // dB/dx_i T
            dB_mat[i_dim].vector_mult_transpose(tmp_vec_n2, tmp_vec_n1); // dB/dx_i^T dB/dx_i T
            Fvec.add(-JxW[qp]*_k_thermal, tmp_vec_n2);
        }
        
        
        if (request_jacobian && c.elem_solution_derivative)
        {
            for (unsigned int i_dim=0; i_dim<dim; i_dim++)
            {
                dB_mat[i_dim].right_multiply_transpose(tmp_mat_n2n2, dB_mat[i_dim]); // dB/dx_i^T dB/dx_i
                Kmat.add(-JxW[qp]*_k_thermal, tmp_mat_n2n2);
            }
        }
    } // end of the quadrature point qp-loop

//    std::cout << "inside element time derivative " << std::endl;
//    c.elem->print_info();
//    std::cout << "sol: " << std::endl; c.elem_solution.print(std::cout);
//    std::cout << "res: " << std::endl; Fvec.print(std::cout);
//    if (request_jacobian && c.elem_solution_derivative)
//        Kmat.print(std::cout);
    
    return request_jacobian;
}



bool ConductionSystem::side_time_derivative (bool request_jacobian,
                                        DiffContext &context)
{
    FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

    // check for the boundary tags to check if the boundary condition needs to be applied to this element
    
    ConductionBoundaryConditionType bc_type;
    
    std::multimap<unsigned int, ConductionBoundaryConditionType>::const_iterator
    bc_it     = this->_boundary_condition.begin(),
    bc_end    = this->_boundary_condition.end();
    
    

    
//    std::cout << "inside side constraint " << std::endl;
//    std::cout << "elem solution" << std::endl; c.elem_solution.print(std::cout);
//    std::cout << if_inf_bc << "  " << if_wall_bc << std::endl;
//    std::cout << "bc vec: " << std::endl; Fvec.print(std::cout);
//    if (request_jacobian && c.elem_solution_derivative)
//        Kmat.print(std::cout);

    return request_jacobian;
}





bool ConductionSystem::mass_residual (bool request_jacobian,
                                 DiffContext &context)
{
    FEMContext &c = libmesh_cast_ref<FEMContext&>(context);
    
    FEBase* elem_fe;
    
    c.get_element_fe(T_var, elem_fe);
    
    // Element Jacobian * quadrature weights for interior integration
    const std::vector<Real> &JxW = elem_fe->get_JxW();
    
    // The number of local degrees of freedom in each variable
    const unsigned int n_dofs = c.get_dof_indices( T_var ).size();
    
    // The subvectors and submatrices we need to fill:
    DenseMatrix<Number>& Kmat = c.get_elem_jacobian();
    DenseVector<Number>& Fvec = c.get_elem_residual();
    
    const unsigned int n_qpoints = (c.get_element_qrule())->n_points(),
    dim = this->get_mesh().mesh_dimension();
    
    FEMOperatorMatrix B_mat;
    
    DenseMatrix<Real> tmp_mat_n1n2, tmp_mat_n2n2;
    DenseVector<Real> tmp_vec_n1, tmp_vec_n2;
    
    tmp_mat_n1n2.resize(1, n_dofs); tmp_mat_n2n2.resize(n_dofs, n_dofs);
    tmp_vec_n1.resize(1); tmp_vec_n2.resize(n_dofs);
    
    
    const std::vector<std::vector<Real> >& phi = elem_fe->get_phi();
    
    for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
        // initialize the FEM operator matrices
        for ( unsigned int i_dim=0; i_dim<dim; i_dim++ )
        {
            for ( unsigned int i_nd=0; i_nd<n_dofs; i_nd++ )
                tmp_vec_n2(i_nd) = phi[i_nd][qp];
            B_mat.reinit(1, tmp_vec_n2);
        }
        
        // now calculate the residual contribution
        B_mat.vector_mult(tmp_vec_n1, c.elem_solution);      // B T
        B_mat.vector_mult_transpose(tmp_vec_n2, tmp_vec_n1); // B^T B T
        Fvec.add(JxW[qp]*_cp_thermal, tmp_vec_n2);
        
        
        if (request_jacobian && c.elem_solution_derivative)
        {
            B_mat.right_multiply_transpose(tmp_mat_n2n2, B_mat); // B^T B_i
            Kmat.add(JxW[qp]*_cp_thermal, tmp_mat_n2n2);
        }
    } // end of the quadrature point qp-loop
    
    //    std::cout << "inside mass_residual " << std::endl;
    //    c.elem->print_info();
    //    std::cout << "sol: " << std::endl; c.elem_solution.print(std::cout);
    //    std::cout << "res: " << std::endl; Fvec.print(std::cout);
    //    if (request_jacobian && c.elem_solution_derivative)
    //        Kmat.print(std::cout);
    
    return request_jacobian;
}


