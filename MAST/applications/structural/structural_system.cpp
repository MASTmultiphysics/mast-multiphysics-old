//
//  Structural_system.cpp
//  FESystem
//
//  Created by Manav Bhatia on 6/24/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//


// FESystem includes
#include "structural/structural_system.h"
#include "numerics/fem_operator_matrix.h"

// Basic include files
#include "libmesh/equation_systems.h"
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/parameters.h"
#include "libmesh/quadrature.h"
#include "libmesh/dof_map.h"
#include "libmesh/zero_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/mesh_function.h"
#include "libmesh/fem_function_base.h"

#ifndef LIBMESH_USE_COMPLEX_NUMBERS


// Bring in everything from the libMesh namespace
using namespace libMesh;


Real structural_solution_value(const Point& p,
                               const Parameters& parameters,
                               const std::string& sys_name,
                               const std::string& var_name)
{
    libmesh_assert_equal_to (sys_name, "StructuralSystem");
    
    // zero initial displacement for the entire structure
    return 0.;
}



void init_structural_variables(EquationSystems& es,
                               const std::string& system_name)
{
    // It is a good idea to make sure we are initializing
    // the proper system.
    libmesh_assert_equal_to (system_name, "StructuralSystem");
    
    // Get a reference to the Convection-Diffusion system object.
    StructuralSystem & system = es.get_system<StructuralSystem>("StructuralSystem");
    
    // Project initial conditions at time 0
    libmesh_assert_equal_to(system.time, 0.);
    
    system.project_solution(structural_solution_value, NULL, es.parameters);
}


void StructuralSystem::init_data ()
{
    this->use_fixed_solution = true;
    
    // Check the input file for Reynolds number, application type,
    // approximation type
    
    GetPot infile("structural.in");
    
    // set parameter values
    Parameters& params = this->get_equation_systems().parameters;

    const unsigned int dim = this->get_mesh().mesh_dimension();
    
    unsigned int o = infile("fe_order", 1);
    std::string fe_family = infile("fe_family", std::string("LAGRANGE"));
    FEFamily fefamily = Utility::string_to_enum<FEFamily>(fe_family);
    
    _vars.resize(dim*2);
    
    if (dim >= 1)
    {
        _vars[0]  = this->add_variable ( "u", static_cast<Order>(o), fefamily);
        _vars[dim]  = this->add_variable ( "tx", static_cast<Order>(o), fefamily);

        if (dim >= 2)
        {
            _vars[1]  = this->add_variable ( "v", static_cast<Order>(o), fefamily);
            _vars[dim+1]  = this->add_variable ( "ty", static_cast<Order>(o), fefamily);

            if (dim >= 3)
            {
                _vars[2]  = this->add_variable ( "w", static_cast<Order>(o), fefamily);
                _vars[dim+2]  = this->add_variable ( "tz", static_cast<Order>(o), fefamily);
            }
        }
        
    }

    // all variables are assumed time varying
    for (unsigned int i_var=0; i_var<dim*2; i_var++)
        this->time_evolving(_vars[i_var]);

    
    // Useful debugging options
    this->verify_analytic_jacobians = infile("verify_analytic_jacobians", 0.);
    this->print_jacobians = infile("print_jacobians", false);
    this->print_element_jacobians = infile("print_element_jacobians", false);
    
}



void StructuralSystem::init_context(DiffContext &context)
{
    FEMContext &c = libmesh_cast_ref<FEMContext&>(context);
    
    FEBase* elem_fe;
    
    const unsigned int dim = this->get_mesh().mesh_dimension();

    for (unsigned int i=0; i<dim+2; i++)
    {
        c.get_element_fe( _vars[i], elem_fe);
        elem_fe->get_JxW();
        elem_fe->get_phi();
        elem_fe->get_dphi();
        elem_fe->get_xyz();
    }
    
    for (unsigned int i=0; i<dim+2; i++)
    {
        c.get_side_fe( _vars[i], elem_fe);
        elem_fe->get_JxW();
        elem_fe->get_phi();
        elem_fe->get_xyz();
    }
}



bool StructuralSystem::element_time_derivative (bool request_jacobian,
                                                DiffContext &context)
{
//    FEMContext &c = libmesh_cast_ref<FEMContext&>(context);
//    
//    FEBase* elem_fe;
//    
//    c.get_element_fe(_vars[0], elem_fe);
//    
//    // shape function derivatives
//    const std::vector<std::vector<RealVectorValue> >& dphi = elem_fe->get_dphi();
//    
//    const unsigned int dim = c.elem->dim(),
//    n_qpoints = c.get_element_qrule()->n_points();
//    
//    // Element Jacobian * quadrature weights for interior integration
//    const std::vector<Real> &JxW = elem_fe->get_JxW();
//    
//    // The number of local degrees of freedom in each variable
//    unsigned int n_dofs = c.get_dof_indices(Tvar).size();
//    
//    // The subvectors and submatrices we need to fill:
//    DenseMatrix<Number>& Kmat = c.get_elem_jacobian();
//    DenseVector<Number>& Fvec = c.get_elem_residual();
//    
//    
//    std::vector<FEMOperatorMatrix> dB_mat(dim);
//    DenseMatrix<Real> tmp_mat_n2n2;               tmp_mat_n2n2.resize(n_dofs, n_dofs);
//    DenseVector<Real> tmp_vec_n1, tmp_vec_n2;     tmp_vec_n1.resize(1); tmp_vec_n2.resize(n_dofs);
//    
//    
//    
//    const unsigned int n_phi = dphi.size();
//    
//    for (unsigned int qp=0; qp != n_qpoints; qp++)
//    {
//        // calculate the operator matrix
//        for ( unsigned int i_dim=0; i_dim<dim; i_dim++ )
//        {
//            for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
//                tmp_vec_n2(i_nd) = dphi[i_nd][qp](i_dim);
//            dB_mat[i_dim].reinit(1, tmp_vec_n2);
//        }
//        
//        
//        for (unsigned int i_dim=0; i_dim<dim; i_dim++)
//        {
//            dB_mat[i_dim].vector_mult(tmp_vec_n1, c.elem_solution); // dB/dx_i T
//            dB_mat[i_dim].vector_mult(tmp_vec_n2, tmp_vec_n1); // dB/dx_i^T dB/dx_i T
//            Fvec.add(-JxW[qp]*_k_thermal, tmp_vec_n2);
//        }
//        
//        if (request_jacobian && c.elem_solution_derivative)
//        {
//            for (unsigned int i_dim=0; i_dim<dim; i_dim++)
//            {
//                dB_mat[i_dim].right_multiply_transpose(tmp_mat_n2n2, dB_mat[i_dim]); // dB/dx_i^T dB/dx_i
//                Kmat.add(-JxW[qp]*_k_thermal, tmp_mat_n2n2);
//            }
//        }
//    } // end of the quadrature point qp-loop
//    
//    //    std::cout << "inside element time derivative " << std::endl;
//    //    c.elem->print_info();
//    //    std::cout << "sol: " << std::endl; c.elem_solution.print(std::cout);
//    //    std::cout << "res: " << std::endl; Fvec.print(std::cout);
//    //    if (request_jacobian && c.elem_solution_derivative)
//    //        Kmat.print(std::cout);
//    
//    return request_jacobian;
}



bool StructuralSystem::side_time_derivative (bool request_jacobian,
                                             DiffContext &context)
{
    
    //    std::cout << "inside side constraint " << std::endl;
    //    std::cout << "elem solution" << std::endl; c.elem_solution.print(std::cout);
    //    std::cout << if_inf_bc << "  " << if_wall_bc << std::endl;
    //    std::cout << "bc vec: " << std::endl; Fvec.print(std::cout);
    //    if (request_jacobian && c.elem_solution_derivative)
    //        Kmat.print(std::cout);
    
    return request_jacobian;
}





bool StructuralSystem::mass_residual (bool request_jacobian,
                                      DiffContext &context)
{
//    FEMContext &c = libmesh_cast_ref<FEMContext&>(context);
//    
//    FEBase* elem_fe;
//    
//    c.get_element_fe(_vars[0], elem_fe);
//    
//    // shape function derivatives
//    const std::vector<std::vector<Real> >& phi = elem_fe->get_phi();
//    
//    const unsigned int dim = c.elem->dim(),
//    n_qpoints = c.get_element_qrule()->n_points();
//    
//    // Element Jacobian * quadrature weights for interior integration
//    const std::vector<Real> &JxW = elem_fe->get_JxW();
//    
//    // The number of local degrees of freedom in each variable
//    unsigned int n_dofs = c.get_dof_indices(_vars[0]).size();
//    
//    // The subvectors and submatrices we need to fill:
//    DenseMatrix<Number>& Kmat = c.get_elem_jacobian();
//    DenseVector<Number>& Fvec = c.get_elem_residual();
//    
//    
//    FEMOperatorMatrix B_mat;
//    DenseMatrix<Real> tmp_mat_n2n2;               tmp_mat_n2n2.resize(n_dofs, n_dofs);
//    DenseVector<Real> tmp_vec_n1, tmp_vec_n2;     tmp_vec_n1.resize(1); tmp_vec_n2.resize(n_dofs);
//    
//    
//    
//    const unsigned int n_phi = phi.size();
//    
//    for (unsigned int qp=0; qp != n_qpoints; qp++)
//    {
//        // calculate the operator matrix
//        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
//            tmp_vec_n2(i_nd) = phi[i_nd][qp];
//        B_mat.reinit(1, tmp_vec_n2);
//        
//        B_mat.vector_mult(tmp_vec_n1, c.elem_solution); // B dT/dt
//        B_mat.vector_mult(tmp_vec_n2, tmp_vec_n1); // B^T B dT/dt
//        Fvec.add(JxW[qp], tmp_vec_n2);
//        
//        if (request_jacobian && c.elem_solution_derivative)
//        {
//            B_mat.right_multiply_transpose(tmp_mat_n2n2, B_mat); // B^T B
//            Kmat.add(JxW[qp], tmp_mat_n2n2);
//        }
//    } // end of the quadrature point qp-loop
//    
//    //    std::cout << "inside mass residual " << std::endl;
//    //    std::cout << "elem velocity" << std::endl; c.elem_solution.print(std::cout);
//    //    std::cout << "elem solution" << std::endl; c.elem_fixed_solution.print(std::cout);
//    //    std::cout << "mass vec: " << std::endl; Fvec.print(std::cout);
//    //    if (request_jacobian && c.elem_solution_derivative)
//    //        Kmat.print(std::cout);
//    
//    return request_jacobian;
}

#endif // LIBMESH_USE_COMPLEX_NUMBERS

