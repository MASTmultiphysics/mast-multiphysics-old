//
//  assembleEuler.cpp
//  FESystem
//
//  Created by Manav Bhatia on 2/21/13.
//
//

// C++ includes
#include <iomanip>

// Basic include files
#include "libmesh/equation_systems.h"
#include "libmesh/error_vector.h"
#include "libmesh/getpot.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/vtk_io.h"
#include "libmesh/tecplot_io.h"
#include "libmesh/kelly_error_estimator.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/uniform_refinement_estimator.h"
#include "libmesh/newton_solver.h"
#include "libmesh/euler_solver.h"
#include "libmesh/steady_solver.h"
#include "libmesh/getpot.h"
#include "libmesh/boundary_info.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fem_context.h"
#include "libmesh/mesh.h"
#include "libmesh/quadrature.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/zero_function.h"

// FESystem includes
#include "euler/assembleEuler.h"
#include "Base/FESystemExceptions.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;


Real solution_value(const Point& p,
                    const Parameters& parameters,
                    const std::string& sys_name,
                    const std::string& var_name)
{
    libmesh_assert_equal_to (sys_name, "EulerSystem");
    
    // since we are initializing the solution, variable values at all points is same
    
    if (var_name == "rho")
    {
        return parameters.get<Real> ("rho_inf");
    }
    else if (var_name == "rhou")
    {
        return parameters.get<Real> ("rhou_inf");
    }
    else if (var_name == "rhov")
    {
        return parameters.get<Real> ("rhov_inf");
    }
    else if (var_name == "rhow")
    {
        return parameters.get<Real> ("rhow_inf");
    }
    else if (var_name == "rhoe")
    {
        return parameters.get<Real> ("rhoe_inf");
    }
    else
        libmesh_assert(false);
}



void init_variables(EquationSystems& es,
                    const std::string& system_name)
{
    // It is a good idea to make sure we are initializing
    // the proper system.
    libmesh_assert_equal_to (system_name, "EulerSystem");
    
    // Get a reference to the Convection-Diffusion system object.
    EulerSystem & system = es.get_system<EulerSystem>("EulerSystem");
    
    // Project initial conditions at time 0
    libmesh_assert_equal_to(system.time, 0.);
    
    system.project_solution(solution_value, NULL, es.parameters);
}




void EulerSystem::init_data ()
{
    // initialize the fluid values
    EulerElemBase::init_data();

    this->use_fixed_solution = true;
    
    dim = this->get_mesh().mesh_dimension();
    const Real pi = acos(-1.);
    
    vars.resize(dim+2);
    

    // Check the input file for Reynolds number, application type,
    // approximation type
    
    GetPot infile("euler.in");

    // set parameter values
    Parameters& params = this->get_equation_systems().parameters;
    
    unsigned int o = infile("fe_order", 1);
    std::string fe_family = infile("fe_family", std::string("LAGRANGE"));
    FEFamily fefamily = Utility::string_to_enum<FEFamily>(fe_family);
    
    vars[0]  = this->add_variable ( "rho", static_cast<Order>(o), fefamily);
    this->time_evolving(vars[0]);
    params.set<Real> ("rho_inf") = rho_inf;

    vars[1] = this->add_variable ("rhou", static_cast<Order>(o), fefamily);
    this->time_evolving(vars[1]);
    params.set<Real> ("rhou_inf") = rho_inf*u1_inf;
    
    if (dim > 1)
    {
        vars[2] = this->add_variable ("rhov", static_cast<Order>(o), fefamily);
        this->time_evolving(vars[2]);
        params.set<Real> ("rhov_inf") = rho_inf*u2_inf;
    }
    
    if (dim > 2)
    {
        vars[3] = this->add_variable ("rhow", static_cast<Order>(o), fefamily);
        this->time_evolving(vars[3]);
        params.set<Real> ("rhow_inf") = rho_inf*u3_inf;
    }
    
    vars[dim+2-1] = this->add_variable ("rhoe", static_cast<Order>(o), fefamily);
    this->time_evolving(vars[dim+2-1]);
    params.set<Real> ("rhoe_inf") = rho_inf*temp_inf*cv + q0_inf;
    
    // Useful debugging options
    // Set verify_analytic_jacobians to 1e-6 to use
    this->verify_analytic_jacobians = infile("verify_analytic_jacobians", 0.);
    this->print_jacobians = infile("print_jacobians", false);
    this->print_element_jacobians = infile("print_element_jacobians", false);
    
    // Do the parent's initialization after variables and boundary constraints are defined
    FEMSystem::init_data();
}



void EulerSystem::init_context(DiffContext &context)
{
    FEMContext &c = libmesh_cast_ref<FEMContext&>(context);
    
    std::vector<FEBase*> elem_fe(dim+2);
    
    for (unsigned int i=0; i<dim+2; i++)
    {
        c.get_element_fe( vars[i], elem_fe[i]);
        elem_fe[i]->get_JxW();
        elem_fe[i]->get_phi();
        elem_fe[i]->get_dphi();
        elem_fe[i]->get_xyz();
    }
    
    std::vector<FEBase*> elem_side_fe(dim+2);
    
    for (unsigned int i=0; i<dim+2; i++)
    {
        c.get_element_fe( vars[i], elem_side_fe[i]);
        elem_side_fe[i]->get_JxW();
        elem_side_fe[i]->get_phi();
        elem_fe[i]->get_xyz();
    }
    
    
}



bool EulerSystem::element_time_derivative (bool request_jacobian,
                                           DiffContext &context)
{
    FEMContext &c = libmesh_cast_ref<FEMContext&>(context);
    
    FEBase* elem_fe;

    c.get_element_fe(vars[0], elem_fe);
    
    // Element Jacobian * quadrature weights for interior integration
    const std::vector<Real> &JxW = elem_fe->get_JxW();
        
    // The number of local degrees of freedom in each variable
    unsigned int n_dofs = 0;
    for (unsigned int i=0; i<dim+2; i++)
        n_dofs += c.get_dof_indices( vars[i] ).size();

    libmesh_assert_equal_to (n_dofs, (dim+2)*c.get_dof_indices( vars[0] ).size());
    
    // The subvectors and submatrices we need to fill:
    DenseMatrix<Number>& Kmat = c.get_elem_jacobian();
    DenseVector<Number>& Fvec = c.get_elem_residual();
    
    // Now we will build the element Jacobian and residual.
    // Constructing the residual requires the solution and its
    // gradient from the previous timestep.  This must be
    // calculated at each quadrature point by summing the
    // solution degree-of-freedom values by the appropriate
    // weight functions.
    const unsigned int n_qpoints = (c.get_element_qrule())->n_points(), n1 = dim+2;

    std::vector<DenseMatrix<Real> > dB_mat(dim), Ai_advection(dim);
    DenseMatrix<Real> LS_mat, B_mat, Ai_Bi_advection, A_entropy, A_inv_entropy, tmp_mat, tmp_mat2;
    DenseVector<Real> flux, tmp_vec1_n1, tmp_vec2_n1, tmp_vec3_n2, conservative_sol;
    LS_mat.resize(n1, n_dofs); B_mat.resize(dim+2, n_dofs); Ai_Bi_advection.resize(dim+2, n_dofs); A_inv_entropy.resize(dim+2, dim+2);
    A_entropy.resize(dim+2, dim+2); tmp_mat.resize(dim+2, dim+2);
    flux.resize(n1); tmp_vec1_n1.resize(n1); tmp_vec2_n1.resize(n1); tmp_vec3_n2.resize(n_dofs); conservative_sol.resize(dim+2);
    for (unsigned int i=0; i<dim; i++)
    {
        dB_mat[i].resize(dim+2, n_dofs);
        Ai_advection[i].resize(dim+2, dim+2);
    }
    
    
    Real diff_val=0.;

    PrimitiveSolution primitive_sol;
    
    for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
        // first update the variables at the current quadrature point
        this->update_solution_at_quadrature_point(vars, qp, c, true, true, conservative_sol, primitive_sol, B_mat);
        this->update_jacobian_at_quadrature_point(vars, qp, c, primitive_sol, dB_mat, Ai_advection, A_entropy, A_inv_entropy );
        
        Ai_Bi_advection.zero();
        
        for (unsigned int i_dim=0; i_dim<dim; i_dim++)
        {
            tmp_mat = Ai_advection[i_dim];
            tmp_mat.right_multiply(dB_mat[i_dim]);
            Ai_Bi_advection.add(1.0, tmp_mat);
        }
                
        this->calculate_differential_dperator_matrix(vars, qp, c, true, primitive_sol, B_mat, dB_mat, Ai_advection, Ai_Bi_advection, A_inv_entropy, LS_mat, diff_val);
        
//        if (this->if_update_discont_values)
//            (*this->discontinuity_capturing_value)[i] = diff_val;
//        else
//            diff_val = (*this->discontinuity_capturing_value)[i];
        
        for (unsigned int i_dim=0; i_dim<dim; i_dim++)
        {
            // Galerkin contribution from the advection flux terms
            this->calculate_advection_flux(i_dim, primitive_sol, flux); // F^adv_i
            dB_mat[i_dim].vector_mult_transpose(tmp_vec3_n2, flux); // dBw/dx_i F^adv_i
            Fvec.add(JxW[qp], tmp_vec3_n2);
            
            // discontinuity capturing operator
            dB_mat[i_dim].vector_mult(flux, c.elem_solution);
            dB_mat[i_dim].vector_mult_transpose(tmp_vec3_n2, flux);
            Fvec.add(-JxW[qp]*diff_val, tmp_vec3_n2);
        }
        
        // Least square contribution from flux
        Ai_Bi_advection.vector_mult(flux, c.elem_solution); // d F^adv_i / dxi
        LS_mat.vector_mult_transpose(tmp_vec3_n2, flux); // LS^T tau F^adv_i
        Fvec.add(-JxW[qp], tmp_vec3_n2);
        
        
        if (request_jacobian && c.elem_solution_derivative)
        {
            for (unsigned int i_dim=0; i_dim<dim; i_dim++)
            {
                // Galerkin contribution from the advection flux terms
                tmp_mat = Ai_advection[i_dim];
                tmp_mat.right_multiply(B_mat);
                dB_mat[i_dim].get_transpose(tmp_mat2);
                tmp_mat2.right_multiply(tmp_mat); // dBw/dx_i^T  dF^adv_i/ dU
                Kmat.add(JxW[qp], tmp_mat2);
                
                // discontinuity capturing term
                dB_mat[i_dim].get_transpose(tmp_mat);
                tmp_mat.right_multiply(dB_mat[i_dim]);
                Kmat.add(-JxW[qp]*diff_val, tmp_mat);
            }
            
            // Lease square contribution of flux gradient
            LS_mat.get_transpose(tmp_mat);
            tmp_mat.right_multiply(Ai_Bi_advection); // LS^T tau d^2F^adv_i / dx dU
            Kmat.add(-JxW[qp], tmp_mat);
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



bool EulerSystem::side_time_derivative (bool request_jacobian,
                                        DiffContext &context)
{
    FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

    // check for the boundary tags to check if the boundary condition needs to be applied to this element
    
    bool if_wall_bc = false, if_inf_bc = false;
    
    if ( this->get_mesh().boundary_info->has_boundary_id(c.elem, c.side, 0) ) // wall bc: both in 2D and 3D, this is the solid wall
        if_wall_bc = true;
    if ( this->get_mesh().boundary_info->has_boundary_id(c.elem, c.side, 1) ||
         this->get_mesh().boundary_info->has_boundary_id(c.elem, c.side, 2) ||
         this->get_mesh().boundary_info->has_boundary_id(c.elem, c.side, 3) ||
         this->get_mesh().boundary_info->has_boundary_id(c.elem, c.side, 4) ||
         this->get_mesh().boundary_info->has_boundary_id(c.elem, c.side, 5) ) // infinite bc (1, 2, 3) for 2D and (1, 2, 3, 4, 5) for 3D
        if_inf_bc = true;
    
    if ( !if_wall_bc && !if_inf_bc )
        return request_jacobian;
    else
        libmesh_assert_equal_to(if_wall_bc, !if_inf_bc); // make sure that only one is active
    
    
    const unsigned int n1 = dim+2;
    
    // The number of local degrees of freedom for this element
    unsigned int n_dofs = 0;
    for (unsigned int i=0; i<dim+2; i++)
        n_dofs += c.get_dof_indices( vars[i] ).size();
    
    libmesh_assert_equal_to (n_dofs, (dim+2)*c.get_dof_indices( vars[0] ).size());

    FEBase * side_fe;
    c.get_side_fe(vars[0], side_fe);
    
    // Element Jacobian * quadrature weights for interior integration
    const std::vector<Real> &JxW = side_fe->get_JxW();
    
    // Physical location of the quadrature points
    const std::vector<Point>& qpoint = side_fe->get_xyz();

    // boundary normals
    const std::vector<Point>& face_normals = side_fe->get_normals();


    Point vel; vel.zero(); // zero surface velocity
    
    
    DenseMatrix<Number>& Kmat = c.get_elem_jacobian();
    DenseVector<Number>& Fvec = c.get_elem_residual();


    
    DenseVector<Real>  normal, normal_local, tmp_vec1, flux, U_vec_interpolated, tmp_vec1_n2, conservative_sol;
    DenseMatrix<Real>  eig_val, l_eig_vec, l_eig_vec_inv_tr, tmp_mat1, tmp_mat2, B_mat, A_mat;

    conservative_sol.resize(dim+2);
    normal.resize(3); normal_local.resize(dim); tmp_vec1.resize(n_dofs); flux.resize(n1); tmp_vec1_n2.resize(n_dofs); U_vec_interpolated.resize(n1);
    eig_val.resize(n1, n1); l_eig_vec.resize(n1, n1); l_eig_vec_inv_tr.resize(n1, n1); tmp_mat1.resize(n1, n1); tmp_mat2.resize(n1, n1);
    B_mat.resize(dim+2, n_dofs); A_mat.resize(dim+2, dim+2);

    
    PrimitiveSolution p_sol;
    
    if ( if_wall_bc )
    {
        Real xini = 0.;

        for (unsigned int qp=0; qp<qpoint.size(); qp++)
        {
            xini = face_normals[qp] * vel;
            
            // first update the variables at the current quadrature point
            this->update_solution_at_quadrature_point(vars, qp, c, true, false, conservative_sol, p_sol, B_mat);
            
            flux.zero();
            
            switch (dim)
            {
                case 3:
                    flux(3) = p_sol.u3*p_sol.rho*xini+p_sol.p*face_normals[qp](2);
                case 2:
                    flux(2) = p_sol.u1*p_sol.rho*xini+p_sol.p*face_normals[qp](1);
                case 1:
                    flux(0) = p_sol.rho*xini;
                    flux(1) = p_sol.u3*p_sol.rho*xini+p_sol.p*face_normals[qp](0);
                    flux(n1-1) = xini*(p_sol.rho*(cv*p_sol.T+p_sol.k)+p_sol.p);
                    break;
            }
            
            B_mat.vector_mult_transpose(tmp_vec1, flux);
            Fvec.add(-JxW[qp], tmp_vec1);
            
            if ( request_jacobian && c.get_elem_solution_derivative() )
            {
                this->calculate_advection_flux_jacobian_for_moving_solid_wall_boundary(p_sol, xini, face_normals[qp], A_mat);
                
                tmp_mat2 = A_mat;
                tmp_mat2.right_multiply(B_mat);
                B_mat.get_transpose(tmp_mat1);
                tmp_mat1.right_multiply(tmp_mat2);
                Kmat.add(-JxW[qp], tmp_mat1);
            }
        }
    }

    
    if ( if_inf_bc )
    {
        for (unsigned int qp=0; qp<qpoint.size(); qp++)
        {
            // first update the variables at the current quadrature point
            this->update_solution_at_quadrature_point(vars, qp, c, true, false, conservative_sol, p_sol, B_mat);

            this->calculate_advection_left_eigenvector_and_inverse_for_normal(p_sol, face_normals[qp], eig_val, l_eig_vec, l_eig_vec_inv_tr);
            
            // for all eigenalues that are less than 0, the characteristics are coming into the domain, hence,
            // evaluate them using the given solution.
            tmp_mat1 = l_eig_vec_inv_tr;
            for (unsigned int j=0; j<n1; j++)
                if (eig_val(j, j) < 0)
                    tmp_mat1.scale_column(j, eig_val(j, j)); // L^-T [omaga]_{-}
                else
                    tmp_mat1.scale_column(j, 0.0);
            
            tmp_mat1.right_multiply_transpose(l_eig_vec); // A_{-} = L^-T [omaga]_{-} L^T
            this->get_infinity_vars( U_vec_interpolated );
            tmp_mat1.vector_mult(flux, U_vec_interpolated);  // f_{-} = A_{-} B U
            
            B_mat.vector_mult_transpose(tmp_vec1_n2, flux); // B^T f_{-}   (this is flux coming into the solution domain)
            Fvec.add(-JxW[qp], tmp_vec1_n2);
            
            // now calculate the flux for eigenvalues greater than 0, the characteristics go out of the domain, so that
            // the flux is evaluated using the local solution
            tmp_mat1 = l_eig_vec_inv_tr;
            for (unsigned int j=0; j<n1; j++)
                if (eig_val(j, j) > 0)
                    tmp_mat1.scale_column(j, eig_val(j, j)); // L^-T [omaga]_{+}
                else
                    tmp_mat1.scale_column(j, 0.0);

            tmp_mat1.right_multiply_transpose(l_eig_vec); // A_{+} = L^-T [omaga]_{+} L^T
            tmp_mat1.vector_mult(flux, conservative_sol); // f_{+} = A_{+} B U
            
            B_mat.vector_mult_transpose(tmp_vec1_n2, flux); // B^T f_{+}   (this is flux going out of the solution domain)
            Fvec.add(-JxW[qp], tmp_vec1_n2);

            
            if ( request_jacobian && c.get_elem_solution_derivative() )
            {
                // terms with negative eigenvalues do not contribute to the Jacobian
                
                // now calculate the Jacobian for eigenvalues greater than 0, the characteristics go out of the domain, so that
                // the flux is evaluated using the local solution
                tmp_mat1 = l_eig_vec_inv_tr;
                for (unsigned int j=0; j<n1; j++)
                    if (eig_val(j, j) > 0)
                        tmp_mat1.scale_column(j, eig_val(j, j)); // L^-T [omaga]_{+}
                    else
                        tmp_mat1.scale_column(j, 0.0);
                tmp_mat1.right_multiply_transpose(l_eig_vec); // A_{+} = L^-T [omaga]_{+} L^T
                tmp_mat1.right_multiply(B_mat);
                B_mat.get_transpose(tmp_mat2);
                tmp_mat2.right_multiply(tmp_mat1); // B^T A_{+} B   (this is flux going out of the solution domain)
                
                Kmat.add(-JxW[qp], tmp_mat2);
            }
        }
        
    }
    
//    std::cout << "inside side constraint " << std::endl;
//    std::cout << "elem solution" << std::endl; c.elem_solution.print(std::cout);
//    std::cout << if_inf_bc << "  " << if_wall_bc << std::endl;
//    std::cout << "bc vec: " << std::endl; Fvec.print(std::cout);
//    if (request_jacobian && c.elem_solution_derivative)
//        Kmat.print(std::cout);

    return request_jacobian;
}





bool EulerSystem::mass_residual (bool request_jacobian,
                                 DiffContext &context)
{
    FEMContext &c = libmesh_cast_ref<FEMContext&>(context);
    
    FEBase* elem_fe;
    
    c.get_element_fe(vars[0], elem_fe);
    
    
    // Element Jacobian * quadrature weights for interior integration
    const std::vector<Real> &JxW = elem_fe->get_JxW();
        
    // The number of local degrees of freedom in each variable
    unsigned int n_dofs = 0;
    for (unsigned int i=0; i<dim+2; i++)
        n_dofs += c.get_dof_indices( vars[i] ).size();
    
    libmesh_assert_equal_to (n_dofs, (dim+2)*c.get_dof_indices( vars[0] ).size());
    
    // The subvectors and submatrices we need to fill:
    DenseMatrix<Number>& Kmat = c.get_elem_jacobian();
    DenseVector<Number>& Fvec = c.get_elem_residual();
    
    // Now we will build the element Jacobian and residual.
    // Constructing the residual requires the solution and its
    // gradient from the previous timestep.  This must be
    // calculated at each quadrature point by summing the
    // solution degree-of-freedom values by the appropriate
    // weight functions.
    const unsigned int n_qpoints = (c.get_element_qrule())->n_points(), n1 = dim+2;
    
    std::vector<DenseMatrix<Real> > dB_mat(dim), Ai_advection(dim);
    DenseMatrix<Real> LS_mat, B_mat, Ai_Bi_advection, A_entropy, A_inv_entropy, tmp_mat;
    DenseVector<Real> flux, tmp_vec1_n1, tmp_vec3_n2, conservative_sol;
    LS_mat.resize(n1, n_dofs); B_mat.resize(dim+2, n_dofs); Ai_Bi_advection.resize(dim+2, n_dofs); A_inv_entropy.resize(dim+2, dim+2);
    A_entropy.resize(dim+2, dim+2); tmp_mat.resize(dim+2, dim+2);
    flux.resize(n1); tmp_vec1_n1.resize(n1); tmp_vec3_n2.resize(n_dofs); conservative_sol.resize(dim+2);
    for (unsigned int i=0; i<dim; i++)
    {
        dB_mat[i].resize(dim+2, n_dofs);
        Ai_advection[i].resize(dim+2, dim+2);
    }
    
    
    Real diff_val=0.;
    
    PrimitiveSolution primitive_sol;
        
    for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
        // first update the variables at the current quadrature point
        this->update_solution_at_quadrature_point(vars, qp, c, false, true, conservative_sol, primitive_sol, B_mat);
        this->update_jacobian_at_quadrature_point(vars, qp, c, primitive_sol, dB_mat, Ai_advection, A_entropy, A_inv_entropy );
        
        Ai_Bi_advection.zero();
        
        for (unsigned int i_dim=0; i_dim<dim; i_dim++)
        {
            tmp_mat = Ai_advection[i_dim];
            tmp_mat.right_multiply(dB_mat[i_dim]);
            Ai_Bi_advection.add(1.0, tmp_mat);
        }
        
        this->calculate_differential_dperator_matrix(vars, qp, c, false, primitive_sol, B_mat, dB_mat, Ai_advection, Ai_Bi_advection, A_inv_entropy, LS_mat, diff_val);
        
        //        if (this->if_update_discont_values)
        //            (*this->discontinuity_capturing_value)[i] = diff_val;
        //        else
        //            diff_val = (*this->discontinuity_capturing_value)[i];
        
        // Galerkin contribution to solution
        B_mat.vector_mult( tmp_vec1_n1, c.elem_solution );
        B_mat.vector_mult_transpose(tmp_vec3_n2, tmp_vec1_n1);
        Fvec.add(JxW[qp], tmp_vec3_n2);
        
        // Galerkin contribution to solution
        B_mat.vector_mult( tmp_vec1_n1, c.elem_solution );
        LS_mat.vector_mult_transpose(tmp_vec3_n2, tmp_vec1_n1);
        Fvec.add(JxW[qp], tmp_vec3_n2);
        
        
        if (request_jacobian && c.get_elem_solution_derivative())
        {
            // contribution from unsteady term
            // Galerkin contribution of velocity
            B_mat.get_transpose(tmp_mat);
            tmp_mat.right_multiply(B_mat);
            Kmat.add(JxW[qp], tmp_mat);
            
            // LS contribution of velocity
            LS_mat.get_transpose(tmp_mat);
            tmp_mat.right_multiply(B_mat); // LS^T tau A Bmat
            Kmat.add(JxW[qp], tmp_mat);
        }
    } // end of the quadrature point qp-loop
    
//    std::cout << "inside mass residual " << std::endl;
//    std::cout << "elem velocity" << std::endl; c.elem_solution.print(std::cout);
//    std::cout << "elem solution" << std::endl; c.elem_fixed_solution.print(std::cout);
//    std::cout << "mass vec: " << std::endl; Fvec.print(std::cout);
//    if (request_jacobian && c.elem_solution_derivative)
//        Kmat.print(std::cout);

    return request_jacobian;
}




void EulerSystem::postprocess()
{
    
}



// The main program.
int libmesh_euler_analysis (int argc, char* const argv[])
{
    // Initialize libMesh.
    LibMeshInit init (argc, argv);
    
    // Parse the input file
    GetPot infile("system_input.in");
    
    // Read in parameters from the input file
    const Real global_tolerance          = infile("global_tolerance", 0.);
    const unsigned int nelem_target      = infile("n_elements", 400);
    const bool transient                 = infile("transient", true);
    const Real deltat                    = infile("deltat", 0.005);
    unsigned int n_timesteps             = infile("n_timesteps", 20);
    const unsigned int write_interval    = infile("write_interval", 5);
    const unsigned int coarsegridsize    = infile("coarsegridsize", 15);
    const unsigned int coarserefinements = infile("coarserefinements", 0);
    const unsigned int max_adaptivesteps = infile("max_adaptivesteps", 0);
    const unsigned int dim               = infile("dimension", 2);
    
    // Skip higher-dimensional examples on a lower-dimensional libMesh build
    libmesh_example_assert(dim <= LIBMESH_DIM, "2D/3D support");
    
    // We have only defined 2 and 3 dimensional problems
    libmesh_assert (dim == 2 || dim == 3);
    
    // Create a mesh.
    Mesh mesh;
    
    // And an object to refine it
    MeshRefinement mesh_refinement(mesh);
    mesh_refinement.coarsen_by_parents() = true;
    mesh_refinement.absolute_global_tolerance() = global_tolerance;
    mesh_refinement.nelem_target() = nelem_target;
    mesh_refinement.refine_fraction() = 0.3;
    mesh_refinement.coarsen_fraction() = 0.3;
    mesh_refinement.coarsen_threshold() = 0.1;
    
    const Real pi = acos(-1.), x_length=4.5, y_length=4.5, z_length=1.0, t_by_c = 0.05, chord = 1.0, span = 1.0, thickness = 0.5*t_by_c*chord,
    x0=x_length*0.5-chord*0.5, x1=x0+chord, y0=y_length*0.5-span*0.5, y1=y0+span ;

    // Use the MeshTools::Generation mesh generator to create a uniform
    // grid on the square [-1,1]^D.  We instruct the mesh generator
    // to build a mesh of 8x8 \p Quad9 elements in 2D, or \p Hex27
    // elements in 3D.  Building these higher-order elements allows
    // us to use higher-order approximation, as in example 3.
    if (dim == 2)
        MeshTools::Generation::build_square (mesh,
                                             coarsegridsize,
                                             coarsegridsize,
                                             0., x_length,
                                             0., y_length,
                                             QUAD4);
    else if (dim == 3)
        MeshTools::Generation::build_cube (mesh,
                                           coarsegridsize,
                                           coarsegridsize,
                                           coarsegridsize,
                                           0., x_length,
                                           0., y_length,
                                           0., z_length,
                                           HEX8);
    
    
    mesh_refinement.uniformly_refine(coarserefinements);
    
    // Print information about the mesh to the screen.
    mesh.print_info();

    
    MeshBase::node_iterator   n_it  = mesh.nodes_begin();
    const Mesh::node_iterator n_end = mesh.nodes_end();
    
    Real x_val, y_val, z_val;
    
    for (; n_it != n_end; n_it++)
    {
        Node& n =  **n_it;
        if (dim == 2)
            if ((n(0) >= x0) && (n(0) <= x1))
            {
                x_val = n(0);
                y_val = n(1);
                
                n(1) += thickness*(1.0-y_val/y_length)*sin(pi*(x_val-x0)/chord);
            }
        
        if (dim == 3)
            if ((n(0) >= x0) && (n(0) <= x1) &&
                (n(1) >= y0) && (n(1) <= y1))
            {
                x_val = n(0);
                y_val = n(1);
                z_val = n(2);
                
                n(2) += thickness*(1.0-z_val/z_length)*sin(pi*(x_val-x0)/chord)*sin(pi*(y_val-y0)/span);
            }
    }

    // Create an equation systems object.
    EquationSystems equation_systems (mesh);
    
    // Declare the system "Navier-Stokes" and its variables.
    EulerSystem & system =
    equation_systems.add_system<EulerSystem> ("EulerSystem");
    
    // Solve this as a time-dependent or steady system
    if (transient)
        system.time_solver =
        AutoPtr<TimeSolver>(new EulerSolver(system));
    else
    {
        system.time_solver =
        AutoPtr<TimeSolver>(new SteadySolver(system));
        libmesh_assert_equal_to (n_timesteps, 1);
    }
    
    system.attach_init_function (init_variables);
    
    // Initialize the system
    equation_systems.init ();

    system.print_residual_norms = infile("print_residual_norms", false);
    system.print_residuals = infile("print_residuals", false);
    system.print_jacobian_norms = infile("print_jacobian_norms", false);
    system.print_jacobians = infile("print_jacobians", false);

    // Set the time stepping options
    system.deltat = deltat;
    
    // And the nonlinear solver options
    NewtonSolver &solver = dynamic_cast<NewtonSolver&>(*(system.time_solver->diff_solver().get()));
    solver.quiet = infile("solver_quiet", true);
    solver.verbose = !solver.quiet;
    solver.max_nonlinear_iterations =
    infile("max_nonlinear_iterations", 15);
    solver.relative_step_tolerance =
    infile("relative_step_tolerance", 1.e-3);
    solver.relative_residual_tolerance =
    infile("relative_residual_tolerance", 0.0);
    solver.absolute_residual_tolerance =
    infile("absolute_residual_tolerance", 0.0);
    solver.continue_after_backtrack_failure =
    infile("continue_after_backtrack_failure", false);
    solver.continue_after_max_iterations =
    infile("continue_after_max_iterations", false);
    solver.require_residual_reduction =
    infile("require_residual_reduction", true);
    
    // And the linear solver options
    solver.max_linear_iterations =
    infile("max_linear_iterations", 50000);
    solver.initial_linear_tolerance =
    infile("initial_linear_tolerance", 1.e-3);
    //solver.brent_line_search = false;
    
    // Print information about the system to the screen.
    equation_systems.print_info();
    
    
    // Now we begin the timestep loop to compute the time-accurate
    // solution of the equations.
    for (unsigned int t_step=0; t_step != n_timesteps; ++t_step)
    {
        // A pretty update message
        std::cout << "\n\nSolving time step " << t_step << ", time = "
        << system.time << std::endl;
        
        // Adaptively solve the timestep
        unsigned int a_step = 0;
        for (; a_step != max_adaptivesteps; ++a_step)
        {
            system.solve();
            
            system.postprocess();
            
            ErrorVector error;
            
            AutoPtr<ErrorEstimator> error_estimator;
            
            // To solve to a tolerance in this problem we
            // need a better estimator than Kelly
            if (global_tolerance != 0.)
            {
                // We can't adapt to both a tolerance and a mesh
                // size at once
                libmesh_assert_equal_to (nelem_target, 0);
                
                UniformRefinementEstimator *u =
                new UniformRefinementEstimator;
                
                // The lid-driven cavity problem isn't in H1, so
                // lets estimate L2 error
                u->error_norm = L2;
                
                error_estimator.reset(u);
            }
            else
            {
                // If we aren't adapting to a tolerance we need a
                // target mesh size
                libmesh_assert_greater (nelem_target, 0);
                
                // Kelly is a lousy estimator to use for a problem
                // not in H1 - if we were doing more than a few
                // timesteps we'd need to turn off or limit the
                // maximum level of our adaptivity eventually
                error_estimator.reset(new KellyErrorEstimator);
            }
            
            // Calculate error based on u and v (and w?) but not p
            std::vector<Real> weights(dim+2,0.0);  // all set to 1.0
            weights[0] = 1.0;
            // Keep the same default norm type.
            std::vector<FEMNormType>
            norms(1, error_estimator->error_norm.type(0));
            error_estimator->error_norm = SystemNorm(norms, weights);
            
            error_estimator->estimate_error(system, error);
            
            // Print out status at each adaptive step.
            Real global_error = error.l2_norm();
            std::cout << "Adaptive step " << a_step << ": " << std::endl;
            if (global_tolerance != 0.)
                std::cout << "Global_error = " << global_error
                << std::endl;
            if (global_tolerance != 0.)
                std::cout << "Worst element error = " << error.maximum()
                << ", mean = " << error.mean() << std::endl;
            
            if (global_tolerance != 0.)
            {
                // If we've reached our desired tolerance, we
                // don't need any more adaptive steps
                if (global_error < global_tolerance)
                    break;
                mesh_refinement.flag_elements_by_error_tolerance(error);
            }
            else
            {
                // If flag_elements_by_nelem_target returns true, this
                // should be our last adaptive step.
                if (mesh_refinement.flag_elements_by_nelem_target(error))
                {
                    mesh_refinement.refine_and_coarsen_elements();
                    equation_systems.reinit();
                    a_step = max_adaptivesteps;
                    break;
                }
            }
            
            // Carry out the adaptive mesh refinement/coarsening
            mesh_refinement.refine_and_coarsen_elements();
            equation_systems.reinit();
            
            std::cout << "Refined mesh to "
            << mesh.n_active_elem()
            << " active elements and "
            << equation_systems.n_active_dofs()
            << " active dofs." << std::endl;
        }
        // Do one last solve if necessary
        if (a_step == max_adaptivesteps)
        {
            system.solve();
            
            system.postprocess();
        }
        
        // Advance to the next timestep in a transient problem
        system.time_solver->advance_timestep();
        
#ifdef LIBMESH_HAVE_EXODUS_API
        // Write out this timestep if we're requested to
        if ((t_step+1)%write_interval == 0)
        {
            std::ostringstream file_name;
            
            // We write the file in the ExodusII format.
            file_name << "out_"
            << std::setw(3)
            << std::setfill('0')
            << std::right
            << t_step+1
            << ".e";
            
            ExodusII_IO(mesh).write_timestep(file_name.str(),
                                             equation_systems,
                                             1, /* This number indicates how many time steps
                                                 are being written to the file */
                                             system.time);
        }
#endif // #ifdef LIBMESH_HAVE_EXODUS_API
    }
    
    // All done.  
    return 0;
}
