//
//  assembleEuler.cpp
//  FESystem
//
//  Created by Manav Bhatia on 2/21/13.
//
//


// FESystem includes
#include "euler/assembleEuler.h"

#ifndef LIBMESH_USE_COMPLEX_NUMBERS

// Basic include files
#include "libmesh/equation_systems.h"
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/parameters.h"
#include "libmesh/quadrature.h"


// Bring in everything from the libMesh namespace
using namespace libMesh;


Real euler_solution_value(const Point& p,
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



void init_euler_variables(EquationSystems& es,
                          const std::string& system_name)
{
    // It is a good idea to make sure we are initializing
    // the proper system.
    libmesh_assert_equal_to (system_name, "EulerSystem");
    
    // Get a reference to the Convection-Diffusion system object.
    EulerSystem & system = es.get_system<EulerSystem>("EulerSystem");
    
    // Project initial conditions at time 0
    libmesh_assert_equal_to(system.time, 0.);
    
    system.project_solution(euler_solution_value, NULL, es.parameters);
}


void EulerSystem::init_data ()
{
    // initialize the system communicator for use by the Euler element base class
    this->system_comm = &this->comm();
    
    this->use_fixed_solution = true;
    
    dim = this->get_mesh().mesh_dimension();
    
    vars.resize(dim+2);
 
    // initialize the fluid values
    EulerElemBase::init_data();

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
    DenseVector<Real> flux, tmp_vec1_n1, tmp_vec2_n1, tmp_vec3_n2, conservative_sol, delta_vals;
    LS_mat.resize(n1, n_dofs); B_mat.resize(dim+2, n_dofs); Ai_Bi_advection.resize(dim+2, n_dofs); A_inv_entropy.resize(dim+2, dim+2);
    A_entropy.resize(dim+2, dim+2); tmp_mat.resize(dim+2, dim+2);
    flux.resize(n1); tmp_vec1_n1.resize(n1); tmp_vec2_n1.resize(n1); tmp_vec3_n2.resize(n_dofs); conservative_sol.resize(dim+2);
    delta_vals.resize(n_qpoints);
    for (unsigned int i=0; i<dim; i++)
    {
        dB_mat[i].resize(dim+2, n_dofs);
        Ai_advection[i].resize(dim+2, dim+2);
    }
    
    
    Real diff_val=0.;

    System& delta_val_system = this->get_equation_systems().get_system<System>("DeltaValSystem");
    NumericVector<Real>& diff_val_vec = (*delta_val_system.solution.get());
    
    if (if_use_stored_delta)
    {
        diff_val = diff_val_vec.el(c.elem->dof_number(delta_val_system.number(), 0, 0));
        for (unsigned int qp=0; qp<n_qpoints; qp++)
            delta_vals(qp) = diff_val;
    }
    

    PrimitiveSolution primitive_sol;
    
    for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
        // first update the variables at the current quadrature point
        this->update_solution_at_quadrature_point(vars, qp, c, true, c.elem_solution, conservative_sol, primitive_sol, B_mat);
        this->update_jacobian_at_quadrature_point(vars, qp, c, primitive_sol, dB_mat, Ai_advection, A_entropy, A_inv_entropy );
        
        Ai_Bi_advection.zero();
        
        for (unsigned int i_dim=0; i_dim<dim; i_dim++)
        {
            tmp_mat = Ai_advection[i_dim];
            tmp_mat.right_multiply(dB_mat[i_dim]);
            Ai_Bi_advection.add(1.0, tmp_mat);
        }
                
        this->calculate_differential_operator_matrix(vars, qp, c, c.elem_solution, primitive_sol, B_mat, dB_mat, Ai_advection, Ai_Bi_advection, A_inv_entropy, LS_mat, diff_val);

        if (if_use_stored_delta)
            diff_val = delta_vals(qp);
        else
            delta_vals(qp) = diff_val;
        
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
    
    if (!if_use_stored_delta)
    {
        diff_val = 0.;
        for (unsigned int qp=0; qp<n_qpoints; qp++)
            diff_val += delta_vals(qp);
        diff_val /= (1.*n_qpoints);
        diff_val_vec.set(c.elem->dof_number(delta_val_system.number(), 0, 0), diff_val);
    }
    
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


    
    DenseVector<Real> tmp_vec1, flux, U_vec_interpolated, tmp_vec1_n2, conservative_sol;
    DenseMatrix<Real>  eig_val, l_eig_vec, l_eig_vec_inv_tr, tmp_mat1, tmp_mat2, B_mat, A_mat;

    conservative_sol.resize(dim+2);
    tmp_vec1.resize(n_dofs); flux.resize(n1); tmp_vec1_n2.resize(n_dofs); U_vec_interpolated.resize(n1);
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
            this->update_solution_at_quadrature_point(vars, qp, c, false, c.elem_solution, conservative_sol, p_sol, B_mat);
            
            flux.zero();
            
            switch (dim)
            {
                case 3:
                    flux(3) = p_sol.u3*p_sol.rho*xini+p_sol.p*face_normals[qp](2);
                case 2:
                    flux(2) = p_sol.u2*p_sol.rho*xini+p_sol.p*face_normals[qp](1);
                case 1:
                    flux(0) = p_sol.rho*xini;
                    flux(1) = p_sol.u1*p_sol.rho*xini+p_sol.p*face_normals[qp](0);
                    flux(n1-1) = xini*(p_sol.rho*(cv*p_sol.T+p_sol.k)+p_sol.p);
                    break;
            }
            
            B_mat.vector_mult_transpose(tmp_vec1, flux);
            Fvec.add(-JxW[qp], tmp_vec1);
            
            // update the force vector
            for (unsigned int i=0; i<dim; i++)
                (*integrated_force)(i) += face_normals[qp](i)*JxW[qp]*p_sol.p;
            
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
            this->update_solution_at_quadrature_point(vars, qp, c, false, c.elem_solution, conservative_sol, p_sol, B_mat);

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
        this->update_solution_at_quadrature_point(vars, qp, c, true, c.elem_fixed_solution, conservative_sol, primitive_sol, B_mat);
        this->update_jacobian_at_quadrature_point(vars, qp, c, primitive_sol, dB_mat, Ai_advection, A_entropy, A_inv_entropy );
        
        Ai_Bi_advection.zero();
        
        for (unsigned int i_dim=0; i_dim<dim; i_dim++)
        {
            tmp_mat = Ai_advection[i_dim];
            tmp_mat.right_multiply(dB_mat[i_dim]);
            Ai_Bi_advection.add(1.0, tmp_mat);
        }
        
        this->calculate_differential_operator_matrix(vars, qp, c, c.elem_fixed_solution, primitive_sol, B_mat, dB_mat, Ai_advection, Ai_Bi_advection, A_inv_entropy, LS_mat, diff_val);
        
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
            tmp_mat.right_multiply(B_mat); // LS^T tau Bmat
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



#endif // LIBMESH_USE_COMPLEX_NUMBERS
