/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2014  Manav Bhatia
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

// C++ includes
#include <iomanip>

// MAST includes
#include "FluidElems/fluid_system.h"



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







Real euler_solution_value(const libMesh::Point& p,
                                   const libMesh::Parameters& parameters,
                                   const std::string& sys_name,
                                   const std::string& var_name)
{
    libmesh_assert_equal_to (sys_name, "FluidSystem");
    
    // since we are initializing the solution, variable values at all points is same
    
    if (var_name == "rho")
        
    {
        return parameters.get<Real> ("rho_inf");
    }
    else if (var_name == "rhoux")
    {
        return parameters.get<Real> ("rhoux_inf");
    }
    else if (var_name == "rhouy")
    {
        return parameters.get<Real> ("rhouy_inf");
    }
    else if (var_name == "rhouz")
    {
        return parameters.get<Real> ("rhouz_inf");
    }
    else if (var_name == "rhoe")
    {
        return parameters.get<Real> ("rhoe_inf");
    }
    else
        libmesh_assert(false);
}



void init_euler_variables(libMesh::EquationSystems& es,
                          const std::string& system_name)
{
    // It is a good idea to make sure we are initializing
    // the proper system.
    libmesh_assert_equal_to (system_name, "FluidSystem");
    
    // Get a reference to the Convection-Diffusion system object.
    FluidSystem & system = es.get_system<FluidSystem>("FluidSystem");
    
    // Project initial conditions at time 0
    libmesh_assert_equal_to(system.time, 0.);
    
    system.project_solution(euler_solution_value, NULL, es.parameters);
}


void FluidSystem::init_data ()
{
    this->use_fixed_solution = true;
    
    dim = this->get_mesh().mesh_dimension();
    
    vars.resize(dim+2);
    
    // initialize the fluid values
    FluidElemBase::init_data();
    
    // Check the input file for Reynolds number, application type,
    // approximation type
    
    // set parameter values
    libMesh::Parameters& params = this->get_equation_systems().parameters;
    
    unsigned int o = _infile("fe_order", 1);
    std::string fe_family = _infile("fe_family", std::string("LAGRANGE"));
    libMesh::FEFamily fefamily = libMesh::Utility::string_to_enum<libMesh::FEFamily>(fe_family);
    
    vars[0]  = this->add_variable ( "rho", static_cast<libMesh::Order>(o), fefamily);
    this->time_evolving(vars[0]);
    params.set<Real> ("rho_inf") = flight_condition->rho();
    
    vars[1] = this->add_variable ("rhoux", static_cast<libMesh::Order>(o), fefamily);
    this->time_evolving(vars[1]);
    params.set<Real> ("rhoux_inf") = flight_condition->rho_u1();
    
    if (dim > 1)
    {
        vars[2] = this->add_variable ("rhouy", static_cast<libMesh::Order>(o), fefamily);
        this->time_evolving(vars[2]);
        params.set<Real> ("rhouy_inf") = flight_condition->rho_u2();
    }
    
    if (dim > 2)
    {
        vars[3] = this->add_variable ("rhouz", static_cast<libMesh::Order>(o), fefamily);
        this->time_evolving(vars[3]);
        params.set<Real> ("rhouz_inf") = flight_condition->rho_u3();
    }
    
    vars[dim+2-1] = this->add_variable ("rhoe", static_cast<libMesh::Order>(o), fefamily);
    this->time_evolving(vars[dim+2-1]);
    params.set<Real> ("rhoe_inf") = flight_condition->rho_e();
    
    // Useful debugging options
    this->verify_analytic_jacobians = _infile("verify_analytic_jacobians", 0.);
    this->print_jacobians = _infile("print_jacobians", false);
    this->print_element_jacobians = _infile("print_element_jacobians", false);
    
    // initialize the Dirichlet boundary conditions
    std::multimap<unsigned, FluidBoundaryConditionType>::const_iterator
    bc_it = _boundary_condition.begin(), bc_end = _boundary_condition.end();
    
    std::set<libMesh::boundary_id_type> no_slip_boundary, isothermal_boundary;
    
    for ( ; bc_it!=bc_end; bc_it++)
    {
        switch (bc_it->second)
        {
            case NO_SLIP_WALL: // rho vi = 0
                no_slip_boundary.insert(bc_it->first);
                break;
                
            default:
                break;
        }
    }
    
    if (_if_viscous && (no_slip_boundary.size() > 0))
    {
        // Dirichlet Boundary condition for the no-slip wall
        libMesh::ZeroFunction<Real> zero_function;
        std::vector<unsigned int> rho_ui_vars(dim);
        for (unsigned int i_dim=0; i_dim<dim; i_dim++)
            rho_ui_vars[i_dim] = vars[i_dim+1];
        
        this->get_dof_map().add_dirichlet_boundary(libMesh::DirichletBoundary(no_slip_boundary, rho_ui_vars,
                                                                     &zero_function));
    }
    
    
    // Do the parent's initialization after variables and boundary constraints are defined
    libMesh::FEMSystem::init_data();
}



void FluidSystem::init_context(libMesh::DiffContext &context)
{
    context.add_localized_vector(*_dc_ref_sol, *this);
    
    libMesh::FEMContext &c = libMesh::libmesh_cast_ref<libMesh::FEMContext&>(context);
    
    std::vector<libMesh::FEBase*> elem_fe(dim+2);
    
    for (unsigned int i=0; i<dim+2; i++)
    {
        c.get_element_fe( vars[i], elem_fe[i]);
        elem_fe[i]->get_JxW();
        elem_fe[i]->get_phi();
        elem_fe[i]->get_dphi();
        elem_fe[i]->get_xyz();
        if (_if_viscous)
            elem_fe[i]->get_d2phi();
    }
    
    std::vector<libMesh::FEBase*> elem_side_fe(dim+2);
    
    for (unsigned int i=0; i<dim+2; i++)
    {
        c.get_side_fe( vars[i], elem_side_fe[i]);
        elem_side_fe[i]->get_JxW();
        elem_side_fe[i]->get_phi();
        elem_side_fe[i]->get_dphi();
        elem_side_fe[i]->get_xyz();
    }
    
    libMesh::FEMSystem::init_context(context);
}



bool FluidSystem::element_time_derivative (bool request_jacobian,
                                           libMesh::DiffContext &context)
{
    libMesh::FEMContext &c = libMesh::libmesh_cast_ref<libMesh::FEMContext&>(context);
    
    libMesh::FEBase* elem_fe;
    
    c.get_element_fe(vars[0], elem_fe);
    
    // Element Jacobian * quadrature weights for interior integration
    const std::vector<Real> &JxW = elem_fe->get_JxW();
    
    // The subvectors and submatrices we need to fill:
    DenseRealMatrix& Kmat = c.get_elem_jacobian();
    DenseRealVector& Fvec = c.get_elem_residual();
    
    // Now we will build the element Jacobian and residual.
    // Constructing the residual requires the solution and its
    // gradient from the previous timestep.  This must be
    // calculated at each quadrature point by summing the
    // solution degree-of-freedom values by the appropriate
    // weight functions.
    const unsigned int n_qpoints = c.get_element_qrule().n_points(), n1 = dim+2,
    n_dofs = n1*elem_fe->n_shape_functions();
    
    FEMOperatorMatrix B_mat;
    std::vector<FEMOperatorMatrix> dB_mat(dim);
    std::vector<DenseRealMatrix >  Ai_advection(dim);
    DenseRealMatrix LS_mat, LS_sens, Ai_Bi_advection, mat_n1n1,
    mat_n1n2, mat2_n2n2, mat3, A_sens, stress_tensor,
    dprim_dcons, dcons_dprim;
    
    DenseRealVector flux, vec1_n1, vec2_n1, vec3_n2,
    conservative_sol, temp_grad, diff_val, dc_ref_sol;
    
    LS_mat.resize(n1, n_dofs); LS_sens.resize(n_dofs, n_dofs);
    Ai_Bi_advection.resize(dim+2, n_dofs);
    mat_n1n2.resize(dim+2, n_dofs); A_sens.resize(n1, n_dofs);
    stress_tensor.resize(dim, dim); dprim_dcons.resize(dim+2, dim+2);
    dcons_dprim.resize(dim+2, dim+2);
    mat2_n2n2.resize(n_dofs, n_dofs);
    
    flux.resize(n1); vec1_n1.resize(n1); vec2_n1.resize(n1);
    vec3_n2.resize(n_dofs); conservative_sol.resize(dim+2);
    temp_grad.resize(dim); diff_val.resize(dim);
    
    for (unsigned int i=0; i<dim; i++)
        Ai_advection[i].resize(dim+2, dim+2);
    
    std::vector<std::vector<DenseRealMatrix > > flux_jacobian_sens;
    flux_jacobian_sens.resize(dim);
    for (unsigned int i_dim=0; i_dim<dim; i_dim++)
    {
        flux_jacobian_sens[i_dim].resize(n1); // number of variables for sensitivity
        for (unsigned int i_cvar=0; i_cvar<n1; i_cvar++)
            flux_jacobian_sens[i_dim][i_cvar].resize(n1, n1);
    }

    if (!if_use_stored_dc_coeff)
        dc_ref_sol = c.get_elem_solution();
    else
        dc_ref_sol = c.get_localized_vector(*_dc_ref_sol);
    
    PrimitiveSolution primitive_sol;
    const std::vector<std::vector<Real> >& phi = elem_fe->get_phi(); // assuming that all variables have the same interpolation
    const unsigned int n_phi = phi.size();
    
    for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
        // first update the variables at the current quadrature point
        this->update_solution_at_quadrature_point
        (vars, qp, c,
         true, c.get_elem_solution(),
         conservative_sol, primitive_sol,
         B_mat, dB_mat);
        
        for ( unsigned int i_dim=0; i_dim < dim; i_dim++)
            this->calculate_advection_flux_jacobian ( i_dim, primitive_sol,
                                                     Ai_advection[i_dim] );
        
        
        
        if (_if_viscous)
        {
            // stress tensor and temperature gradient
            this->calculate_conservative_variable_jacobian(primitive_sol,
                                                           dcons_dprim,
                                                           dprim_dcons);
            
            this->calculate_diffusion_tensors(c.get_elem_solution(), dB_mat,
                                              dprim_dcons, primitive_sol,
                                              stress_tensor, temp_grad);
        }
        
        Ai_Bi_advection.zero();
        
        for (unsigned int i_dim=0; i_dim<dim; i_dim++)
        {
            dB_mat[i_dim].left_multiply(mat_n1n2, Ai_advection[i_dim]);
            Ai_Bi_advection.add(1.0, mat_n1n2);
        }
        
        if (_if_full_linearization)
        {
            for (unsigned int i_dim=0; i_dim<dim; i_dim++)
                this->calculate_advection_flux_jacobian_sensitivity_for_conservative_variable
                (i_dim, primitive_sol, flux_jacobian_sens[i_dim]);
        }
        if (_if_update_stabilization_per_quadrature_point || (qp == 0)) {
            this->calculate_differential_operator_matrix
            (vars, qp, c,
             c.get_elem_solution(), primitive_sol,
             B_mat, dB_mat, Ai_advection,
             Ai_Bi_advection, flux_jacobian_sens,
             LS_mat, LS_sens);
            
            this->calculate_aliabadi_discontinuity_operator(vars, qp, c,
                                                            primitive_sol,
                                                            dc_ref_sol,
                                                            dB_mat,
                                                            Ai_Bi_advection,
                                                            diff_val);
        }
        
        for (unsigned int i_dim=0; i_dim<dim; i_dim++)
        {
            // Galerkin contribution from the advection flux terms
            this->calculate_advection_flux(i_dim, primitive_sol, flux); // F^adv_i
            dB_mat[i_dim].vector_mult_transpose(vec3_n2, flux); // dBw/dx_i F^adv_i
            Fvec.add(JxW[qp], vec3_n2);
            
            if (_if_viscous)
            {
                // Galerkin contribution from the diffusion flux terms
                this->calculate_diffusion_flux(i_dim, primitive_sol,
                                               stress_tensor, temp_grad,
                                               flux); // F^diff_i
                dB_mat[i_dim].vector_mult_transpose(vec3_n2, flux); // dBw/dx_i F^diff_i
                Fvec.add(-JxW[qp], vec3_n2);
            }
            
            // discontinuity capturing operator
            dB_mat[i_dim].vector_mult(flux, c.get_elem_solution());
            dB_mat[i_dim].vector_mult_transpose(vec3_n2, flux);
            Fvec.add(-JxW[qp]*diff_val(i_dim), vec3_n2);
        }
        
        // Least square contribution from divergence of advection flux
        Ai_Bi_advection.vector_mult(flux, c.get_elem_solution()); // d F^adv_i / dxi
        LS_mat.vector_mult_transpose(vec3_n2, flux); // LS^T tau F^adv_i
        Fvec.add(-JxW[qp], vec3_n2);
        
        // Least square contribution from divergence of diffusion flux
        // TODO: this requires a 2nd order differential of the flux
        
        
        if (request_jacobian && c.elem_solution_derivative)
        {
            A_sens.zero();
            
            for (unsigned int i_dim=0; i_dim<dim; i_dim++)
            {
                // Galerkin contribution from the advection flux terms
                B_mat.left_multiply(mat_n1n2, Ai_advection[i_dim]);
                dB_mat[i_dim].right_multiply_transpose(mat2_n2n2,
                                                       mat_n1n2);// dBw/dx_i^T  dF^adv_i/ dU
                Kmat.add(JxW[qp], mat2_n2n2);
                
                // discontinuity capturing term
                dB_mat[i_dim].right_multiply_transpose(mat2_n2n2,
                                                       dB_mat[i_dim]);
                Kmat.add(-JxW[qp]*diff_val(i_dim), mat2_n2n2);

                if (_if_full_linearization)
                {
                    // sensitivity of Ai_Bi with respect to U:   [dAi/dUj.Bi.U  ...  dAi/dUn.Bi.U]
                    dB_mat[i_dim].vector_mult(vec2_n1,
                                              c.get_elem_solution());
                    for (unsigned int i_cvar=0; i_cvar<n1; i_cvar++)
                    {
                        flux_jacobian_sens[i_dim][i_cvar].vector_mult(vec1_n1,
                                                                      vec2_n1);
                        for (unsigned int i_phi=0; i_phi<n_phi; i_phi++)
                            A_sens.add_column((n_phi*i_cvar)+i_phi,
                                              phi[i_phi][qp],
                                              vec1_n1); // assuming that all variables have same n_phi
                    }
                }
                
                if (_if_viscous)
                {
                    for (unsigned int deriv_dim=0; deriv_dim<dim; deriv_dim++)
                    {
                        this->calculate_diffusion_flux_jacobian
                        (i_dim, deriv_dim,
                         primitive_sol, mat_n1n1); // Kij
                        dB_mat[deriv_dim].left_multiply(mat_n1n2,
                                                        mat_n1n1); // Kij dB/dx_j
                        dB_mat[i_dim].right_multiply_transpose(mat2_n2n2,
                                                               mat_n1n2); // dB/dx_i^T Kij dB/dx_j
                        Kmat.add(-JxW[qp], mat2_n2n2);
                    }
                    
                }
            }
            
            // Least square contribution of flux gradient
            LS_mat.get_transpose(mat3);
            mat3.right_multiply(Ai_Bi_advection); // LS^T tau d^2F^adv_i / dx dU   (Ai constant)
            Kmat.add(-JxW[qp], mat3);
            
            if (_if_full_linearization)
            {
                LS_mat.get_transpose(mat3);
                mat3.right_multiply(A_sens); // LS^T tau d^2F^adv_i / dx dU  (Ai sensitivity)
                Kmat.add(-JxW[qp], mat3);
                
                // contribution sensitivity of the LS.tau matrix
                Kmat.add(-JxW[qp], LS_sens);
            }
        }
    } // end of the quadrature point qp-loop
    
    //    libMesh::out << "inside element time derivative " << std::endl;
    //    c.elem->print_info();
    //    libMesh::out << "sol: " << std::endl; c.elem_solution.print(libMesh::out);
    //    libMesh::out << "res: " << std::endl; Fvec.print(libMesh::out);
    //    if (request_jacobian && c.elem_solution_derivative)
    //        Kmat.print(libMesh::out);
    
    return request_jacobian;
}



bool FluidSystem::side_time_derivative (bool request_jacobian,
                                        libMesh::DiffContext &context)
{
    libMesh::FEMContext &c = libMesh::libmesh_cast_ref<libMesh::FEMContext&>(context);
    
    // check for the boundary tags to check if the boundary condition needs to be applied to this element
    
    FluidBoundaryConditionType mechanical_bc_type, thermal_bc_type;
    
    unsigned int n_mechanical_bc = 0, n_thermal_bc = 0;
    
    std::multimap<unsigned int, FluidBoundaryConditionType>::const_iterator
    bc_it     = this->_boundary_condition.begin(),
    bc_end    = this->_boundary_condition.end();
    
    
    for ( ; bc_it != bc_end; bc_it++)
    {
        if ( c.has_side_boundary_id(bc_it->first))
        {
            switch (bc_it->second)
            {
                case SLIP_WALL:
                    mechanical_bc_type = SLIP_WALL;
                    n_mechanical_bc++;
                    break;
             
                case SYMMETRY_WALL:
                    mechanical_bc_type = SYMMETRY_WALL;
                    n_mechanical_bc++;
                    break;
                    
                case NO_SLIP_WALL:
                    libmesh_assert(_if_viscous);
                    mechanical_bc_type = NO_SLIP_WALL;
                    n_mechanical_bc++;
                    break;
                    
                case FAR_FIELD:
                    mechanical_bc_type = FAR_FIELD;
                    n_mechanical_bc++;
                    break;

                case EXHAUST:
                    mechanical_bc_type = EXHAUST;
                    n_mechanical_bc++;
                    break;

                case ISOTHERMAL:
                    libmesh_assert(_if_viscous);
                    libmesh_assert_equal_to(mechanical_bc_type, NO_SLIP_WALL);
                    thermal_bc_type = ISOTHERMAL;
                    n_thermal_bc++;
                    break;
                    
                case ADIABATIC:
                    libmesh_assert(_if_viscous);
                    libmesh_assert_equal_to(mechanical_bc_type, NO_SLIP_WALL);
                    thermal_bc_type = ADIABATIC;
                    n_thermal_bc++;
                    break;
                    
            }
        }
    }
    
    // return if no boundary condition is applied
    if (n_mechanical_bc == 0)
        return request_jacobian;
    
    
    libmesh_assert_equal_to(n_mechanical_bc, 1); // make sure that only one is active
    libmesh_assert(n_thermal_bc <= 1);           // make sure that only one is active
    if (n_thermal_bc == 1)                       // thermal bc must be accompanied with mechanical bc
        libmesh_assert_equal_to(mechanical_bc_type, NO_SLIP_WALL);
    
    
    
    libMesh::FEBase * side_fe;
    c.get_side_fe(vars[0], side_fe); // assuming all variables have the same FE

    const unsigned int n1 = dim+2, n_dofs = n1*side_fe->n_shape_functions(),
    spatial_dim = this->get_mesh().spatial_dimension();
    
    
    
    // Element Jacobian * quadrature weights for interior integration
    const std::vector<Real> &JxW = side_fe->get_JxW();
    
    // Physical location of the quadrature points
    const std::vector<libMesh::Point>& qpoint = side_fe->get_xyz();
    
    // boundary normals
    const std::vector<libMesh::Point>& face_normals = side_fe->get_normals();
    
    DenseRealMatrix& Kmat = c.get_elem_jacobian();
    DenseRealVector& Fvec = c.get_elem_residual();
    
    
    
    FEMOperatorMatrix B_mat;
    DenseRealVector vec1_n2, flux, U_vec_interpolated, vec2_n1,
    conservative_sol, temp_grad, dnormal, surface_def, surface_vel, local_normal, uvec;
    DenseRealMatrix  eig_val, l_eig_vec, l_eig_vec_inv_tr, mat_n1n1,
    mat1_n1n2, mat2_n2n2, mat3_n1n1, A_mat, dcons_dprim, dprim_dcons,
    stress_tensor;
    
    conservative_sol.resize(dim+2); temp_grad.resize(dim);
    vec1_n2.resize(n_dofs); flux.resize(n1); vec2_n1.resize(n1);
    U_vec_interpolated.resize(n1); dnormal.resize(spatial_dim);
    uvec.resize(spatial_dim); surface_def.resize(spatial_dim);
    surface_vel.resize(spatial_dim); local_normal.resize(spatial_dim);
    eig_val.resize(n1, n1); l_eig_vec.resize(n1, n1);
    l_eig_vec_inv_tr.resize(n1, n1);
    mat1_n1n2.resize(n1, n_dofs); mat2_n2n2.resize(n_dofs, n_dofs);
    A_mat.resize(dim+2, dim+2); mat_n1n1.resize(n1, n1);
    dcons_dprim.resize(n1, n1); dprim_dcons.resize(n1, n1);
    stress_tensor.resize(dim, dim); mat3_n1n1.resize(n1, n1);
    
    std::vector<FEMOperatorMatrix> dB_mat(dim);
    std::vector<DenseRealMatrix > Ai_advection(dim);
    for (unsigned int i_dim=0; i_dim<dim; i_dim++)
        Ai_advection[i_dim].resize(n1, n1);
    
    PrimitiveSolution p_sol;
    
    
    switch (mechanical_bc_type)
    // adiabatic and isothermal are handled in the no-slip wall BC
    {
        case NO_SLIP_WALL: // only for viscous flows
                           // conditions enforced are
                           // vi ni = 0, vi = 0 = rho.vi      (no-slip wall, Dirichlet BC)
                           // thermal BC is handled here
        {
            
            Real ui_ni = 0.;
            
            for (unsigned int qp=0; qp<qpoint.size(); qp++)
            {
                // first update the variables at the current quadrature point
                this->update_solution_at_quadrature_point
                (vars, qp, c,
                 false, c.get_elem_solution(),
                 conservative_sol, p_sol,
                 B_mat, dB_mat);
                
                // stress tensor and temperature gradient
                this->calculate_conservative_variable_jacobian
                (p_sol, dcons_dprim, dprim_dcons);
                this->calculate_diffusion_tensors
                (c.get_elem_solution(), dB_mat, dprim_dcons,
                 p_sol, stress_tensor, temp_grad);
                
                if (thermal_bc_type == ADIABATIC)
                    temp_grad.zero();

                
                // copy the surface normal
                for (unsigned int i_dim=0; i_dim<dim; i_dim++)
                    local_normal(i_dim) = face_normals[qp](i_dim);
                
                // now check if the surface deformation is defined and
                // needs to be applied through transpiration boundary
                // condition
                ui_ni = 0.;
                p_sol.get_uvec(uvec);
                
                if (surface_motion) // get the surface motion data
                {
                    surface_motion->surface_velocity
                    (this->time, qpoint[qp], face_normals[qp],
                     surface_def, surface_vel, dnormal);
                    
                    // update the normal with the deformation
                    // this defines the normal of the surface that has been
                    // deformed, although the geometry of the flow mesh does
                    // not conform to that deformation
                    local_normal.add(1., dnormal);

                    //    ui * (ni + dni) = wi_dot * (ni + dni)
                    // => ui * ni = wi_dot * (ni + dni) - ui * dni
                    
                    // now initialize the surface velocity
                    // note that the perturbed local normal is used here
                    // since it resembles the normal of the surface that
                    // has undergone a static deformation, even though the
                    // surface_vel might be identically zero for a static
                    // body.
                    ui_ni  = surface_vel.dot(local_normal);
                    ui_ni -= uvec.dot(dnormal);
                }
                
                flux.zero();
                flux.add(ui_ni, conservative_sol);
                flux(n1-1) += p_sol.p*ui_ni;
                for (unsigned int i_dim=0; i_dim<dim; i_dim++)
                    flux(i_dim+1) += p_sol.p * face_normals[qp](i_dim);
                
                // contribution from advection flux
                B_mat.vector_mult_transpose(vec1_n2, flux);
                Fvec.add(-JxW[qp], vec1_n2);
                
                // contribution from diffusion flux
                flux.zero();
                for (unsigned int i_dim=0; i_dim<dim; i_dim++)
                {
                    this->calculate_diffusion_flux(i_dim, p_sol, stress_tensor,
                                                   temp_grad, vec2_n1);
                    flux.add(face_normals[qp](i_dim), vec2_n1); // fi ni
                }
                B_mat.vector_mult_transpose(vec1_n2, flux);
                Fvec.add(JxW[qp], vec1_n2);
                
                
                if ( request_jacobian && c.get_elem_solution_derivative() )
                {
                    this->calculate_advection_flux_jacobian_for_moving_solid_wall_boundary
                    (p_sol, ui_ni, face_normals[qp], dnormal, A_mat);
                    
                    // contribution from advection flux
                    B_mat.left_multiply(mat1_n1n2, A_mat);
                    B_mat.right_multiply_transpose(mat2_n2n2, mat1_n1n2);
                    Kmat.add(-JxW[qp], mat2_n2n2);
                    
                    // contribution from diffusion flux
                    for (unsigned int i_dim=0; i_dim<dim; i_dim++)
                        for (unsigned int deriv_dim=0; deriv_dim<dim; deriv_dim++)
                        {
                            this->calculate_diffusion_flux_jacobian
                            (i_dim, deriv_dim, p_sol, mat_n1n1); // Kij
                            dB_mat[deriv_dim].left_multiply
                            (mat1_n1n2, mat_n1n1); // Kij dB/dx_j
                            B_mat.right_multiply_transpose
                            (mat2_n2n2, mat1_n1n2); // B^T Kij dB/dx_j
                            Kmat.add(JxW[qp], mat2_n2n2);
                        }
                }
            }
        }
            break;
            
        case SYMMETRY_WALL: // inviscid boundary condition without any diffusive component
                        // conditions enforced are
                        // vi ni = 0       (slip wall)
                        // tau_ij nj = 0   (because velocity gradient at wall = 0)
                        // qi ni = 0       (since heat flux occurs only on no-slip wall and far-field bc)
        {
            for (unsigned int qp=0; qp<qpoint.size(); qp++)
            {
                // first update the variables at the current quadrature point
                this->update_solution_at_quadrature_point
                (vars, qp, c,
                 false, c.get_elem_solution(),
                 conservative_sol, p_sol,
                 B_mat, dB_mat);
                
                flux.zero();
                // since vi=0, vi ni = 0, so the advection flux gets evaluated
                // using the slip wall condition
                for (unsigned int i_dim=0; i_dim<dim; i_dim++)
                    flux(i_dim+1) += p_sol.p * face_normals[qp](i_dim);
                
                B_mat.vector_mult_transpose(vec1_n2, flux);
                Fvec.add(-JxW[qp], vec1_n2);
                
                if ( request_jacobian && c.get_elem_solution_derivative() )
                {
                    this->calculate_advection_flux_jacobian_for_moving_solid_wall_boundary
                    (p_sol, 0., face_normals[qp], dnormal, A_mat);
                    
                    B_mat.left_multiply(mat1_n1n2, A_mat);
                    B_mat.right_multiply_transpose(mat2_n2n2, mat1_n1n2);
                    Kmat.add(-JxW[qp], mat2_n2n2);
                }
            }
        }
            break;

        
        case SLIP_WALL: // inviscid boundary condition without any diffusive component
                        // conditions enforced are
                        // vi ni = wi_dot (ni + dni) - ui dni   (moving slip wall with deformation)
                        // tau_ij nj = 0   (because velocity gradient at wall = 0)
                        // qi ni = 0       (since heat flux occurs only on no-slip wall and far-field bc)
        {
            Real ui_ni = 0.;
            
            for (unsigned int qp=0; qp<qpoint.size(); qp++)
            {
                // first update the variables at the current quadrature point
                this->update_solution_at_quadrature_point
                (vars, qp, c,
                 false, c.get_elem_solution(),
                 conservative_sol, p_sol,
                 B_mat, dB_mat);
                
                // copy the surface normal
                for (unsigned int i_dim=0; i_dim<dim; i_dim++)
                    local_normal(i_dim) = face_normals[qp](i_dim);

                // now check if the surface deformation is defined and
                // needs to be applied through transpiration boundary
                // condition
                ui_ni = 0.;
                p_sol.get_uvec(uvec);
                
                if (surface_motion) // get the surface motion data
                {
                    surface_motion->surface_velocity
                    (this->time, qpoint[qp], face_normals[qp],
                     surface_def, surface_vel, dnormal);
                    
                    // update the normal with the deformation
                    // this defines the normal of the surface that has been
                    // deformed, although the geometry of the flow mesh does
                    // not conform to that deformation
                    local_normal.add(1., dnormal);
                    
                    //    ui * (ni + dni) = wi_dot * (ni + dni)
                    // => ui * ni = wi_dot * (ni + dni) - ui * dni
                    
                    // now initialize the surface velocity
                    // note that the perturbed local normal is used here
                    // since it resembles the normal of the surface that
                    // has undergone a static deformation, even though the
                    // surface_vel might be identically zero for a static
                    // body.
                    ui_ni  = surface_vel.dot(local_normal);
                    ui_ni -= uvec.dot(dnormal);
                }
                
                flux.zero();
                flux.add(ui_ni, conservative_sol);
                flux(n1-1) += p_sol.p*ui_ni;
                for (unsigned int i_dim=0; i_dim<dim; i_dim++)
                    flux(i_dim+1) += p_sol.p * face_normals[qp](i_dim);
                
                B_mat.vector_mult_transpose(vec1_n2, flux);
                Fvec.add(-JxW[qp], vec1_n2);
                
                if ( request_jacobian && c.get_elem_solution_derivative() )
                {
                    this->calculate_advection_flux_jacobian_for_moving_solid_wall_boundary
                    (p_sol, ui_ni, face_normals[qp], dnormal, A_mat);
                    
                    B_mat.left_multiply(mat1_n1n2, A_mat);
                    B_mat.right_multiply_transpose(mat2_n2n2, mat1_n1n2);
                    Kmat.add(-JxW[qp], mat2_n2n2);
                }
            }
        }
            break;
            
            
        case FAR_FIELD:
            // conditions enforced are:
            // -- f_adv_i ni =  f_adv = f_adv(+) + f_adv(-)     (flux vector splitting for advection)
            // -- f_diff_i ni  = f_diff                         (evaluation of diffusion flux based on domain solution)
        {
            for (unsigned int qp=0; qp<qpoint.size(); qp++)
            {
                // first update the variables at the current quadrature point
                this->update_solution_at_quadrature_point
                (vars, qp, c,
                 false, c.get_elem_solution(),
                 conservative_sol, p_sol,
                 B_mat, dB_mat);
                
                this->calculate_advection_left_eigenvector_and_inverse_for_normal
                (p_sol, face_normals[qp], eig_val,
                 l_eig_vec, l_eig_vec_inv_tr);
                
                // for all eigenalues that are less than 0, the characteristics are coming into the domain, hence,
                // evaluate them using the given solution.
                mat_n1n1 = l_eig_vec_inv_tr;
                for (unsigned int j=0; j<n1; j++)
                    if (eig_val(j, j) < 0)
                        mat_n1n1.scale_column(j, eig_val(j, j)); // L^-T [omaga]_{-}
                    else
                        mat_n1n1.scale_column(j, 0.0);
                
                mat_n1n1.right_multiply_transpose(l_eig_vec); // A_{-} = L^-T [omaga]_{-} L^T
                
                this->get_infinity_vars( U_vec_interpolated );
                mat_n1n1.vector_mult(flux, U_vec_interpolated);  // f_{-} = A_{-} B U
                
                B_mat.vector_mult_transpose(vec1_n2, flux); // B^T f_{-}   (this is flux coming into the solution domain)
                Fvec.add(-JxW[qp], vec1_n2);
                
                // now calculate the flux for eigenvalues greater than 0,
                // the characteristics go out of the domain, so that
                // the flux is evaluated using the local solution
                mat_n1n1 = l_eig_vec_inv_tr;
                for (unsigned int j=0; j<n1; j++)
                    if (eig_val(j, j) > 0)
                        mat_n1n1.scale_column(j, eig_val(j, j)); // L^-T [omaga]_{+}
                    else
                        mat_n1n1.scale_column(j, 0.0);
                
                mat_n1n1.right_multiply_transpose(l_eig_vec); // A_{+} = L^-T [omaga]_{+} L^T
                mat_n1n1.vector_mult(flux, conservative_sol); // f_{+} = A_{+} B U
                
                B_mat.vector_mult_transpose(vec1_n2, flux); // B^T f_{+}   (this is flux going out of the solution domain)
                Fvec.add(-JxW[qp], vec1_n2);
                
                
                if (_if_viscous) // evaluate the viscous flux using the domain solution
                {
                    // stress tensor and temperature gradient
                    this->calculate_conservative_variable_jacobian
                    (p_sol, dcons_dprim, dprim_dcons);
                    this->calculate_diffusion_tensors
                    (c.get_elem_solution(), dB_mat, dprim_dcons,
                     p_sol, stress_tensor, temp_grad);
                    
                    // contribution from diffusion flux
                    flux.zero();
                    for (unsigned int i_dim=0; i_dim<dim; i_dim++)
                    {
                        this->calculate_diffusion_flux
                        (i_dim, p_sol, stress_tensor,
                         temp_grad, vec2_n1);
                        flux.add(face_normals[qp](i_dim), vec2_n1); // fi ni
                    }
                    B_mat.vector_mult_transpose(vec1_n2, flux);
                    Fvec.add(JxW[qp], vec1_n2);
                }
                
                
                if ( request_jacobian && c.get_elem_solution_derivative() )
                {
                    // terms with negative eigenvalues do not contribute to the Jacobian
                    
                    // now calculate the Jacobian for eigenvalues greater than 0,
                    // the characteristics go out of the domain, so that
                    // the flux is evaluated using the local solution
                    mat_n1n1 = l_eig_vec_inv_tr;
                    for (unsigned int j=0; j<n1; j++)
                        if (eig_val(j, j) > 0)
                            mat_n1n1.scale_column(j, eig_val(j, j)); // L^-T [omaga]_{+}
                        else
                            mat_n1n1.scale_column(j, 0.0);
                    mat_n1n1.right_multiply_transpose(l_eig_vec); // A_{+} = L^-T [omaga]_{+} L^T
                    B_mat.left_multiply(mat1_n1n2, mat_n1n1);
                    B_mat.right_multiply_transpose(mat2_n2n2, mat1_n1n2); // B^T A_{+} B   (this is flux going out of the solution domain)
                    
                    Kmat.add(-JxW[qp], mat2_n2n2);
                    
                    if (_if_viscous)
                    {
                        // contribution from diffusion flux
                        for (unsigned int i_dim=0; i_dim<dim; i_dim++)
                            for (unsigned int deriv_dim=0; deriv_dim<dim; deriv_dim++)
                            {
                                this->calculate_diffusion_flux_jacobian
                                (i_dim, deriv_dim, p_sol, mat_n1n1); // Kij
                                dB_mat[deriv_dim].left_multiply
                                (mat1_n1n2, mat_n1n1); // Kij dB/dx_j
                                B_mat.right_multiply_transpose
                                (mat2_n2n2, mat1_n1n2); // B^T Kij dB/dx_j
                                Kmat.add(JxW[qp], mat2_n2n2);
                            }
                    }
                }
            }
        }
            break;

            
        case EXHAUST:
            // conditions enforced are:
            // -- f_adv_i ni =  f_adv = f_adv(+) + f_adv(-)     (flux vector splitting for advection)
            // -- f_adv(+) is evaluated using interior solution
            // -- f_adv(-) is evaluated using a combined solution where one variable is provided and the others are used from the interior
            // -- f_diff_i ni  = f_diff                         (evaluation of diffusion flux based on domain solution)
        {
            for (unsigned int qp=0; qp<qpoint.size(); qp++)
            {
                // first update the variables at the current quadrature point
                this->update_solution_at_quadrature_point
                (vars, qp, c,
                 false, c.get_elem_solution(),
                 conservative_sol, p_sol,
                 B_mat, dB_mat);
                
                this->calculate_advection_left_eigenvector_and_inverse_for_normal
                (p_sol, face_normals[qp], eig_val,
                 l_eig_vec, l_eig_vec_inv_tr);
                
                // for all eigenalues that are less than 0, the characteristics are coming into the domain, hence,
                // evaluate them using the given solution.
                mat_n1n1 = l_eig_vec_inv_tr;
                for (unsigned int j=0; j<n1; j++)
                    if (eig_val(j, j) < 0)
                        mat_n1n1.scale_column(j, eig_val(j, j)); // L^-T [omaga]_{-}
                    else
                        mat_n1n1.scale_column(j, 0.0);
                
                mat_n1n1.right_multiply_transpose(l_eig_vec); // A_{-} = L^-T [omaga]_{-} L^T
                
                // instead of using the infinity vars, the constrained variable is used
                // and the others are used from the local solution
                // use constrained value at exhaust to create the U{-} solution vector
                Real rho_constrained = 0.7519;
                U_vec_interpolated = conservative_sol;
                U_vec_interpolated.scale(rho_constrained/conservative_sol(0));  // value is not constrained to density at exhaust

                // use the constrained vector to calculate the flux value
                mat_n1n1.vector_mult(flux, U_vec_interpolated);  // f_{-} = A_{-} B U{-}
                
                B_mat.vector_mult_transpose(vec1_n2, flux); // B^T f_{-}   (this is flux coming into the solution domain)
                Fvec.add(-JxW[qp], vec1_n2);
                
                // now calculate the flux for eigenvalues greater than 0,
                // the characteristics go out of the domain, so that
                // the flux is evaluated using the local solution
                mat_n1n1 = l_eig_vec_inv_tr;
                for (unsigned int j=0; j<n1; j++)
                    if (eig_val(j, j) > 0)
                        mat_n1n1.scale_column(j, eig_val(j, j)); // L^-T [omaga]_{+}
                    else
                        mat_n1n1.scale_column(j, 0.0);
                
                mat_n1n1.right_multiply_transpose(l_eig_vec); // A_{+} = L^-T [omaga]_{+} L^T
                mat_n1n1.vector_mult(flux, U_vec_interpolated); // f_{+} = A_{+} B U{+}
                
                B_mat.vector_mult_transpose(vec1_n2, flux); // B^T f_{+}   (this is flux going out of the solution domain)
                Fvec.add(-JxW[qp], vec1_n2);
                
                
                if (_if_viscous) // evaluate the viscous flux using the domain solution
                {
                    // stress tensor and temperature gradient
                    this->calculate_conservative_variable_jacobian
                    (p_sol, dcons_dprim, dprim_dcons);
                    this->calculate_diffusion_tensors
                    (c.get_elem_solution(), dB_mat, dprim_dcons,
                     p_sol, stress_tensor, temp_grad);
                    
                    // contribution from diffusion flux
                    flux.zero();
                    for (unsigned int i_dim=0; i_dim<dim; i_dim++)
                    {
                        this->calculate_diffusion_flux
                        (i_dim, p_sol, stress_tensor,
                         temp_grad, vec2_n1);
                        flux.add(face_normals[qp](i_dim), vec2_n1); // fi ni
                    }
                    B_mat.vector_mult_transpose(vec1_n2, flux);
                    Fvec.add(JxW[qp], vec1_n2);
                }
                
                
                if ( request_jacobian && c.get_elem_solution_derivative() )
                {
                    // terms with negative eigenvalues will contribute to the Jacobian, but only for
                    // the unconstrained variables
                    mat_n1n1 = l_eig_vec_inv_tr;
                    for (unsigned int j=0; j<n1; j++)
                        if (eig_val(j, j) < 0)
                            mat_n1n1.scale_column(j, eig_val(j, j)); // L^-T [omaga]_{-}
                        else
                            mat_n1n1.scale_column(j, 0.0);
                    mat_n1n1.right_multiply_transpose(l_eig_vec); // A_{-} = L^-T [omaga]_{-} L^T
                    
                    // now project this matrix to the constraint at the exhause
                    // df{-}/dU_c = A{-} dU{-}/dU_c
                    // U{-} = U{-} (U_p, U_p^constrained)
                    // dU{-}/dU_c = dU{-}/dU_p dU_p/dU_c + dU{-}/dU_p^constrained dU_p^constrained/dU_c
                    // the last term is zero because U_p^constrained is not a function of the interior conservative variables
                    // the first term requires the jacobian of conservative and primitive variables
                    this->calculate_conservative_variable_jacobian
                    (p_sol, dcons_dprim, dprim_dcons);
                    
                    // the term dU{-}/dU_p is calculated by scaling the value from interior solution
                    mat3_n1n1 = dcons_dprim;
                    mat3_n1n1.scale(rho_constrained/conservative_sol(0));
                    // zero the first column since the U{-} is not a function of interior density value
                    for (unsigned int i=0; i<n1; i++)
                        mat3_n1n1(i,0) = 0.;
                    
                    // now multiply the matrices and get the
                    mat3_n1n1.right_multiply(dprim_dcons);
                    mat_n1n1.right_multiply(mat3_n1n1); // this is the jacobian matrix projected to the constraint of the exhaust values
                    
                    B_mat.left_multiply(mat1_n1n2, mat_n1n1);
                    B_mat.right_multiply_transpose(mat2_n2n2, mat1_n1n2); // B^T A_{-} B   (this is flux coming into the solution domain)
                    
                    Kmat.add(-JxW[qp], mat2_n2n2);
                    
                    
                    // now calculate the Jacobian for eigenvalues greater than 0,
                    // the characteristics go out of the domain, so that
                    // the flux is evaluated using the local solution
                    mat_n1n1 = l_eig_vec_inv_tr;
                    for (unsigned int j=0; j<n1; j++)
                        if (eig_val(j, j) > 0)
                            mat_n1n1.scale_column(j, eig_val(j, j)); // L^-T [omaga]_{+}
                        else
                            mat_n1n1.scale_column(j, 0.0);
                    mat_n1n1.right_multiply_transpose(l_eig_vec); // A_{+} = L^-T [omaga]_{+} L^T

                    // now do the same operation on the outoging characteristics
                    // the term dU{-}/dU_p is calculated by scaling the value from interior solution
                    mat3_n1n1 = dcons_dprim;
                    mat3_n1n1.scale(rho_constrained/conservative_sol(0));
                    // zero the first column since the U{-} is not a function of interior density value
                    for (unsigned int i=0; i<n1; i++)
                        mat3_n1n1(i,0) = 0.;
                    
                    // now multiply the matrices and get the
                    mat3_n1n1.right_multiply(dprim_dcons);
                    mat_n1n1.right_multiply(mat3_n1n1); // this is the jacobian matrix projected to the constraint of the exhaust values

                    B_mat.left_multiply(mat1_n1n2, mat_n1n1);
                    B_mat.right_multiply_transpose(mat2_n2n2, mat1_n1n2); // B^T A_{+} B   (this is flux going out of the solution domain)

                    Kmat.add(-JxW[qp], mat2_n2n2);
                    
                    if (_if_viscous)
                    {
                        // contribution from diffusion flux
                        for (unsigned int i_dim=0; i_dim<dim; i_dim++)
                            for (unsigned int deriv_dim=0; deriv_dim<dim; deriv_dim++)
                            {
                                this->calculate_diffusion_flux_jacobian
                                (i_dim, deriv_dim, p_sol, mat_n1n1); // Kij
                                dB_mat[deriv_dim].left_multiply
                                (mat1_n1n2, mat_n1n1); // Kij dB/dx_j
                                B_mat.right_multiply_transpose
                                (mat2_n2n2, mat1_n1n2); // B^T Kij dB/dx_j
                                Kmat.add(JxW[qp], mat2_n2n2);
                            }
                    }
                }
            }
        }
            break;

    }
    
    
    //    libMesh::out << "inside side constraint " << std::endl;
    //    libMesh::out << "elem solution" << std::endl; c.elem_solution.print(libMesh::out);
    //    libMesh::out << if_inf_bc << "  " << if_wall_bc << std::endl;
    //    libMesh::out << "bc vec: " << std::endl; Fvec.print(libMesh::out);
    //    if (request_jacobian && c.elem_solution_derivative)
    //        Kmat.print(libMesh::out);
    
    return request_jacobian;
}





bool FluidSystem::mass_residual (bool request_jacobian,
                                 libMesh::DiffContext &context)
{
    libMesh::FEMContext &c = libMesh::libmesh_cast_ref<libMesh::FEMContext&>(context);
    
    libMesh::FEBase* elem_fe;
    
    c.get_element_fe(vars[0], elem_fe);
    
    
    // Element Jacobian * quadrature weights for interior integration
    const std::vector<Real> &JxW = elem_fe->get_JxW();
    
    // The subvectors and submatrices we need to fill:
    DenseRealMatrix& Kmat = c.get_elem_jacobian();
    DenseRealVector& Fvec = c.get_elem_residual();
    
    // Now we will build the element Jacobian and residual.
    // Constructing the residual requires the solution and its
    // gradient from the previous timestep.  This must be
    // calculated at each quadrature point by summing the
    // solution degree-of-freedom values by the appropriate
    // weight functions.
    const unsigned int n_qpoints = c.get_element_qrule().n_points(), n1 = dim+2,
    n_dofs = n1 * elem_fe->n_shape_functions();
    
    FEMOperatorMatrix B_mat;
    std::vector<FEMOperatorMatrix> dB_mat(dim);
    std::vector<DenseRealMatrix > Ai_advection(dim);
    DenseRealMatrix LS_mat, LS_sens, Ai_Bi_advection, mat_n1n2,
    mat2_n2n2, mat3;
    DenseRealVector flux, vec1_n1, vec3_n2,
    conservative_sol, dc_ref_sol;
    LS_mat.resize(n1, n_dofs); LS_sens.resize(n_dofs, n_dofs);
    mat2_n2n2.resize(n_dofs, n_dofs);
    Ai_Bi_advection.resize(dim+2, n_dofs);  mat_n1n2.resize(dim+2, n_dofs);
    flux.resize(n1); vec1_n1.resize(n1); vec3_n2.resize(n_dofs);
    conservative_sol.resize(dim+2);
    for (unsigned int i=0; i<dim; i++)
        Ai_advection[i].resize(dim+2, dim+2);
    
    std::vector<std::vector<DenseRealMatrix > > flux_jacobian_sens;
    flux_jacobian_sens.resize(dim);
    for (unsigned int i_dim=0; i_dim<dim; i_dim++)
    {
        flux_jacobian_sens[i_dim].resize(n1); // number of variables for sensitivity
        for (unsigned int i_cvar=0; i_cvar<n1; i_cvar++)
            flux_jacobian_sens[i_dim][i_cvar].resize(n1, n1);
    }

    if (!if_use_stored_dc_coeff)
        dc_ref_sol = c.get_elem_fixed_solution();
    else
        dc_ref_sol = c.get_localized_vector(*_dc_ref_sol);

    PrimitiveSolution primitive_sol;

    for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
        // first update the variables at the current quadrature point
        this->update_solution_at_quadrature_point
        (vars, qp, c,
         true, c.get_elem_fixed_solution(),
         conservative_sol, primitive_sol,
         B_mat, dB_mat);
        
        for ( unsigned int i_dim=0; i_dim < dim; i_dim++)
            this->calculate_advection_flux_jacobian ( i_dim, primitive_sol,
                                                     Ai_advection[i_dim] );
        
        Ai_Bi_advection.zero();
        
        for (unsigned int i_dim=0; i_dim<dim; i_dim++)
        {
            dB_mat[i_dim].left_multiply(mat_n1n2, Ai_advection[i_dim]);
            Ai_Bi_advection.add(1.0, mat_n1n2);
        }
        
        if (_if_update_stabilization_per_quadrature_point || (qp == 0)) {
            this->calculate_differential_operator_matrix
            (vars, qp, c, c.get_elem_fixed_solution(),
             primitive_sol, B_mat, dB_mat,
             Ai_advection, Ai_Bi_advection,
             flux_jacobian_sens, LS_mat, LS_sens);
        }
        // Galerkin contribution to velocity
        B_mat.vector_mult( vec1_n1, c.get_elem_solution() );
        B_mat.vector_mult_transpose(vec3_n2, vec1_n1);
        Fvec.add(JxW[qp], vec3_n2);
        
        // LS contribution to velocity
        B_mat.vector_mult( vec1_n1, c.get_elem_solution() );
        LS_mat.vector_mult_transpose(vec3_n2, vec1_n1);
        Fvec.add(JxW[qp], vec3_n2);
        
        
        if (request_jacobian && c.get_elem_solution_derivative())
        {
            // contribution from unsteady term
            // Galerkin contribution of velocity
            B_mat.right_multiply_transpose(mat2_n2n2, B_mat);
            Kmat.add(JxW[qp], mat2_n2n2);
            
            // LS contribution of velocity
            LS_mat.get_transpose(mat3);
            B_mat.left_multiply(mat2_n2n2, mat3); // LS^T tau Bmat
            Kmat.add(JxW[qp], mat2_n2n2);
        }
    } // end of the quadrature point qp-loop
    
    //    libMesh::out << "inside mass residual " << std::endl;
    //    libMesh::out << "elem velocity" << std::endl; c.elem_solution.print(libMesh::out);
    //    libMesh::out << "elem solution" << std::endl; c.elem_fixed_solution.print(libMesh::out);
    //    libMesh::out << "mass vec: " << std::endl; Fvec.print(libMesh::out);
    //    if (request_jacobian && c.elem_solution_derivative)
    //        Kmat.print(libMesh::out);
    
    return request_jacobian;
}



void FluidSystem::postprocess()
{
    for (unsigned int i=0; i<this->qoi.size(); i++)
        if (i == 0)
            libMesh::out << "Lift: " << std::setw(25)
            << std::setprecision(14) << qoi[i] << std::endl;
        else if (i == 1)
            libMesh::out << "Drag: " << std::setw(25)
            << std::setprecision(14) << qoi[i] << std::endl;
        else if (i == 2)
            libMesh::out << "Total V: " << std::setw(25)
            << std::setprecision(14) << qoi[i] << std::endl;
        else if (i == 3)
            libMesh::out << "Entropy Error: " << std::setw(25)
            << std::setprecision(14) << sqrt(qoi[i]/qoi[i-1]) << std::endl;
}





void FluidSystem::evaluate_recalculate_dc_flag()
{
    Real norm = this->calculate_norm(*(this->solution), this->vars[0], libMesh::L2),
    relative_change0 = fabs(_rho_norm_curr - _rho_norm_old)/_rho_norm_curr,
    relative_change1 = fabs(norm - _rho_norm_curr)/norm;
    
    libMesh::out
    << "Rho L2-norm delta: " << relative_change0
    << " , " << relative_change1 << std::endl;
    
    // if the relative change in the density norm is less than the threshold, then hold
    // the dc_coeffs constant
    if ((relative_change0 < this->dc_recalculate_tolerance) &&
        (relative_change1 < this->dc_recalculate_tolerance))
    {
        // if the previous flag was false, then update the solution
        if (!this->if_use_stored_dc_coeff) {
            _dc_ref_sol.reset(libMesh::NumericVector<Real>::build(this->comm()).release());
            _dc_ref_sol->init(this->n_dofs(), this->n_local_dofs(),
                              this->get_dof_map().get_send_list(),
                              false, libMesh::GHOSTED);
            this->solution->localize(*_dc_ref_sol, this->get_dof_map().get_send_list());
        }
        
        // use the last updated dc_ref_sol
        this->if_use_stored_dc_coeff = true;
        libMesh::out << "Using stored dc coeff" << std::endl;
    }
    else
    {
        // update the dc_ref_sol
        _dc_ref_sol.reset(libMesh::NumericVector<Real>::build(this->comm()).release());
        _dc_ref_sol->init(this->n_dofs(), this->n_local_dofs(),
                          this->get_dof_map().get_send_list(),
                          false, libMesh::GHOSTED);
        this->solution->localize(*_dc_ref_sol, this->get_dof_map().get_send_list());
        
        if_use_stored_dc_coeff = false;
        libMesh::out << "Recalculating dc coeff" << std::endl;
    }
    
    _rho_norm_old = _rho_norm_curr;
    _rho_norm_curr = norm;
}





Real get_var_val(const std::string& var_name, const PrimitiveSolution& p_sol,
                 Real p0, Real q0)
{
    if (var_name == "ux")
        return p_sol.u1;
    else if (var_name == "uy")
        return p_sol.u2;
    else if (var_name == "uz")
        return p_sol.u3;
    else if (var_name == "T")
        return p_sol.T;
    else if (var_name == "s")
        return p_sol.entropy;
    else if (var_name == "p")
        return p_sol.p;
    else if (var_name == "cp")
        return p_sol.c_pressure(p0, q0);
    else if (var_name == "a")
        return p_sol.a;
    else if (var_name == "M")
        return p_sol.mach;
    else
        libmesh_assert(false);
}



class PrimitiveFEMFunction : public libMesh::FEMFunctionBase<Real>
{
public:
    // Constructor
    PrimitiveFEMFunction(libMesh::AutoPtr<libMesh::FunctionBase<Real> > fluid_func,
                         std::vector<std::string>& vars,
                         Real cp, Real cv, Real p0, Real q0):
    libMesh::FEMFunctionBase<Real>(),
    _fluid_function(fluid_func.release()),
    _vars(vars), _cp(cp), _cv(cv), _p0(p0), _q0(q0)
    {
        _fluid_function->init();
    }
    
    // Destructor
    virtual ~PrimitiveFEMFunction () {}
    
    virtual libMesh::AutoPtr<libMesh::FEMFunctionBase<Real> > clone () const
    {return libMesh::AutoPtr<libMesh::FEMFunctionBase<Real> >( new PrimitiveFEMFunction
                                              (_fluid_function->clone(), _vars, _cp, _cv,
                                               _p0, _q0) ); }
    
    virtual void operator() (const libMesh::FEMContext& c, const libMesh::Point& p,
                             const Real t, DenseRealVector& val)
    {
        DenseRealVector fluid_sol;
        (*_fluid_function)(p, t, fluid_sol);
        PrimitiveSolution p_sol;
        
        p_sol.init(c.get_dim(), fluid_sol, _cp, _cv, false);
        
        for (unsigned int i=0; i<_vars.size(); i++)
            val(i) = get_var_val(_vars[i], p_sol, _p0, _q0);
    }
    
    
    virtual Real component(const libMesh::FEMContext& c, unsigned int i_comp,
                             const libMesh::Point& p, Real t=0.)
    {
        DenseRealVector fluid_sol;
        (*_fluid_function)(p, t, fluid_sol);
        PrimitiveSolution p_sol;
        
        p_sol.init(c.get_dim(), fluid_sol, _cp, _cv, false);
        
        return get_var_val(_vars[i_comp], p_sol, _p0, _q0);
    }
    
    
    virtual Real operator() (const libMesh::FEMContext&, const libMesh::Point& p,
                               const Real time = 0.)
    {libmesh_error();}
    
private:
    
    libMesh::AutoPtr<libMesh::FunctionBase<Real> > _fluid_function;
    std::vector<std::string>& _vars;
    Real _cp, _cv, _p0, _q0;
};




void FluidPostProcessSystem::init_data()
{
    const unsigned int dim = this->get_mesh().mesh_dimension();
    
    const libMesh::FEFamily fefamily = libMesh::LAGRANGE;
    const libMesh::Order order = libMesh::FIRST;
    
    u = this->add_variable("ux", order, fefamily);
    if (dim > 1)
        v = this->add_variable("uy", order, fefamily);
    if (dim > 2)
        w = this->add_variable("uz", order, fefamily);
    T = this->add_variable("T", order, fefamily);
    s = this->add_variable("s", order, fefamily);
    p = this->add_variable("p", order, fefamily);
    cp = this->add_variable("cp", order, fefamily);
    a = this->add_variable("a", order, fefamily);
    M = this->add_variable("M", order, fefamily);
    
    libMesh::System::init_data();
}





void FluidPostProcessSystem::postprocess()
{
    // get the solution vector from
    const FluidSystem& fluid =
    this->get_equation_systems().get_system<FluidSystem>("FluidSystem");
    fluid.get_mesh().sub_point_locator();
    
    std::vector<unsigned int> fluid_vars(fluid.n_vars());
    std::vector<std::string> post_process_var_names(this->n_vars());
    
    for (unsigned int i=0; i<fluid_vars.size(); i++)
        fluid_vars[i] = i;
    for (unsigned int i=0; i<this->n_vars(); i++)
        post_process_var_names[i] = this->variable_name(i);
    
    
    libMesh::AutoPtr<libMesh::NumericVector<Real> > soln =
    libMesh::NumericVector<Real>::build(this->get_equation_systems().comm());
    soln->init(fluid.solution->size(), true, libMesh::SERIAL);
    fluid.solution->localize(*soln,
                             fluid.get_dof_map().get_send_list());
    
    
    libMesh::AutoPtr<libMesh::MeshFunction> mesh_function
    (new libMesh::MeshFunction(this->get_equation_systems(), *soln,
                      fluid.get_dof_map(), fluid_vars));
    mesh_function->init();
    
    libMesh::AutoPtr<libMesh::FEMFunctionBase<Real> > post_process_function
    (new PrimitiveFEMFunction(mesh_function->clone(), post_process_var_names,
                              flight_condition->gas_property.cp,
                              flight_condition->gas_property.cv,
                              flight_condition->p0(),
                              flight_condition->q0()));
    
    this->project_solution(post_process_function.get());
    
    this->update();
}






