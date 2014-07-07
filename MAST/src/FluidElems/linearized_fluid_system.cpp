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


// MAST includes
#include "FluidElems/linearized_fluid_system.h"


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






void LinearizedFluidSystem::init_data ()
{
    this->use_fixed_solution = true;
    
    dim = this->get_mesh().mesh_dimension();
    
    vars.resize(dim+2);
    
    // initialize the fluid values
    FluidElemBase::init_data();
    
    // Check the input file for Reynolds number, application type,
    // approximation type
    
    // set parameter values
    unsigned int o = _infile("fe_order", 1);
    std::string fe_family = _infile("fe_family", std::string("LAGRANGE"));
    libMesh::FEFamily fefamily = libMesh::Utility::string_to_enum<libMesh::FEFamily>(fe_family);
    
    vars[0]  = this->add_variable ( "drho", static_cast<libMesh::Order>(o), fefamily);
    this->time_evolving(vars[0]);
    
    vars[1] = this->add_variable ("drhoux", static_cast<libMesh::Order>(o), fefamily);
    this->time_evolving(vars[1]);
    
    if (dim > 1)
    {
        vars[2] = this->add_variable ("drhouy", static_cast<libMesh::Order>(o), fefamily);
        this->time_evolving(vars[2]);
    }
    
    if (dim > 2)
    {
        vars[3] = this->add_variable ("drhouz", static_cast<libMesh::Order>(o), fefamily);
        this->time_evolving(vars[3]);
    }
    
    vars[dim+2-1] = this->add_variable ("drhoe", static_cast<libMesh::Order>(o), fefamily);
    this->time_evolving(vars[dim+2-1]);
    
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



void LinearizedFluidSystem::init_context(libMesh::DiffContext &context)
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



void LinearizedFluidSystem::localize_fluid_solution()
{
    libmesh_assert(!_if_localized_sol);
    
    _local_fluid_solution =
    libMesh::NumericVector<Real>::build(this->get_equation_systems().comm());
    
    
    libMesh::System& fluid =
    this->get_equation_systems().get_system<libMesh::System>("FluidSystem");
    
    _local_fluid_solution->init(fluid.solution->size(), true, libMesh::SERIAL);
    fluid.solution->localize(*_local_fluid_solution,
                             fluid.get_dof_map().get_send_list());
    
    _if_localized_sol = true;
}


bool LinearizedFluidSystem::element_time_derivative (bool request_jacobian,
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
    conservative_sol, temp_grad, diff_val, diff_val2, dc_ref_sol,
    ref_sol, total_sol, conservative_deltasol;
    
    LS_mat.resize(n1, n_dofs); LS_sens.resize(n_dofs, n_dofs);
    Ai_Bi_advection.resize(dim+2, n_dofs);
    mat_n1n2.resize(dim+2, n_dofs); A_sens.resize(n1, n_dofs);
    stress_tensor.resize(dim, dim); dprim_dcons.resize(dim+2, dim+2);
    dcons_dprim.resize(dim+2, dim+2);
    mat2_n2n2.resize(n_dofs, n_dofs);
    
    flux.resize(n1); vec1_n1.resize(n1); vec2_n1.resize(n1);
    vec3_n2.resize(n_dofs); conservative_sol.resize(dim+2);
    conservative_deltasol.resize(dim+2);
    temp_grad.resize(dim); diff_val.resize(dim); diff_val2.resize(dim);
    ref_sol.resize(n_dofs); total_sol.resize(n_dofs);
    
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
    
    // reference solution for small-disturbance
    libMesh::System& fluid =
    this->get_equation_systems().get_system<libMesh::System>("FluidSystem");

    std::vector<libMesh::dof_id_type> fluid_dof_indices;
    fluid.get_dof_map().dof_indices(&c.get_elem(), fluid_dof_indices);
    
    for (unsigned int i=0; i<n_dofs; i++)
        ref_sol(i) = (*_local_fluid_solution)(fluid_dof_indices[i]);
    total_sol = ref_sol;
    total_sol += c.get_elem_solution();

    PrimitiveSolution primitive_sol;
    SmallPerturbationPrimitiveSolution<Real> delta_p_sol;
    
    // assuming that all variables have the same interpolation
    const std::vector<std::vector<Real> >& phi = elem_fe->get_phi();
    const unsigned int n_phi = phi.size();
    
    for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
        // first update the variables at the current quadrature point
        this->update_solution_at_quadrature_point(vars, qp, c,
                                                  true,
                                                  ref_sol,
                                                  conservative_sol,
                                                  primitive_sol,
                                                  B_mat,
                                                  dB_mat);
        
        for ( unsigned int i_dim=0; i_dim < dim; i_dim++)
            this->calculate_advection_flux_jacobian ( i_dim, primitive_sol,
                                                     Ai_advection[i_dim] );
        
        Ai_Bi_advection.zero();
        
        for (unsigned int i_dim=0; i_dim<dim; i_dim++)
        {
            dB_mat[i_dim].left_multiply(mat_n1n2, Ai_advection[i_dim]);
            Ai_Bi_advection.add(1.0, mat_n1n2);
        }
        
        // initialize the delta_p_sol
        B_mat.vector_mult(conservative_deltasol, c.get_elem_solution());
        
        delta_p_sol.zero();
        delta_p_sol.init(primitive_sol, conservative_deltasol);

        
        if (_if_full_linearization)
        {
            for (unsigned int i_dim=0; i_dim<dim; i_dim++)
                this->calculate_advection_flux_jacobian_sensitivity_for_conservative_variable
                (i_dim, primitive_sol, flux_jacobian_sens[i_dim]);
        }
        if (_if_update_stabilization_per_quadrature_point || (qp == 0)) {
            this->calculate_differential_operator_matrix
            (vars, qp, c,
             ref_sol, primitive_sol,
             B_mat, dB_mat, Ai_advection,
             Ai_Bi_advection, flux_jacobian_sens,
             LS_mat, LS_sens);
            
            this->calculate_aliabadi_discontinuity_operator(vars, qp, c,
                                                            primitive_sol,
                                                            ref_sol,
                                                            dB_mat,
                                                            Ai_Bi_advection,
                                                            diff_val);
            
            calculate_small_disturbance_aliabadi_discontinuity_operator(vars, qp, c,
                                                                        primitive_sol,
                                                                        delta_p_sol,
                                                                        total_sol,
                                                                        dB_mat,
                                                                        Ai_Bi_advection,
                                                                        diff_val2);
            diff_val += diff_val2;
        }
        
        for (unsigned int i_dim=0; i_dim<dim; i_dim++)
        {
            // Galerkin contribution from the advection flux terms
            // linearized advection flux is obtained using the Ai dU  product
            Ai_advection[i_dim].vector_mult(flux, conservative_deltasol); // dF^adv_i
            dB_mat[i_dim].vector_mult_transpose(vec3_n2, flux); // dBw/dx_i dF^adv_i
            Fvec.add(JxW[qp], vec3_n2);
            
            
            if (_if_full_linearization)
            {
                // contribution from sensitivity of A matrix
                A_sens.zero();
                dB_mat[i_dim].vector_mult(vec2_n1, ref_sol);
                for (unsigned int i_cvar=0; i_cvar<n1; i_cvar++)
                {
                    flux_jacobian_sens[i_dim][i_cvar].vector_mult(vec1_n1,
                                                                  vec2_n1);
                    for (unsigned int i_phi=0; i_phi<n_phi; i_phi++)
                        A_sens.add_column((n_phi*i_cvar)+i_phi,
                                          phi[i_phi][qp], vec1_n1); // assuming that all variables have same n_phi
                }
                
                A_sens.vector_mult(flux, c.get_elem_solution());
                LS_mat.vector_mult_transpose(vec3_n2 , flux);
                Fvec.add(-JxW[qp], vec3_n2);// contribution from sensitivity of Ai Jacobians
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
        
        if (_if_full_linearization)
        {
            // sensitivity of LS term
            LS_sens.vector_mult(vec3_n2, c.get_elem_solution());
            Fvec.add(-JxW[qp], vec3_n2); // contribution from sensitivity of LS matrix
        }

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
                    dB_mat[i_dim].vector_mult(vec2_n1, ref_sol);
                    for (unsigned int i_cvar=0; i_cvar<n1; i_cvar++)
                    {
                        flux_jacobian_sens[i_dim][i_cvar].vector_mult
                        (vec1_n1, vec2_n1);
                        for (unsigned int i_phi=0; i_phi<n_phi; i_phi++)
                            A_sens.add_column((n_phi*i_cvar)+i_phi,
                                              phi[i_phi][qp], vec1_n1); // assuming that all variables have same n_phi
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
    
    //    std::cout << "inside element time derivative " << std::endl;
    //    c.elem->print_info();
    //    std::cout << "sol: " << std::endl; c.elem_solution.print(std::cout);
    //    std::cout << "res: " << std::endl; Fvec.print(std::cout);
    //    if (request_jacobian && c.elem_solution_derivative)
    //        Kmat.print(std::cout);
    
    return request_jacobian;
}



bool LinearizedFluidSystem::side_time_derivative (bool request_jacobian,
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
    DenseRealVector vec1_n2, flux, U_vec_interpolated, vec2_n1, ref_sol,
    conservative_sol, temp_grad, surface_steady_disp, surface_steady_vel,
    local_normal, duvec, uvec, conservative_deltasol, dnormal_steady,
    surface_unsteady_disp, surface_unsteady_vel, dnormal_unsteady;
    DenseRealMatrix  eig_val, l_eig_vec, l_eig_vec_inv_tr, mat_n1n1,
    mat1_n1n2, mat2_n2n2, mat3_n1n1, A_mat, dcons_dprim, dprim_dcons,
    stress_tensor;
    
    conservative_sol.resize(dim+2); conservative_deltasol.resize(dim+2);
    surface_unsteady_disp.resize(spatial_dim); surface_unsteady_vel.resize(spatial_dim);
    temp_grad.resize(dim); ref_sol.resize(n_dofs); dnormal_steady.resize(spatial_dim);
    vec1_n2.resize(n_dofs); flux.resize(n1); vec2_n1.resize(n1);
    U_vec_interpolated.resize(n1); dnormal_unsteady.resize(spatial_dim);
    uvec.resize(spatial_dim); duvec.resize(spatial_dim); surface_steady_disp.resize(spatial_dim);
    surface_steady_vel.resize(spatial_dim); local_normal.resize(spatial_dim);
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
    
    // reference solution for small-disturbance
    libMesh::System& fluid =
    this->get_equation_systems().get_system<libMesh::System>("FluidSystem");
    
    std::vector<libMesh::dof_id_type> fluid_dof_indices;
    fluid.get_dof_map().dof_indices(&c.get_elem(), fluid_dof_indices);
    
    for (unsigned int i=0; i<n_dofs; i++)
        ref_sol(i) = (*_local_fluid_solution)(fluid_dof_indices[i]);
    
    PrimitiveSolution p_sol;
    SmallPerturbationPrimitiveSolution<Real> delta_p_sol;
    
    switch (mechanical_bc_type)
    // adiabatic and isothermal are handled in the no-slip wall BC
    {
        case NO_SLIP_WALL: // only for viscous flows
                           // conditions enforced are
                           // vi ni = 0, vi = 0 = rho.vi      (no-slip wall, Dirichlet BC)
                           // thermal BC is handled here
        {
            
            Real dui_ni_unsteady = 0., ui_ni_steady = 0.;
            
            for (unsigned int qp=0; qp<qpoint.size(); qp++)
            {
                // first update the variables at the current quadrature point
                this->update_solution_at_quadrature_point(vars, qp, c,
                                                          false,
                                                          ref_sol,
                                                          conservative_sol,
                                                          p_sol,
                                                          B_mat,
                                                          dB_mat);
                
                // stress tensor and temperature gradient
                this->calculate_conservative_variable_jacobian(p_sol,
                                                               dcons_dprim,
                                                               dprim_dcons);
                
                this->calculate_diffusion_tensors(c.get_elem_solution(),
                                                  dB_mat,
                                                  dprim_dcons,
                                                  p_sol,
                                                  stress_tensor,
                                                  temp_grad);
                
                if (thermal_bc_type == ADIABATIC)
                    temp_grad.zero();
                
                
                // copy the surface normal
                for (unsigned int i_dim=0; i_dim<dim; i_dim++)
                    local_normal(i_dim) = face_normals[qp](i_dim);
                
                // now check if the surface deformation is defined and
                // needs to be applied through transpiration boundary
                // condition
                B_mat.vector_mult(conservative_deltasol, c.get_elem_solution());

                delta_p_sol.zero();
                delta_p_sol.init(p_sol, conservative_deltasol);

                dui_ni_unsteady = 0.; ui_ni_steady = 0.;
                p_sol.get_uvec(uvec);
                delta_p_sol.get_duvec(duvec);
                
                // calculate the surface velocity perturbations
                perturbed_surface_motion->surface_velocity(this->time,
                                                           qpoint[qp],
                                                           face_normals[qp],
                                                           surface_unsteady_disp,
                                                           surface_unsteady_vel,
                                                           dnormal_unsteady);
                
                if (surface_motion) // get the surface motion data
                {
                    surface_motion->surface_velocity(this->time,
                                                     qpoint[qp],
                                                     face_normals[qp],
                                                     surface_steady_disp,
                                                     surface_steady_vel,
                                                     dnormal_steady);
                    
                    // now initialize the surface velocity
                    // note that the perturbed local normal is used here
                    // since it resembles the normal of the surface that
                    // has undergone a static deformation, even though the
                    // surface_vel might be identically zero for a static
                    // body.
                    dui_ni_unsteady  =
                    surface_steady_vel.dot(dnormal_unsteady); // wi_dot * delta_ni
                    dui_ni_unsteady +=
                    surface_unsteady_vel.dot(dnormal_steady); // delta_wi_dot * dni
                    dui_ni_unsteady -= duvec.dot(dnormal_steady);
                    
                    // also the steady part
                    ui_ni_steady  = surface_steady_vel.dot(dnormal_steady);
                    ui_ni_steady -= uvec.dot(dnormal_steady);
                    
                    for (unsigned int i_dim=0; i_dim<dim; i_dim++)
                        ui_ni_steady    += // wi_dot * ni
                        surface_steady_vel(i_dim) * face_normals[qp](i_dim);
                }
                
                // now add the contribution from unsteady normal perturbation
                for (unsigned int i_dim=0; i_dim<dim; i_dim++)
                    dui_ni_unsteady += // delta_wi_dot * ni
                    surface_unsteady_vel(i_dim) * face_normals[qp](i_dim);
                dui_ni_unsteady -= uvec.dot(dnormal_unsteady); // ui delta_ni
                
                /*// add the contribution from divergence of ui_ni
                 for (unsigned int i_dim=0; i_dim<dim; i_dim++) {
                 dB_mat[i_dim].vector_mult(vec2_n1, ref_sol); // dU/dx_i
                 dprim_dcons.vector_mult(vec3_n1, vec2_n1);  // dU_primitive / dx_i
                 for (unsigned int j_dim=0; j_dim<dim; j_dim++)
                 dui_ni_unsteady -= vec3_n1(j_dim+1) * surface_unsteady_disp(i_dim) *
                 (face_normals[qp](j_dim) + std::real(dnormal_steady(j_dim)));
                 }*/
                
                // now prepare the flux vector
                flux.zero();
                flux.add(dui_ni_unsteady, conservative_sol);
                flux.add(ui_ni_steady, conservative_deltasol);
                flux(n1-1) += dui_ni_unsteady*p_sol.p + ui_ni_steady*delta_p_sol.dp;
                for (unsigned int i_dim=0; i_dim<dim; i_dim++)
                    flux(i_dim+1) += delta_p_sol.dp * face_normals[qp](i_dim);

                B_mat.vector_mult_transpose(vec1_n2, flux);
                Fvec.add(JxW[qp], vec1_n2);
                
                
                if ( request_jacobian)
                {
                    this->calculate_advection_flux_jacobian_for_moving_solid_wall_boundary
                    (p_sol, ui_ni_steady, face_normals[qp], dnormal_steady, A_mat);
                    
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
            dnormal_steady.zero();

            for (unsigned int qp=0; qp<qpoint.size(); qp++)
            {
                // first update the variables at the current quadrature point
                this->update_solution_at_quadrature_point(vars, qp, c,
                                                          false,
                                                          ref_sol,
                                                          conservative_sol,
                                                          p_sol,
                                                          B_mat,
                                                          dB_mat);
                
                // initialize flux to interpolated sol for initialization
                // of perturbed vars
                B_mat.vector_mult(conservative_deltasol, c.get_elem_solution());
                
                delta_p_sol.zero();
                delta_p_sol.init(p_sol, conservative_deltasol);

                
                flux.zero();
                // since vi=0, vi ni = 0, so the advection flux gets evaluated
                // using the slip wall condition
                for (unsigned int i_dim=0; i_dim<dim; i_dim++)
                    flux(i_dim+1) += delta_p_sol.dp * face_normals[qp](i_dim);
                
                B_mat.vector_mult_transpose(vec1_n2, flux);
                Fvec.add(-JxW[qp], vec1_n2);
                
                if ( request_jacobian)
                {
                    this->calculate_advection_flux_jacobian_for_moving_solid_wall_boundary
                    (p_sol, 0., face_normals[qp], dnormal_steady, A_mat);
                    
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
            Real dui_ni_unsteady = 0., ui_ni_steady = 0.;
            
            for (unsigned int qp=0; qp<qpoint.size(); qp++)
            {
                // first update the variables at the current quadrature point
                this->update_solution_at_quadrature_point(vars, qp, c,
                                                          false,
                                                          ref_sol,
                                                          conservative_sol,
                                                          p_sol,
                                                          B_mat,
                                                          dB_mat);
                
                // copy the surface normal
                for (unsigned int i_dim=0; i_dim<dim; i_dim++)
                    local_normal(i_dim) = face_normals[qp](i_dim);
                
                // calculate the conservative flow variable Jacobians
                this->calculate_conservative_variable_jacobian
                (p_sol, dcons_dprim, dprim_dcons);

                // initialize flux to interpolated sol for initialization
                // of perturbed vars
                B_mat.vector_mult(conservative_deltasol, c.get_elem_solution());
                
                delta_p_sol.zero();
                delta_p_sol.init(p_sol, conservative_deltasol);

                // now check if the surface deformation is defined and
                // needs to be applied through transpiration boundary
                // condition
                dui_ni_unsteady = 0.; ui_ni_steady = 0.;
                p_sol.get_uvec(uvec);
                delta_p_sol.get_duvec(duvec);

                // calculate the surface velocity perturbations
                if (perturbed_surface_motion)
                    perturbed_surface_motion->surface_velocity(this->time,
                                                               qpoint[qp],
                                                               face_normals[qp],
                                                               surface_unsteady_disp,
                                                               surface_unsteady_vel,
                                                               dnormal_unsteady);

                if (surface_motion) // get the surface motion data
                {
                    surface_motion->surface_velocity(this->time,
                                                     qpoint[qp],
                                                     face_normals[qp],
                                                     surface_steady_disp,
                                                     surface_steady_vel,
                                                     dnormal_steady);
                    
                    // now initialize the surface velocity
                    // note that the perturbed local normal is used here
                    // since it resembles the normal of the surface that
                    // has undergone a static deformation, even though the
                    // surface_vel might be identically zero for a static
                    // body.
                    dui_ni_unsteady  =
                    surface_steady_vel.dot(dnormal_unsteady); // wi_dot * delta_ni
                    dui_ni_unsteady +=
                    surface_unsteady_vel.dot(dnormal_steady); // delta_wi_dot * dni
                    dui_ni_unsteady -= duvec.dot(dnormal_steady);
                    
                    // also the steady part
                    ui_ni_steady  = surface_steady_vel.dot(dnormal_steady);
                    ui_ni_steady -= uvec.dot(dnormal_steady);
                    
                    for (unsigned int i_dim=0; i_dim<dim; i_dim++)
                        ui_ni_steady    += // wi_dot * ni
                        surface_steady_vel(i_dim) * face_normals[qp](i_dim);
                }
                
                // now add the contribution from unsteady normal perturbation
                for (unsigned int i_dim=0; i_dim<dim; i_dim++) {
                    dui_ni_unsteady += // delta_wi_dot * ni
                    surface_unsteady_vel(i_dim) * face_normals[qp](i_dim);
                    dui_ni_unsteady -= uvec(i_dim) * face_normals[qp](i_dim);
                }
                dui_ni_unsteady -= uvec.dot(dnormal_unsteady); // ui delta_ni
                
                /*// add the contribution from divergence of ui_ni
                 for (unsigned int i_dim=0; i_dim<dim; i_dim++) {
                 dB_mat[i_dim].vector_mult(vec2_n1, ref_sol); // dU/dx_i
                 dprim_dcons.vector_mult(vec3_n1, vec2_n1);  // dU_primitive / dx_i
                 for (unsigned int j_dim=0; j_dim<dim; j_dim++)
                 dui_ni_unsteady -= vec3_n1(j_dim+1) * surface_unsteady_disp(i_dim) *
                 (face_normals[qp](j_dim) + std::real(dnormal_steady(j_dim)));
                 }*/
                
                // now prepare the flux vector
                flux.zero();
                flux.add(dui_ni_unsteady, conservative_sol);
                flux.add(ui_ni_steady, conservative_deltasol);
                flux(n1-1) += dui_ni_unsteady*p_sol.p + ui_ni_steady*delta_p_sol.dp;
                for (unsigned int i_dim=0; i_dim<dim; i_dim++)
                    flux(i_dim+1) += delta_p_sol.dp * face_normals[qp](i_dim);
                
                B_mat.vector_mult_transpose(vec1_n2, flux);
                Fvec.add(-JxW[qp], vec1_n2);

                if ( request_jacobian )
                {
                    this->calculate_advection_flux_jacobian_for_moving_solid_wall_boundary
                    (p_sol, ui_ni_steady, face_normals[qp], dnormal_steady, A_mat);
                    
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
                this->update_solution_at_quadrature_point(vars, qp, c,
                                                          false,
                                                          ref_sol,
                                                          conservative_sol,
                                                          p_sol,
                                                          B_mat, dB_mat);
                
                this->calculate_advection_left_eigenvector_and_inverse_for_normal
                (p_sol, face_normals[qp], eig_val,
                 l_eig_vec, l_eig_vec_inv_tr);
                
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
                
                B_mat.vector_mult(conservative_deltasol, c.get_elem_solution()); // B dU
                
                mat_n1n1.vector_mult(flux, conservative_deltasol); // f_{+} = A_{+} B U
                
                B_mat.vector_mult_transpose(vec1_n2, flux); // B^T f_{+}   (this is flux going out of the solution domain)
                Fvec.add(-JxW[qp], vec1_n2);
                
                
                if ( request_jacobian )
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
            libmesh_error(); // to be implemented
        }
            break;
            
    }
    
    
    //    std::cout << "inside side constraint " << std::endl;
    //    std::cout << "elem solution" << std::endl; c.elem_solution.print(std::cout);
    //    std::cout << if_inf_bc << "  " << if_wall_bc << std::endl;
    //    std::cout << "bc vec: " << std::endl; Fvec.print(std::cout);
    //    if (request_jacobian && c.elem_solution_derivative)
    //        Kmat.print(std::cout);
    
    return request_jacobian;
}





bool LinearizedFluidSystem::mass_residual (bool request_jacobian,
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
    DenseRealVector flux, vec1_n1, vec3_n2, conservative_sol, ref_sol;
    LS_mat.resize(n1, n_dofs); LS_sens.resize(n_dofs, n_dofs);
    mat2_n2n2.resize(n_dofs, n_dofs);
    Ai_Bi_advection.resize(dim+2, n_dofs);  mat_n1n2.resize(dim+2, n_dofs);
    flux.resize(n1); vec1_n1.resize(n1); vec3_n2.resize(n_dofs);
    conservative_sol.resize(dim+2);
    ref_sol.resize(n_dofs);
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
    
    // reference solution for small-disturbance
    libMesh::System& fluid =
    this->get_equation_systems().get_system<libMesh::System>("FluidSystem");
    
    std::vector<libMesh::dof_id_type> fluid_dof_indices;
    fluid.get_dof_map().dof_indices(&c.get_elem(), fluid_dof_indices);
    
    for (unsigned int i=0; i<n_dofs; i++)
        ref_sol(i) = (*_local_fluid_solution)(fluid_dof_indices[i]);
    
    PrimitiveSolution primitive_sol;
    
    for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
        // first update the variables at the current quadrature point
        this->update_solution_at_quadrature_point(vars, qp, c,
                                                  true, ref_sol,
                                                  conservative_sol,
                                                  primitive_sol,
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
            (vars, qp, c, ref_sol,
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
        
        
        if (request_jacobian)
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
    
    //    std::cout << "inside mass residual " << std::endl;
    //    std::cout << "elem velocity" << std::endl; c.elem_solution.print(std::cout);
    //    std::cout << "elem solution" << std::endl; c.elem_fixed_solution.print(std::cout);
    //    std::cout << "mass vec: " << std::endl; Fvec.print(std::cout);
    //    if (request_jacobian && c.elem_solution_derivative)
    //        Kmat.print(std::cout);
    
    return request_jacobian;
}



void LinearizedFluidSystem::postprocess()
{ }





void LinearizedFluidSystem::evaluate_recalculate_dc_flag()
{
    Real norm = this->calculate_norm(*(this->solution), this->vars[0], libMesh::L2),
    relative_change0 = fabs(_rho_norm_curr - _rho_norm_old)/_rho_norm_curr,
    relative_change1 = fabs(norm - _rho_norm_curr)/norm;
    
    libMesh::out
    << "Rho libMesh::L2-norm delta: " << relative_change0
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





