//
//  frequency_domain_linearized_euler.cpp
//  FESystem
//
//  Created by Manav Bhatia on 3/8/13.
//
//

// FEsystem includes
#include "FluidElems/frequency_domain_linearized_fluid_system.h"


#ifdef LIBMESH_USE_COMPLEX_NUMBERS

// FESystem includes
#include "FluidElems/fluid_system.h"
#include "Numerics/fem_operator_matrix.h"

// Basic include files
#include "libmesh/equation_systems.h"
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/parameters.h"
#include "libmesh/quadrature.h"
#include "libmesh/mesh_function.h"
#include "libmesh/fem_function_base.h"
#include "libmesh/dof_map.h"



void FrequencyDomainLinearizedFluidSystem::init_data()
{
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
    
    vars[0]  = this->add_variable ( "drho", static_cast<Order>(o), fefamily);
    this->time_evolving(vars[0]);
    params.set<Real> ("rho_inf") = flight_condition->rho();
    
    vars[1] = this->add_variable ("drhoux", static_cast<Order>(o), fefamily);
    this->time_evolving(vars[1]);
    params.set<Real> ("rhoux_inf") = flight_condition->rho_u1();
    
    if (dim > 1)
    {
        vars[2] = this->add_variable ("drhouy", static_cast<Order>(o), fefamily);
        this->time_evolving(vars[2]);
        params.set<Real> ("rhouy_inf") = flight_condition->rho_u2();
    }
    
    if (dim > 2)
    {
        vars[3] = this->add_variable ("drhouz", static_cast<Order>(o), fefamily);
        this->time_evolving(vars[3]);
        params.set<Real> ("rhouz_inf") = flight_condition->rho_u3();
    }
    
    vars[dim+2-1] = this->add_variable ("drhoe", static_cast<Order>(o), fefamily);
    this->time_evolving(vars[dim+2-1]);
    params.set<Real> ("rhoe_inf") = flight_condition->rho_e();
    
    // Useful debugging options
    // Set verify_analytic_jacobians to 1e-6 to use
    this->verify_analytic_jacobians = infile("verify_analytic_jacobians", 0.);
    this->print_jacobians = infile("print_jacobians", false);
    this->print_element_jacobians = infile("print_element_jacobians", false);
    
    // Do the parent's initialization after variables and boundary constraints are defined
    FEMSystem::init_data();
}



void FrequencyDomainLinearizedFluidSystem::init_context(libMesh::DiffContext &context)
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
        if (_if_viscous)
            elem_fe[i]->get_d2phi();
    }
    
    std::vector<FEBase*> elem_side_fe(dim+2);
    
    for (unsigned int i=0; i<dim+2; i++)
    {
        c.get_side_fe( vars[i], elem_side_fe[i]);
        elem_side_fe[i]->get_JxW();
        elem_side_fe[i]->get_phi();
        elem_side_fe[i]->get_xyz();
        if (_if_viscous)
            elem_side_fe[i]->get_dphi();
    }
}



void FrequencyDomainLinearizedFluidSystem::localize_fluid_solution()
{
    libmesh_assert(!_if_localized_sol);
    
    _local_fluid_solution =
    NumericVector<Number>::build(this->get_equation_systems().comm());
    
    
    System& fluid =
    this->get_equation_systems().get_system<System>("EulerSystem");

    _local_fluid_solution->init(fluid.solution->size(), true, SERIAL);
    fluid.solution->localize(*_local_fluid_solution,
                             fluid.get_dof_map().get_send_list());
    
    _if_localized_sol = true;
}




bool FrequencyDomainLinearizedFluidSystem::element_time_derivative
(bool request_jacobian, DiffContext &context)
{
    libmesh_assert(_if_localized_sol);
    
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
    const unsigned int n_qpoints = c.get_element_qrule().n_points(), n1 = dim+2;
    
    std::vector<FEMOperatorMatrix> dB_mat(dim);
    std::vector<DenseMatrix<Real> > Ai_advection(dim);
    FEMOperatorMatrix B_mat;
    DenseMatrix<Real> LS_mat, LS_sens, Ai_Bi_advection, A_entropy,
    A_inv_entropy, tmp_mat_n1n2, tmp_mat_n2n2,
    tmp_mat_n2n1, tmp_mat_n1n1, A_sens, stress_tensor, dcons_dprim, dprim_dcons;
    DenseVector<Real> tmp_vec1_n1, tmp_vec2_n1, conservative_sol, ref_sol,
    temp_grad, flux_real, tmp_vec3_n2_real;
    DenseVector<Number> elem_interpolated_sol, flux, tmp_vec3_n2, tmp_vec4_n1,
    tmp_vec5_n1;
    
    LS_mat.resize(n1, n_dofs); LS_sens.resize(n_dofs, n_dofs),
    Ai_Bi_advection.resize(dim+2, n_dofs); A_sens.resize(n1, n_dofs);
    stress_tensor.resize(dim, dim);
    A_inv_entropy.resize(dim+2, dim+2); A_entropy.resize(dim+2, dim+2);
    flux.resize(n1); tmp_vec1_n1.resize(n1); tmp_vec2_n1.resize(n1);
    tmp_vec3_n2.resize(n_dofs);
    conservative_sol.resize(dim+2); temp_grad.resize(dim); tmp_vec4_n1.resize(dim+2);
    tmp_vec5_n1.resize(dim+2); tmp_mat_n1n2.resize(dim+2, n_dofs);
    tmp_mat_n2n2.resize(n_dofs, n_dofs);
    tmp_mat_n2n1.resize(n_dofs, dim+2); tmp_mat_n1n1.resize(dim+2, dim+2);
    dcons_dprim.resize(dim+2, dim+2); dprim_dcons.resize(dim+2, dim+2);
    elem_interpolated_sol.resize(n1); ref_sol.resize(n_dofs);
    tmp_vec3_n2_real.resize(n_dofs);
    
    for (unsigned int i=0; i<dim; i++)
        Ai_advection[i].resize(dim+2, dim+2);
    
    std::vector<std::vector<DenseMatrix<Real> > > flux_jacobian_sens;
    flux_jacobian_sens.resize(dim);
    for (unsigned int i_dim=0; i_dim<dim; i_dim++)
    {
        flux_jacobian_sens[i_dim].resize(n1); // number of variables for sensitivity
        for (unsigned int i_cvar=0; i_cvar<n1; i_cvar++)
            flux_jacobian_sens[i_dim][i_cvar].resize(n1, n1);
    }
    
    Real diff_val=0.;
    
    PrimitiveSolution primitive_sol;

    // element dofs from steady solution to calculate the linearized quantities
    System& fluid =
    this->get_equation_systems().get_system<System>("EulerSystem");
    
    std::vector<dof_id_type> fluid_dof_indices;
    fluid.get_dof_map().dof_indices(&c.get_elem(), fluid_dof_indices);
    
    for (unsigned int i=0; i<c.get_dof_indices().size(); i++)
        ref_sol(i) = real((*_local_fluid_solution)(fluid_dof_indices[i]));
    
    Number iota(0, 1.), scaling = iota*perturbed_surface_motion->frequency;

    System& delta_val_system =
    this->get_equation_systems().get_system<System>("DeltaValSystem");
    NumericVector<Number>& diff_val_vec = (*delta_val_system.solution.get());
    
    const std::vector<std::vector<Real> >& phi = elem_fe->get_phi(); // assuming that all variables have the same interpolation
    const unsigned int n_phi = phi.size();
    

    for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
        // first update the variables at the current quadrature point
        this->update_solution_at_quadrature_point(vars, qp, c, true, ref_sol,
                                                  conservative_sol,
                                                  primitive_sol,
                                                  B_mat, dB_mat);

        for ( unsigned int i_dim=0; i_dim < dim; i_dim++)
            this->calculate_advection_flux_jacobian( i_dim, primitive_sol,
                                                    Ai_advection[i_dim] );
        
        this->calculate_entropy_variable_jacobian ( primitive_sol, A_entropy,
                                                   A_inv_entropy );

        Ai_Bi_advection.zero();
        
        for (unsigned int i_dim=0; i_dim<dim; i_dim++)
        {
            dB_mat[i_dim].left_multiply(tmp_mat_n1n2, Ai_advection[i_dim]);
            Ai_Bi_advection.add(1.0, tmp_mat_n1n2);
        }
        
        if (_if_full_linearization)
        {
            for (unsigned int i_dim=0; i_dim<dim; i_dim++)
                this->calculate_advection_flux_jacobian_sensitivity_for_conservative_variable
                (i_dim, primitive_sol, flux_jacobian_sens[i_dim]);
        }
        
        if (_if_update_stabilization_per_quadrature_point || (qp == 0))
            this->calculate_differential_operator_matrix
            (vars, qp, c, ref_sol, primitive_sol, B_mat,
             dB_mat, Ai_advection, Ai_Bi_advection,
             A_inv_entropy, flux_jacobian_sens,
             LS_mat, LS_sens, diff_val);
        
        // the discontinuity capturing term is reused from the steady solution
        diff_val = real(diff_val_vec.el(c.get_elem().dof_number
                                        (delta_val_system.number(), 0, 0)));

        // calculate the interpolated solution value at this quadrature point
        B_mat.vector_mult(elem_interpolated_sol, c.get_elem_solution());  // B dU

        // Galerkin contribution to solution
        B_mat.vector_mult_transpose(tmp_vec3_n2, elem_interpolated_sol); // B^T B dU
        Fvec.add(JxW[qp]*scaling, tmp_vec3_n2); // mass term
        
        // Galerkin contribution to solution
        LS_mat.vector_mult_transpose(tmp_vec3_n2, elem_interpolated_sol); // LS^T tau B dU
        Fvec.add(JxW[qp]*scaling, tmp_vec3_n2); // mass term

        
        
        for (unsigned int i_dim=0; i_dim<dim; i_dim++)
        {
            // Galerkin contribution from the advection flux terms
            // linearized advection flux is obtained using the Ai dU  product
            Ai_advection[i_dim].vector_mult(flux, elem_interpolated_sol); // dF^adv_i
            dB_mat[i_dim].vector_mult_transpose(tmp_vec3_n2, flux); // dBw/dx_i dF^adv_i
            Fvec.add(-JxW[qp], tmp_vec3_n2);

            if (_if_full_linearization)
            {
                // contribution from sensitivity of A matrix
                A_sens.zero();
                dB_mat[i_dim].vector_mult(tmp_vec2_n1, ref_sol);
                for (unsigned int i_cvar=0; i_cvar<n1; i_cvar++)
                {
                    flux_jacobian_sens[i_dim][i_cvar].vector_mult(tmp_vec1_n1,
                                                                  tmp_vec2_n1);
                    for (unsigned int i_phi=0; i_phi<n_phi; i_phi++)
                        A_sens.add_column((n_phi*i_cvar)+i_phi,
                                          phi[i_phi][qp], tmp_vec1_n1); // assuming that all variables have same n_phi
                }
                A_sens.vector_mult(flux, c.get_elem_solution());
                LS_mat.vector_mult_transpose(tmp_vec3_n2 , flux);
                Fvec.add(JxW[qp], tmp_vec3_n2); // contribution from sensitivity of Ai Jacobians
            }
            
            if (_if_viscous)
            {
                // Galerkin contribution from the diffusion flux terms
                for (unsigned int deriv_dim=0; deriv_dim<dim; deriv_dim++)
                {
                    dB_mat[deriv_dim].vector_mult(tmp_vec4_n1,
                                                  c.get_elem_solution()); // dU/dx_j
                    this->calculate_diffusion_flux_jacobian
                    (i_dim, deriv_dim, primitive_sol, tmp_mat_n1n1); // Kij
                    tmp_mat_n1n1.vector_mult(tmp_vec5_n1, tmp_vec4_n1); // Kij dB/dx_j
                    dB_mat[i_dim].vector_mult_transpose
                    (tmp_vec3_n2, tmp_vec5_n1); // dB/dx_i Kij dB/dx_j
                    Fvec.add(JxW[qp], tmp_vec3_n2);
                }
            }

            // discontinuity capturing operator
            dB_mat[i_dim].vector_mult(flux, c.get_elem_solution());  // d dU/ dx_i
            dB_mat[i_dim].vector_mult_transpose(tmp_vec3_n2, flux); // dB/dx_i d dU/ dx_i
            Fvec.add(JxW[qp]*diff_val, tmp_vec3_n2);
        }
        
        // Least square contribution from flux
        Ai_Bi_advection.vector_mult(flux, c.get_elem_solution()); // d dF^adv_i / dxi
        LS_mat.vector_mult_transpose(tmp_vec3_n2, flux); // LS^T tau dF^adv_i
        Fvec.add(JxW[qp], tmp_vec3_n2);

        if (_if_full_linearization)
        {
            // sensitivity of LS term
            LS_sens.vector_mult(tmp_vec3_n2, c.get_elem_solution());
            Fvec.add(JxW[qp], tmp_vec3_n2); // contribution from sensitivity of LS matrix
        }

        // Least square contribution from divergence of diffusion flux
        // TODO: this requires a 2nd order differential of the flux

        if (request_jacobian && c.get_elem_solution_derivative())
        {
            // contribution from unsteady term
            // Galerkin contribution of velocity
            B_mat.right_multiply_transpose(tmp_mat_n2n2, B_mat); // B^T B
            Kmat.add(JxW[qp]*scaling, tmp_mat_n2n2); // mass term
            
            // LS contribution of velocity
            LS_mat.get_transpose(tmp_mat_n2n1);
            B_mat.left_multiply(tmp_mat_n2n2, tmp_mat_n2n1); // LS^T tau Bmat
            Kmat.add(JxW[qp]*scaling, tmp_mat_n2n2); // mass term

            
            A_sens.zero();

            for (unsigned int i_dim=0; i_dim<dim; i_dim++)
            {
                // Galerkin contribution from the advection flux terms
                B_mat.left_multiply(tmp_mat_n1n2, Ai_advection[i_dim]);
                dB_mat[i_dim].right_multiply_transpose
                (tmp_mat_n2n2, tmp_mat_n1n2); // dBw/dx_i^T  d dF^adv_i/ dU
                Kmat.add(-JxW[qp], tmp_mat_n2n2);
                
                // discontinuity capturing term
                dB_mat[i_dim].right_multiply_transpose(tmp_mat_n2n2,
                                                       dB_mat[i_dim]);
                Kmat.add(JxW[qp]*diff_val, tmp_mat_n2n2);

                if (_if_full_linearization)
                {
                    // contribution from sensitivity of A matrix
                    dB_mat[i_dim].vector_mult(tmp_vec2_n1, ref_sol);
                    for (unsigned int i_cvar=0; i_cvar<n1; i_cvar++)
                    {
                        flux_jacobian_sens[i_dim][i_cvar].vector_mult
                        (tmp_vec1_n1, tmp_vec2_n1);
                        for (unsigned int i_phi=0; i_phi<n_phi; i_phi++)
                            A_sens.add_column((n_phi*i_cvar)+i_phi,
                                              phi[i_phi][qp], tmp_vec1_n1); // assuming that all variables have same n_phi
                    }
                }

                if (_if_viscous)
                {
                    for (unsigned int deriv_dim=0; deriv_dim<dim; deriv_dim++)
                    {
                        this->calculate_diffusion_flux_jacobian
                        (i_dim, deriv_dim, primitive_sol, tmp_mat_n1n1); // Kij
                        dB_mat[deriv_dim].left_multiply(tmp_mat_n1n2,
                                                        tmp_mat_n1n1); // Kij dB/dx_j
                        dB_mat[i_dim].right_multiply_transpose
                        (tmp_mat_n2n2, tmp_mat_n1n2); // dB/dx_i^T Kij dB/dx_j
                        Kmat.add(JxW[qp], tmp_mat_n2n2);
                    }
                }
            }
            
            // Lease square contribution of flux gradient
            LS_mat.get_transpose(tmp_mat_n2n1);
            tmp_mat_n2n1.right_multiply(Ai_Bi_advection); // LS^T tau d^2 dF^adv_i / dx dU
            Kmat.add(JxW[qp], tmp_mat_n2n1);

            if (_if_full_linearization)
            {
                LS_mat.get_transpose(tmp_mat_n2n1);
                tmp_mat_n2n1.right_multiply(A_sens); // LS^T tau d^2F^adv_i / dx dU  (Ai sensitivity)
                Kmat.add(JxW[qp], tmp_mat_n2n1);
                
                // contribution from sensitivity of the LS.tau matrix
                Kmat.add(JxW[qp], LS_sens);
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



bool FrequencyDomainLinearizedFluidSystem::side_time_derivative
(bool request_jacobian, DiffContext &context)
{
    libmesh_assert(_if_localized_sol);

    FEMContext &c = libmesh_cast_ref<FEMContext&>(context);
    
    // check for the boundary tags to check if the boundary condition needs to be applied to this element
    
    FluidBoundaryConditionType mechanical_bc_type, thermal_bc_type;
    
    unsigned int n_mechanical_bc = 0, n_thermal_bc = 0;
    
    std::multimap<unsigned int, FluidBoundaryConditionType>::const_iterator
    bc_it     = this->_boundary_condition.begin(),
    bc_end    = this->_boundary_condition.end();
    
    
    for ( ; bc_it != bc_end; bc_it++)
    {
        if ( this->get_mesh().boundary_info->has_boundary_id
            (&c.get_elem(), c.side, bc_it->first))
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
    
    
    const unsigned int n1 = dim+2;
    
    // The number of local degrees of freedom for this element
    unsigned int n_dofs = 0;
    for (unsigned int i=0; i<dim+2; i++)
        n_dofs += c.get_dof_indices( vars[i] ).size();
    
    libmesh_assert_equal_to (n_dofs, (dim+2)*c.get_dof_indices( vars[0] ).size());
    
    // The subvectors and submatrices we need to fill:
    DenseMatrix<Number>& Kmat = c.get_elem_jacobian();
    DenseVector<Number>& Fvec = c.get_elem_residual();
    
    FEBase * side_fe;
    c.get_side_fe(vars[0], side_fe);
    
    // Element Jacobian * quadrature weights for interior integration
    const std::vector<Real> &JxW = side_fe->get_JxW();
    
    // Physical location of the quadrature points
    const std::vector<Point>& qpoint = side_fe->get_xyz();
    
    // boundary normals
    const std::vector<Point>& face_normals = side_fe->get_normals();
    
    
    Point vel; vel.zero(); // zero surface velocity
    

    std::vector<FEMOperatorMatrix> dB_mat(dim);
    FEMOperatorMatrix B_mat;

    DenseMatrix<Real>  eig_val, l_eig_vec, l_eig_vec_inv_tr, tmp_mat_n1n1,
    tmp_mat_n1n2, tmp_mat_n2n2, A_mat,
    dcons_dprim, dprim_dcons, stress_tensor, Kmat_viscous;
    DenseVector<Real>  normal, normal_local, tmp_vec1, U_vec_interpolated,
    conservative_sol, ref_sol, temp_grad, uvec, dnormal_steady_real;
    DenseVector<Number> elem_interpolated_sol, flux, tmp_vec1_n2,
    surface_unsteady_vel, dnormal_unsteady, duvec, conservative_deltasol,
    dnormal_steady, surface_steady_vel;
    
    conservative_sol.resize(dim+2);
    normal.resize(3); normal_local.resize(dim); tmp_vec1.resize(n_dofs);
    flux.resize(n1); tmp_vec1_n2.resize(n_dofs);
    conservative_deltasol.resize(dim+2);
    U_vec_interpolated.resize(n1); temp_grad.resize(dim);
    elem_interpolated_sol.resize(n1); ref_sol.resize(n_dofs);
    dnormal_steady.resize(3); dnormal_steady_real.resize(3);
    surface_steady_vel.resize(3); uvec.resize(3); duvec.resize(3);
    dnormal_unsteady.resize(3); surface_unsteady_vel.resize(3);
    
    
    eig_val.resize(n1, n1); l_eig_vec.resize(n1, n1);
    l_eig_vec_inv_tr.resize(n1, n1);
    tmp_mat_n1n1.resize(n1, n1); tmp_mat_n1n2.resize(n1, n_dofs);
    tmp_mat_n2n2.resize(n_dofs, n_dofs);
    dcons_dprim.resize(n1, n1); dprim_dcons.resize(n1, n1);
    Kmat_viscous.resize(n1, n_dofs);
    A_mat.resize(dim+2, dim+2); stress_tensor.resize(dim, dim);


    
    // element dofs from steady solution to calculate the linearized quantities
    System& fluid =
    this->get_equation_systems().get_system<System>("EulerSystem");
    
    std::vector<dof_id_type> fluid_dof_indices;
    fluid.get_dof_map().dof_indices(&c.get_elem(), fluid_dof_indices);
    
    for (unsigned int i=0; i<c.get_dof_indices().size(); i++)
        ref_sol(i) = real((*_local_fluid_solution)(fluid_dof_indices[i]));
    
    
    PrimitiveSolution p_sol;
    SmallPerturbationPrimitiveSolution<Number> delta_p_sol;
    
    const DenseVector<Number>& elem_sol = c.get_elem_solution();
    
    switch (mechanical_bc_type)
    // adiabatic and isothermal are handled in the no-slip wall BC
    {
        case NO_SLIP_WALL: // only for viscous flows
            // conditions enforced are
            // vi ni = 0, vi = 0 = rho.vi      (no-slip wall, Dirichlet BC)
            // thermal BC is handled here
        {
            
            Number dui_ni_unsteady = 0., ui_ni_steady = 0.;
            
            for (unsigned int qp=0; qp<qpoint.size(); qp++)
            {
                // first update the variables at the current quadrature point
                this->update_solution_at_quadrature_point
                (vars, qp, c, false,
                 ref_sol, conservative_sol,
                 p_sol, B_mat, dB_mat);
                
                // stress tensor and temperature gradient
                this->calculate_conservative_variable_jacobian
                (p_sol, dcons_dprim, dprim_dcons);
                this->calculate_diffusion_tensors
                (ref_sol, dB_mat, dprim_dcons,
                 p_sol, stress_tensor, temp_grad);
                
                if (thermal_bc_type == ADIABATIC)
                    temp_grad.zero();
                
                // initialize flux to interpolated sol for initialization
                // of perturbed vars
                B_mat.vector_mult(conservative_deltasol, elem_sol);
                
                delta_p_sol.zero();
                delta_p_sol.init(p_sol, conservative_deltasol);

                dui_ni_unsteady = 0.; ui_ni_steady = 0.;
                p_sol.get_uvec(uvec);
                delta_p_sol.get_duvec(duvec);
                
                // calculate the surface velocity perturbations
                perturbed_surface_motion->surface_velocity_frequency_domain
                (qpoint[qp], face_normals[qp],
                 surface_unsteady_vel, dnormal_unsteady);

                // now check if the surface deformation is defined and
                // needs to be applied through transpiration boundary
                // condition
                if (surface_motion) // get the surface motion data
                {
                    surface_motion->surface_velocity_time_domain
                    (0., qpoint[qp], face_normals[qp],
                     surface_steady_vel, dnormal_steady);
                    
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
                    {
                        dui_ni_unsteady += // delta_wi_dot * ni
                        surface_unsteady_vel(i_dim) * face_normals[qp](i_dim);
                        ui_ni_steady    += // wi_dot * ni
                        surface_steady_vel(i_dim) * face_normals[qp](i_dim);
                    }
                }


                // now add the contribution from unsteady normal perturbation
                dui_ni_unsteady -= uvec.dot(dnormal_unsteady); // ui delta_ni                
                
                // now prepare the flux vector
                flux.zero();
                flux.add(dui_ni_unsteady, conservative_sol);
                flux.add(ui_ni_steady, conservative_deltasol);
                flux(n1-1) += dui_ni_unsteady*p_sol.p + ui_ni_steady*delta_p_sol.dp;
                for (unsigned int i_dim=0; i_dim<dim; i_dim++)
                    flux(i_dim+1) += delta_p_sol.dp * face_normals[qp](i_dim);

                B_mat.vector_mult_transpose(tmp_vec1_n2, flux); // inviscid flux contribution
                Fvec.add(JxW[qp], tmp_vec1_n2);
                
                // for viscous flux contribution, multiply the viscous flux Jacobian with the
                // solution, and add it to the force vector
                Kmat_viscous.zero();
                for (unsigned int i_dim=0; i_dim<dim; i_dim++)
                    for (unsigned int deriv_dim=0; deriv_dim<dim; deriv_dim++)
                    {
                        this->calculate_diffusion_flux_jacobian
                        (i_dim, deriv_dim, p_sol, tmp_mat_n1n1); // Kij
                        dB_mat[deriv_dim].left_multiply
                        (tmp_mat_n1n2, tmp_mat_n1n1); // Kij dB/dx_j
                        Kmat_viscous.add(1.0, tmp_mat_n1n2);
                    }
                
                Kmat_viscous.vector_mult(flux, elem_sol);  // Kij dB/dx_j dU
                B_mat.vector_mult_transpose(tmp_vec1_n2, flux);   // B^T Kij dB/dx_j dU
                Fvec.add(-JxW[qp], tmp_vec1_n2);

                
                if ( request_jacobian && c.get_elem_solution_derivative() )
                {
                    for (unsigned int i=0; i<dim; i++)
                        dnormal_steady_real(i) = std::real(dnormal_steady(i));
                    
                    this->calculate_advection_flux_jacobian_for_moving_solid_wall_boundary
                    (p_sol, std::real(ui_ni_steady), face_normals[qp],
                     dnormal_steady_real, A_mat);
                    
                    B_mat.left_multiply(tmp_mat_n1n2, A_mat); // A B
                    B_mat.right_multiply_transpose(tmp_mat_n2n2, tmp_mat_n1n2); // B^T A B
                    Kmat.add(JxW[qp], tmp_mat_n2n2);
                    
                    // contribution from diffusion flux
                    for (unsigned int i_dim=0; i_dim<dim; i_dim++)
                        for (unsigned int deriv_dim=0;
                             deriv_dim<dim; deriv_dim++)
                        {
                            this->calculate_diffusion_flux_jacobian
                            (i_dim, deriv_dim, p_sol, tmp_mat_n1n1); // Kij
                            dB_mat[deriv_dim].left_multiply
                            (tmp_mat_n1n2, tmp_mat_n1n1); // Kij dB/dx_j
                            B_mat.right_multiply_transpose
                            (tmp_mat_n2n2, tmp_mat_n1n2); // B^T Kij dB/dx_j
                            Kmat.add(-JxW[qp], tmp_mat_n2n2);
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
            dnormal_steady_real.zero();
            
            for (unsigned int qp=0; qp<qpoint.size(); qp++)
            {
                // first update the variables at the current quadrature point
                this->update_solution_at_quadrature_point
                (vars, qp, c, false,
                 ref_sol, conservative_sol,
                 p_sol, B_mat, dB_mat);
                
                // initialize flux to interpolated sol for initialization
                // of perturbed vars
                B_mat.vector_mult(conservative_deltasol, elem_sol);
                
                delta_p_sol.zero();
                delta_p_sol.init(p_sol, conservative_deltasol);
                
                // now prepare the flux vector
                flux.zero();
                for (unsigned int i_dim=0; i_dim<dim; i_dim++)
                    flux(i_dim+1) += delta_p_sol.dp * face_normals[qp](i_dim);
                
                B_mat.vector_mult_transpose(tmp_vec1_n2, flux);
                Fvec.add(JxW[qp], tmp_vec1_n2);
                
                if ( request_jacobian && c.get_elem_solution_derivative() )
                {
                    this->calculate_advection_flux_jacobian_for_moving_solid_wall_boundary
                    (p_sol, 0., face_normals[qp], dnormal_steady_real, A_mat);
                    
                    B_mat.left_multiply(tmp_mat_n1n2, A_mat); // A B
                    B_mat.right_multiply_transpose(tmp_mat_n2n2,
                                                   tmp_mat_n1n2); // B^T A B
                    Kmat.add(JxW[qp], tmp_mat_n2n2);
                }
            }
        }
            break;

        
        case SLIP_WALL: // inviscid boundary condition without any diffusive component
            // conditions enforced are
            // vi ni = wi_dot (ni + dni) - ui dni       (slip wall)
            // tau_ij nj = 0   (because velocity gradient at wall = 0)
            // qi ni = 0       (since heat flux occurs only on no-slip wall and far-field bc)
        {
            Number dui_ni_unsteady = 0., ui_ni_steady = 0.;
            
            for (unsigned int qp=0; qp<qpoint.size(); qp++)
            {
                // first update the variables at the current quadrature point
                this->update_solution_at_quadrature_point
                (vars, qp, c, false,
                 ref_sol, conservative_sol,
                 p_sol, B_mat, dB_mat);
                
                // initialize flux to interpolated sol for initialization
                // of perturbed vars
                B_mat.vector_mult(conservative_deltasol, elem_sol);
                
                delta_p_sol.zero();
                delta_p_sol.init(p_sol, conservative_deltasol);
                
                dui_ni_unsteady = 0.; ui_ni_steady = 0.;
                p_sol.get_uvec(uvec);
                delta_p_sol.get_duvec(duvec);

                // calculate the surface velocity perturbations
                perturbed_surface_motion->surface_velocity_frequency_domain
                (qpoint[qp], face_normals[qp],
                 surface_unsteady_vel, dnormal_unsteady);

                // now check if the surface deformation is defined and
                // needs to be applied through transpiration boundary
                // condition
                if (surface_motion) // get the surface motion data
                {
                    surface_motion->surface_velocity_time_domain
                    (0., qpoint[qp], face_normals[qp],
                     surface_steady_vel, dnormal_steady);
                    
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
                    {
                        dui_ni_unsteady += // delta_wi_dot * ni
                        surface_unsteady_vel(i_dim) * face_normals[qp](i_dim);
                        ui_ni_steady    += // wi_dot * ni
                        surface_steady_vel(i_dim) * face_normals[qp](i_dim);
                    }
                }
                
                
                // now add the contribution from unsteady normal perturbation
                dui_ni_unsteady -= uvec.dot(dnormal_unsteady); // ui delta_ni

                
                // now prepare the flux vector
                flux.zero();
                flux.add(dui_ni_unsteady, conservative_sol);
                flux.add(ui_ni_steady, conservative_deltasol);
                flux(n1-1) += dui_ni_unsteady*p_sol.p + ui_ni_steady*delta_p_sol.dp;
                for (unsigned int i_dim=0; i_dim<dim; i_dim++)
                    flux(i_dim+1) += delta_p_sol.dp * face_normals[qp](i_dim);

                B_mat.vector_mult_transpose(tmp_vec1_n2, flux);
                Fvec.add(JxW[qp], tmp_vec1_n2);
                
                if ( request_jacobian && c.get_elem_solution_derivative() )
                {
                    for (unsigned int i=0; i<dim; i++)
                        dnormal_steady_real(i) = std::real(dnormal_steady(i));

                    this->calculate_advection_flux_jacobian_for_moving_solid_wall_boundary
                    (p_sol, std::real(ui_ni_steady), face_normals[qp],
                     dnormal_steady_real, A_mat);
                    
                    B_mat.left_multiply(tmp_mat_n1n2, A_mat); // A B
                    B_mat.right_multiply_transpose(tmp_mat_n2n2,
                                                   tmp_mat_n1n2); // B^T A B
                    Kmat.add(JxW[qp], tmp_mat_n2n2);
                }
            }
        }
            break;
            
            
        case FAR_FIELD:
            // conditions enforced are:
            // -- f_adv_i ni =  f_adv = f_adv(+) + f_adv(-)     (flux vector splitting for advection)
            // -- f_diff_i ni  = f_diff                         (evaluation of diffusion flux based on domain solution)
        {

            U_vec_interpolated.zero(); // that is, the perturbation in this solution is zero
            
            for (unsigned int qp=0; qp<qpoint.size(); qp++)
            {
                // first update the variables at the current quadrature point
                this->update_solution_at_quadrature_point
                (vars, qp, c, false,
                 ref_sol, conservative_sol, p_sol, B_mat, dB_mat);
                
                this->calculate_advection_left_eigenvector_and_inverse_for_normal
                (p_sol, face_normals[qp], eig_val, l_eig_vec, l_eig_vec_inv_tr);
                
                // for all eigenalues that are less than 0, the characteristics are coming into the domain, hence,
                // evaluate them using the given solution.
                
                // perturbation in all incoming characteristics is zero. Hence, that boundary condition can be zeroed out.
                
                // now calculate the flux for eigenvalues greater than 0, the characteristics go out of the domain, so that
                // the flux is evaluated using the local solution
                tmp_mat_n1n1 = l_eig_vec_inv_tr;
                for (unsigned int j=0; j<n1; j++)
                    if (eig_val(j, j) > 0)
                        tmp_mat_n1n1.scale_column(j, eig_val(j, j)); // L^-T [omaga]_{+}
                    else
                        tmp_mat_n1n1.scale_column(j, 0.0);
                
                tmp_mat_n1n1.right_multiply_transpose(l_eig_vec); // A_{+} = L^-T [omaga]_{+} L^T
                
                B_mat.vector_mult(elem_interpolated_sol, c.get_elem_solution()); // B dU
                tmp_mat_n1n1.vector_mult(flux, elem_interpolated_sol); // f_{+} = A_{+} B dU
                
                B_mat.vector_mult_transpose(tmp_vec1_n2, flux); // B^T f_{+}   (this is flux going out of the solution domain)
                Fvec.add(JxW[qp], tmp_vec1_n2);
                
                
                if (_if_viscous) // evaluate the viscous flux using the domain solution
                {
                    // stress tensor and temperature gradient
                    this->calculate_conservative_variable_jacobian
                    (p_sol, dcons_dprim, dprim_dcons);
                    this->calculate_diffusion_tensors
                    (ref_sol, dB_mat, dprim_dcons,
                     p_sol, stress_tensor, temp_grad);
                    
                    // contribution from diffusion flux
                    flux.zero();
                    for (unsigned int i_dim=0; i_dim<dim; i_dim++)
                    {
                        this->calculate_diffusion_flux(i_dim, p_sol, stress_tensor,
                                                       temp_grad, tmp_vec1);
                        flux.add(face_normals[qp](i_dim), tmp_vec1); // fi ni
                    }
                    B_mat.vector_mult_transpose(tmp_vec1_n2, flux);
                    Fvec.add(-JxW[qp], tmp_vec1_n2);
                }

                if ( request_jacobian && c.get_elem_solution_derivative() )
                {
                    // terms with negative eigenvalues do not contribute to the Jacobian
                    
                    // now calculate the Jacobian for eigenvalues greater than 0, the characteristics
                    // go out of the domain, so that the flux is evaluated using the local solution
                    tmp_mat_n1n1 = l_eig_vec_inv_tr;
                    for (unsigned int j=0; j<n1; j++)
                        if (eig_val(j, j) > 0)
                            tmp_mat_n1n1.scale_column(j, eig_val(j, j)); // L^-T [omaga]_{+}
                        else
                            tmp_mat_n1n1.scale_column(j, 0.0);
                    tmp_mat_n1n1.right_multiply_transpose(l_eig_vec); // A_{+} = L^-T [omaga]_{+} L^T
                    B_mat.left_multiply(tmp_mat_n1n2, tmp_mat_n1n1);
                    B_mat.right_multiply_transpose(tmp_mat_n2n2, tmp_mat_n1n2); // B^T A_{+} B   (this is flux going out of the solution domain)
                    
                    Kmat.add(JxW[qp], tmp_mat_n2n2);

                    if (_if_viscous)
                    {
                        // contribution from diffusion flux
                        for (unsigned int i_dim=0; i_dim<dim; i_dim++)
                            for (unsigned int deriv_dim=0; deriv_dim<dim; deriv_dim++)
                            {
                                this->calculate_diffusion_flux_jacobian
                                (i_dim, deriv_dim, p_sol, tmp_mat_n1n1); // Kij
                                dB_mat[deriv_dim].left_multiply
                                (tmp_mat_n1n2, tmp_mat_n1n1); // Kij dB/dx_j
                                B_mat.right_multiply_transpose
                                (tmp_mat_n2n2, tmp_mat_n1n2); // B^T Kij dB/dx_j
                                Kmat.add(-JxW[qp], tmp_mat_n2n2);
                            }
                    }
                }
            }
        }
    }

    
//    std::cout << "inside side constraint " << std::endl;
//    std::cout << "elem solution" << std::endl; c.elem_solution.print(std::cout);
//    std::cout << mechanical_bc_type << std::endl;
//    std::cout << "bc vec: " << std::endl; Fvec.print(std::cout);
//    if (request_jacobian && c.elem_solution_derivative)
//        Kmat.print(std::cout);
    
    return request_jacobian;
}






Real get_complex_var_val(const std::string& var_name, const SmallPerturbationPrimitiveSolution<Number>& delta_p_sol, Real q0)
{
    if (var_name == "du_re")
        return std::real(delta_p_sol.du1);
    else if (var_name == "dv_re")
        return std::real(delta_p_sol.du2);
    else if (var_name == "dw_re")
        return std::real(delta_p_sol.du3);
    else if (var_name == "dT_re")
        return std::real(delta_p_sol.dT);
    else if (var_name == "ds_re")
        return std::real(delta_p_sol.dentropy);
    else if (var_name == "dp_re")
        return std::real(delta_p_sol.dp);
    else if (var_name == "dcp_re")
        return std::real(delta_p_sol.c_pressure(q0));
    else if (var_name == "da_re")
        return std::real(delta_p_sol.da);
    else if (var_name == "dM_re")
        return std::real(delta_p_sol.dmach);
    else if (var_name == "du_im")
        return std::imag(delta_p_sol.du1);
    else if (var_name == "dv_im")
        return std::imag(delta_p_sol.du2);
    else if (var_name == "dw_im")
        return std::imag(delta_p_sol.du3);
    else if (var_name == "dT_im")
        return std::imag(delta_p_sol.dT);
    else if (var_name == "ds_im")
        return std::imag(delta_p_sol.dentropy);
    else if (var_name == "dp_im")
        return std::imag(delta_p_sol.dp);
    else if (var_name == "dcp_im")
        return std::imag(delta_p_sol.c_pressure(q0));
    else if (var_name == "da_im")
        return std::imag(delta_p_sol.da);
    else if (var_name == "dM_im")
        return std::imag(delta_p_sol.dmach);
    else
        libmesh_assert(false);
}



class FrequencyDomainPrimitiveFEMFunction : public FEMFunctionBase<Number>
{
public:
    // Constructor
    FrequencyDomainPrimitiveFEMFunction(MeshFunction& fluid_func,
                                        MeshFunction& sd_fluid_func,
                                        std::vector<std::string>& vars,
                                        Real cp, Real cv, Real q0):
    FEMFunctionBase<Number>(),
    _fluid_function(fluid_func),
    _sd_fluid_function(sd_fluid_func),
    _vars(vars), _cp(cp), _cv(cv), _q0(q0)
    {}
    
    // Destructor
    virtual ~FrequencyDomainPrimitiveFEMFunction () {}
    
    virtual AutoPtr<FEMFunctionBase<Number> > clone () const
    {return AutoPtr<FEMFunctionBase<Number> >( new FrequencyDomainPrimitiveFEMFunction
                                              (_fluid_function,
                                               _sd_fluid_function,
                                               _vars, _cp, _cv,
                                               _q0) ); }
    
    virtual void operator() (const FEMContext& c, const Point& p,
                             const Real t, DenseVector<Number>& val)
    {
        DenseVector<Number> fluid_sol_complex, sd_fluid_sol;
        DenseVector<Real> fluid_sol;
        
        // get the solution values at the point
        _fluid_function(p, t, fluid_sol_complex);
        _sd_fluid_function(p, t, sd_fluid_sol);
        
        // since the system returns complex, convert it into real numbers
        for (unsigned int i_var=0; i_var<fluid_sol_complex.size(); i_var++)
            fluid_sol(i_var) = std::real(fluid_sol_complex(i_var));

        PrimitiveSolution p_sol;
        SmallPerturbationPrimitiveSolution<Number> delta_p_sol;

        // now initialize the primitive variable contexts
        p_sol.init(c.get_dim(), fluid_sol, _cp, _cv, false);
        delta_p_sol.init(p_sol, sd_fluid_sol);
        
        for (unsigned int i=0; i<_vars.size(); i++)
            val(i) = get_complex_var_val(_vars[i], delta_p_sol, _q0);
    }
    
    
    virtual Number component(const FEMContext& c, unsigned int i_comp,
                             const Point& p, Real t=0.)
    {
        DenseVector<Number> fluid_sol_complex, sd_fluid_sol;
        DenseVector<Real> fluid_sol;
        
        // get the solution values at the point
        _fluid_function(p, t, fluid_sol_complex);
        _sd_fluid_function(p, t, sd_fluid_sol);
        
        // since the system returns complex, convert it into real numbers
        fluid_sol.resize(fluid_sol_complex.size());
        for (unsigned int i_var=0; i_var<fluid_sol_complex.size(); i_var++)
            fluid_sol(i_var) = std::real(fluid_sol_complex(i_var));
        
        PrimitiveSolution p_sol;
        SmallPerturbationPrimitiveSolution<Number> delta_p_sol;
        
        // now initialize the primitive variable contexts
        p_sol.init(c.get_dim(), fluid_sol, _cp, _cv, false);
        delta_p_sol.init(p_sol, sd_fluid_sol);
        
        return get_complex_var_val(_vars[i_comp], delta_p_sol, _q0);
    }
    
    
    virtual Number operator() (const FEMContext&, const Point& p,
                               const Real time = 0.)
    {libmesh_error();}
    
private:
    
    MeshFunction& _fluid_function;
    MeshFunction& _sd_fluid_function;
    std::vector<std::string>& _vars;
    Real _cp, _cv, _q0;
};





void FrequencyDomainFluidPostProcessSystem::init_data()
{
    const unsigned int dim = this->get_mesh().mesh_dimension();
    
    GetPot infile("euler.in");
        
    unsigned int o = infile("fe_order", 1);
    std::string fe_family = infile("fe_family", std::string("LAGRANGE"));
    FEFamily fefamily = Utility::string_to_enum<FEFamily>(fe_family);
    
    
    const Order order = static_cast<Order>(fmin(2, o));
    
    u = this->add_variable("du_re", order, fefamily);
    if (dim > 1)
        v = this->add_variable("dv_re", order, fefamily);
    if (dim > 2)
        w = this->add_variable("dw_re", order, fefamily);
    T = this->add_variable("dT_re", order, fefamily);
    s = this->add_variable("ds_re", order, fefamily);
    p = this->add_variable("dp_re", order, fefamily);
    cp = this->add_variable("dcp_re", order, fefamily);
    a = this->add_variable("da_re", order, fefamily);
    M = this->add_variable("dM_re", order, fefamily);
    
    u_im = this->add_variable("du_im", order, fefamily);
    if (dim > 1)
        v_im = this->add_variable("dv_im", order, fefamily);
    if (dim > 2)
        w_im = this->add_variable("dw_im", order, fefamily);
    T_im = this->add_variable("dT_im", order, fefamily);
    s_im = this->add_variable("ds_im", order, fefamily);
    p_im = this->add_variable("dp_im", order, fefamily);
    cp_im = this->add_variable("dcp_im", order, fefamily);
    a_im = this->add_variable("da_im", order, fefamily);
    M_im = this->add_variable("dM_im", order, fefamily);
    
    System::init_data();
}





void FrequencyDomainFluidPostProcessSystem::postprocess()
{

    // initialize the mesh function for the fluid solution,
    // this will be used for calculation of the element solution
    const System& fluid = this->get_equation_systems().get_system<System>
    ("EulerSystem");
    const FrequencyDomainLinearizedFluidSystem& sd_fluid =
    this->get_equation_systems().get_system<FrequencyDomainLinearizedFluidSystem>
    ("FrequencyDomainLinearizedFluidSystem");
    
    std::vector<std::string> post_process_var_names(this->n_vars());
    for (unsigned int i=0; i<this->n_vars(); i++)
        post_process_var_names[i] = this->variable_name(i);
    
    
    AutoPtr<NumericVector<Number> >
    fluid_sol =
    NumericVector<Number>::build(this->get_equation_systems().comm()),  // for the nonlinear system
    sd_fluid_sol =
    NumericVector<Number>::build(this->get_equation_systems().comm());  // for the small-disturbance system

    // ask systems to localize their solutions
    fluid_sol->init(fluid.solution->size(), true, SERIAL);
    fluid.solution->localize(*fluid_sol,
                             fluid.get_dof_map().get_send_list());
    
    // now the small-disturbance solution
    sd_fluid_sol->init(sd_fluid.solution->size(), true, SERIAL);
    sd_fluid.solution->localize(*sd_fluid_sol,
                                sd_fluid.get_dof_map().get_send_list());

    // create the mesh functions to interpolate solutions
    AutoPtr<MeshFunction>
    fluid_mesh_function
    (new MeshFunction(this->get_equation_systems(), *fluid_sol,
                      fluid.get_dof_map(), sd_fluid.vars)),
    sd_fluid_mesh_function
    (new MeshFunction(this->get_equation_systems(), *sd_fluid_sol,
                      sd_fluid.get_dof_map(), sd_fluid.vars));

    fluid_mesh_function->init();
    sd_fluid_mesh_function->init();

    
    AutoPtr<FEMFunctionBase<Number> > post_process_function
    (new FrequencyDomainPrimitiveFEMFunction
     (*fluid_mesh_function, *sd_fluid_mesh_function,
      post_process_var_names,
      flight_condition->gas_property.cp,
      flight_condition->gas_property.cv,
      flight_condition->q0()));
    
    this->project_solution(post_process_function.get());
    
    this->update();
}





#endif // LIBMESH_USE_COMPLEX_NUMBERS