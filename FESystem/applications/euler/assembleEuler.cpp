//
//  assembleEuler.cpp
//  FESystem
//
//  Created by Manav Bhatia on 2/21/13.
//
//

// C++ includes
#include <iomanip>

// FESystem includes
#include "euler/assembleEuler.h"

#ifndef LIBMESH_USE_COMPLEX_NUMBERS

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
    this->verify_analytic_jacobians = infile("verify_analytic_jacobians", 0.);
    this->print_jacobians = infile("print_jacobians", false);
    this->print_element_jacobians = infile("print_element_jacobians", false);
    
    // initialize the Dirichlet boundary conditions
    std::multimap<unsigned, FluidBoundaryConditionType>::const_iterator
    bc_it = _boundary_condition.begin(), bc_end = _boundary_condition.end();
    
    std::set<boundary_id_type> no_slip_boundary, isothermal_boundary;
    
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
        ZeroFunction<Real> zero_function;
        std::vector<unsigned int> rho_ui_vars(dim);
        for (unsigned int i_dim=0; i_dim<dim; i_dim++)
            rho_ui_vars[i_dim] = vars[i_dim+1];
        
        this->get_dof_map().add_dirichlet_boundary(DirichletBoundary(no_slip_boundary, rho_ui_vars,
                                                                     &zero_function));
    }

    
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
        if (_if_viscous)
            elem_fe[i]->get_d2phi();
    }
    
    std::vector<FEBase*> elem_side_fe(dim+2);
    
    for (unsigned int i=0; i<dim+2; i++)
    {
        c.get_element_fe( vars[i], elem_side_fe[i]);
        elem_side_fe[i]->get_JxW();
        elem_side_fe[i]->get_phi();
        elem_fe[i]->get_xyz();
        if (_if_viscous)
            elem_fe[i]->get_dphi();
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

    FEMOperatorMatrix B_mat;
    std::vector<FEMOperatorMatrix> dB_mat(dim);
    std::vector<DenseMatrix<Real> >  Ai_advection(dim);
    DenseMatrix<Real> LS_mat, LS_sens, Ai_Bi_advection, A_entropy, A_inv_entropy, tmp_mat_n1n1,
    tmp_mat_n1n2, tmp_mat2_n2n2, tmp_mat3, A_sens, stress_tensor, dprim_dcons, dcons_dprim;
    
    DenseVector<Real> flux, tmp_vec1_n1, tmp_vec2_n1, tmp_vec3_n2, conservative_sol, delta_vals, temp_grad;
    
    LS_mat.resize(n1, n_dofs); LS_sens.resize(n_dofs, n_dofs); Ai_Bi_advection.resize(dim+2, n_dofs);
    A_inv_entropy.resize(dim+2, dim+2); A_entropy.resize(dim+2, dim+2); tmp_mat_n1n2.resize(dim+2, n_dofs); A_sens.resize(n1, n_dofs);
    stress_tensor.resize(dim, dim); dprim_dcons.resize(dim+2, dim+2); dcons_dprim.resize(dim+2, dim+2);
    tmp_mat2_n2n2.resize(n_dofs, n_dofs);
    
    flux.resize(n1); tmp_vec1_n1.resize(n1); tmp_vec2_n1.resize(n1); tmp_vec3_n2.resize(n_dofs); conservative_sol.resize(dim+2);
    delta_vals.resize(n_qpoints); temp_grad.resize(dim);
    
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

    System& delta_val_system = this->get_equation_systems().get_system<System>("DeltaValSystem");
    NumericVector<Real>& diff_val_vec = (*delta_val_system.solution.get());
    
    if (if_use_stored_dc_coeff)
    {
        diff_val = diff_val_vec.el(c.elem->dof_number(delta_val_system.number(), 0, 0));
        for (unsigned int qp=0; qp<n_qpoints; qp++)
            delta_vals(qp) = diff_val;
    }
    
    PrimitiveSolution primitive_sol;
    const std::vector<std::vector<Real> >& phi = elem_fe->get_phi(); // assuming that all variables have the same interpolation
    const unsigned int n_phi = phi.size();
    
    for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
        // first update the variables at the current quadrature point
        this->update_solution_at_quadrature_point(vars, qp, c, true, c.elem_solution,
                                                  conservative_sol, primitive_sol,
                                                  B_mat, dB_mat);

        for ( unsigned int i_dim=0; i_dim < dim; i_dim++)
            this->calculate_advection_flux_jacobian ( i_dim, primitive_sol, Ai_advection[i_dim] );
        
        this->calculate_entropy_variable_jacobian ( primitive_sol, A_entropy, A_inv_entropy );


        if (_if_viscous)
        {
            // stress tensor and temperature gradient
            this->calculate_conservative_variable_jacobian(primitive_sol, dcons_dprim, dprim_dcons);
            this->calculate_diffusion_tensors(c.elem_solution, dB_mat, dprim_dcons,
                                              primitive_sol, stress_tensor, temp_grad);
        }
        
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
            this->calculate_differential_operator_matrix(vars, qp, c, c.elem_solution, primitive_sol, B_mat,
                                                         dB_mat, Ai_advection, Ai_Bi_advection,
                                                         A_inv_entropy, flux_jacobian_sens, LS_mat, LS_sens, diff_val);

        if (if_use_stored_dc_coeff)
            diff_val = delta_vals(qp);
        else
            delta_vals(qp) = diff_val;
        
        for (unsigned int i_dim=0; i_dim<dim; i_dim++)
        {
            // Galerkin contribution from the advection flux terms
            this->calculate_advection_flux(i_dim, primitive_sol, flux); // F^adv_i
            dB_mat[i_dim].vector_mult_transpose(tmp_vec3_n2, flux); // dBw/dx_i F^adv_i
            Fvec.add(JxW[qp], tmp_vec3_n2);
            
            if (_if_viscous)
            {
                // Galerkin contribution from the diffusion flux terms
                this->calculate_diffusion_flux(i_dim, primitive_sol, stress_tensor, temp_grad, flux); // F^diff_i
                dB_mat[i_dim].vector_mult_transpose(tmp_vec3_n2, flux); // dBw/dx_i F^diff_i
                Fvec.add(-JxW[qp], tmp_vec3_n2);
            }
            
            // discontinuity capturing operator
            dB_mat[i_dim].vector_mult(flux, c.elem_solution);
            dB_mat[i_dim].vector_mult_transpose(tmp_vec3_n2, flux);
            Fvec.add(-JxW[qp]*diff_val, tmp_vec3_n2);
        }
        
        // Least square contribution from divergence of advection flux
        Ai_Bi_advection.vector_mult(flux, c.elem_solution); // d F^adv_i / dxi
        LS_mat.vector_mult_transpose(tmp_vec3_n2, flux); // LS^T tau F^adv_i
        Fvec.add(-JxW[qp], tmp_vec3_n2);
        
        // Least square contribution from divergence of diffusion flux
        // TODO: this requires a 2nd order differential of the flux
        
        
        if (request_jacobian && c.elem_solution_derivative)
        {
            A_sens.zero();

            for (unsigned int i_dim=0; i_dim<dim; i_dim++)
            {
                // Galerkin contribution from the advection flux terms
                B_mat.left_multiply(tmp_mat_n1n2, Ai_advection[i_dim]);
                dB_mat[i_dim].right_multiply_transpose(tmp_mat2_n2n2, tmp_mat_n1n2);// dBw/dx_i^T  dF^adv_i/ dU
                Kmat.add(JxW[qp], tmp_mat2_n2n2);
                
                // discontinuity capturing term
                dB_mat[i_dim].right_multiply_transpose(tmp_mat2_n2n2, dB_mat[i_dim]);
                Kmat.add(-JxW[qp]*diff_val, tmp_mat2_n2n2);
                
                if (_if_full_linearization)
                {
                    // sensitivity of Ai_Bi with respect to U:   [dAi/dUj.Bi.U  ...  dAi/dUn.Bi.U]
                    dB_mat[i_dim].vector_mult(tmp_vec2_n1, c.elem_solution);
                    for (unsigned int i_cvar=0; i_cvar<n1; i_cvar++)
                    {
                        flux_jacobian_sens[i_dim][i_cvar].vector_mult(tmp_vec1_n1, tmp_vec2_n1);
                        for (unsigned int i_phi=0; i_phi<n_phi; i_phi++)
                            A_sens.add_column((n_phi*i_cvar)+i_phi, phi[i_phi][qp], tmp_vec1_n1); // assuming that all variables have same n_phi
                    }
                }
                
                if (_if_viscous)
                {
                    for (unsigned int deriv_dim=0; deriv_dim<dim; deriv_dim++)
                    {
                        this->calculate_diffusion_flux_jacobian(i_dim, deriv_dim, primitive_sol, tmp_mat_n1n1); // Kij
                        dB_mat[deriv_dim].left_multiply(tmp_mat_n1n2, tmp_mat_n1n1); // Kij dB/dx_j
                        dB_mat[i_dim].right_multiply_transpose(tmp_mat2_n2n2, tmp_mat_n1n2); // dB/dx_i^T Kij dB/dx_j
                        Kmat.add(-JxW[qp], tmp_mat2_n2n2);
                    }

                }
            }
            
            // Least square contribution of flux gradient
            LS_mat.get_transpose(tmp_mat3);
            tmp_mat3.right_multiply(Ai_Bi_advection); // LS^T tau d^2F^adv_i / dx dU   (Ai constant)
            Kmat.add(-JxW[qp], tmp_mat3);

            if (_if_full_linearization)
            {
                LS_mat.get_transpose(tmp_mat3);
                tmp_mat3.right_multiply(A_sens); // LS^T tau d^2F^adv_i / dx dU  (Ai sensitivity)
                Kmat.add(-JxW[qp], tmp_mat3);
                
                // contribution sensitivity of the LS.tau matrix
                Kmat.add(-JxW[qp], LS_sens);
            }
        }
    } // end of the quadrature point qp-loop

    if (!if_use_stored_dc_coeff)
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

    FEBase * side_fe;
    c.get_side_fe(vars[0], side_fe); // assuming all variables have the same FE
    
    // Element Jacobian * quadrature weights for interior integration
    const std::vector<Real> &JxW = side_fe->get_JxW();
    
    // Physical location of the quadrature points
    const std::vector<Point>& qpoint = side_fe->get_xyz();

    // boundary normals
    const std::vector<Point>& face_normals = side_fe->get_normals();


    Point vel; vel.zero(); // zero surface velocity
    
    
    DenseMatrix<Number>& Kmat = c.get_elem_jacobian();
    DenseVector<Number>& Fvec = c.get_elem_residual();


    
    FEMOperatorMatrix B_mat;
    DenseVector<Real> tmp_vec1_n2, flux, U_vec_interpolated, tmp_vec2_n1, conservative_sol, temp_grad;
    DenseMatrix<Real>  eig_val, l_eig_vec, l_eig_vec_inv_tr, tmp_mat_n1n1, tmp_mat1_n1n2, tmp_mat2_n2n2,
    A_mat, dcons_dprim, dprim_dcons, stress_tensor;

    conservative_sol.resize(dim+2); temp_grad.resize(dim);
    tmp_vec1_n2.resize(n_dofs); flux.resize(n1); tmp_vec2_n1.resize(n1); U_vec_interpolated.resize(n1);
    eig_val.resize(n1, n1); l_eig_vec.resize(n1, n1); l_eig_vec_inv_tr.resize(n1, n1);
    tmp_mat1_n1n2.resize(n1, n_dofs); tmp_mat2_n2n2.resize(n_dofs, n_dofs);
    A_mat.resize(dim+2, dim+2); tmp_mat_n1n1.resize(n1, n1);
    dcons_dprim.resize(n1, n1); dprim_dcons.resize(n1, n1); stress_tensor.resize(dim, dim);

    std::vector<FEMOperatorMatrix> dB_mat(dim);
    std::vector<DenseMatrix<Real> > Ai_advection(dim);
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

            Real xini = 0.;
            
            for (unsigned int qp=0; qp<qpoint.size(); qp++)
            {
                xini = face_normals[qp] * vel;
                
                // first update the variables at the current quadrature point
                this->update_solution_at_quadrature_point(vars, qp, c, false, c.elem_solution, conservative_sol, p_sol,
                                                          B_mat, dB_mat);

                // stress tensor and temperature gradient
                this->calculate_conservative_variable_jacobian(p_sol, dcons_dprim, dprim_dcons);
                this->calculate_diffusion_tensors(c.elem_solution, dB_mat, dprim_dcons,
                                                  p_sol, stress_tensor, temp_grad);
                
                if (thermal_bc_type == ADIABATIC)
                    temp_grad.zero();

                flux.zero();
                
                // since vi=0, vi ni = 0, so the advection flux gets evaluated using the slip wall condition
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
                
                // contribution from advection flux
                B_mat.vector_mult_transpose(tmp_vec1_n2, flux);
                Fvec.add(-JxW[qp], tmp_vec1_n2);
                
                // contribution from diffusion flux
                flux.zero();
                for (unsigned int i_dim=0; i_dim<dim; i_dim++)
                {
                    this->calculate_diffusion_flux(i_dim, p_sol, stress_tensor,
                                                   temp_grad, tmp_vec2_n1);
                    flux.add(face_normals[qp](i_dim), tmp_vec2_n1); // fi ni
                }
                B_mat.vector_mult_transpose(tmp_vec1_n2, flux);
                Fvec.add(JxW[qp], tmp_vec1_n2);
                
                
                if ( request_jacobian && c.get_elem_solution_derivative() )
                {
                    this->calculate_advection_flux_jacobian_for_moving_solid_wall_boundary(p_sol, xini, face_normals[qp], A_mat);
                    
                    // contribution from advection flux
                    B_mat.left_multiply(tmp_mat1_n1n2, A_mat);
                    B_mat.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat1_n1n2);
                    Kmat.add(-JxW[qp], tmp_mat2_n2n2);

                    // contribution from diffusion flux
                    for (unsigned int i_dim=0; i_dim<dim; i_dim++)
                        for (unsigned int deriv_dim=0; deriv_dim<dim; deriv_dim++)
                        {
                            this->calculate_diffusion_flux_jacobian(i_dim, deriv_dim, p_sol, tmp_mat_n1n1); // Kij
                            dB_mat[deriv_dim].left_multiply(tmp_mat1_n1n2, tmp_mat_n1n1); // Kij dB/dx_j
                            B_mat.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat1_n1n2); // B^T Kij dB/dx_j
                            Kmat.add(JxW[qp], tmp_mat2_n2n2);
                        }
                }
            }
        }
            break;
            
        case SLIP_WALL: // inviscid boundary condition without any diffusive component
            // conditions enforced are
            // vi ni = 0       (slip wall)
            // tau_ij nj = 0   (because velocity gradient at wall = 0)
            // qi ni = 0       (since heat flux occurs only on no-slip wall and far-field bc)
        {
            Real xini = 0.;
            
            for (unsigned int qp=0; qp<qpoint.size(); qp++)
            {
                xini = face_normals[qp] * vel;
                
                // first update the variables at the current quadrature point
                this->update_solution_at_quadrature_point(vars, qp, c, false, c.elem_solution, conservative_sol, p_sol,
                                                          B_mat, dB_mat);
                
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
                
                B_mat.vector_mult_transpose(tmp_vec1_n2, flux);
                Fvec.add(-JxW[qp], tmp_vec1_n2);
                
                if ( request_jacobian && c.get_elem_solution_derivative() )
                {
                    this->calculate_advection_flux_jacobian_for_moving_solid_wall_boundary(p_sol, xini, face_normals[qp], A_mat);
                    
                    B_mat.left_multiply(tmp_mat1_n1n2, A_mat);
                    B_mat.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat1_n1n2);
                    Kmat.add(-JxW[qp], tmp_mat2_n2n2);
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
                this->update_solution_at_quadrature_point(vars, qp, c, false, c.elem_solution, conservative_sol, p_sol,
                                                          B_mat, dB_mat);
                
                this->calculate_advection_left_eigenvector_and_inverse_for_normal(p_sol, face_normals[qp], eig_val,
                                                                                  l_eig_vec, l_eig_vec_inv_tr);
                
                // for all eigenalues that are less than 0, the characteristics are coming into the domain, hence,
                // evaluate them using the given solution.
                tmp_mat_n1n1 = l_eig_vec_inv_tr;
                for (unsigned int j=0; j<n1; j++)
                    if (eig_val(j, j) < 0)
                        tmp_mat_n1n1.scale_column(j, eig_val(j, j)); // L^-T [omaga]_{-}
                    else
                        tmp_mat_n1n1.scale_column(j, 0.0);
                
                tmp_mat_n1n1.right_multiply_transpose(l_eig_vec); // A_{-} = L^-T [omaga]_{-} L^T
                this->get_infinity_vars( U_vec_interpolated );
                tmp_mat_n1n1.vector_mult(flux, U_vec_interpolated);  // f_{-} = A_{-} B U
                
                B_mat.vector_mult_transpose(tmp_vec1_n2, flux); // B^T f_{-}   (this is flux coming into the solution domain)
                Fvec.add(-JxW[qp], tmp_vec1_n2);
                
                // now calculate the flux for eigenvalues greater than 0,
                // the characteristics go out of the domain, so that
                // the flux is evaluated using the local solution
                tmp_mat_n1n1 = l_eig_vec_inv_tr;
                for (unsigned int j=0; j<n1; j++)
                    if (eig_val(j, j) > 0)
                        tmp_mat_n1n1.scale_column(j, eig_val(j, j)); // L^-T [omaga]_{+}
                    else
                        tmp_mat_n1n1.scale_column(j, 0.0);
                
                tmp_mat_n1n1.right_multiply_transpose(l_eig_vec); // A_{+} = L^-T [omaga]_{+} L^T
                tmp_mat_n1n1.vector_mult(flux, conservative_sol); // f_{+} = A_{+} B U
                
                B_mat.vector_mult_transpose(tmp_vec1_n2, flux); // B^T f_{+}   (this is flux going out of the solution domain)
                Fvec.add(-JxW[qp], tmp_vec1_n2);
                
                
                if (_if_viscous) // evaluate the viscous flux using the domain solution
                {
                    // stress tensor and temperature gradient
                    this->calculate_conservative_variable_jacobian(p_sol, dcons_dprim, dprim_dcons);
                    this->calculate_diffusion_tensors(c.elem_solution, dB_mat, dprim_dcons,
                                                      p_sol, stress_tensor, temp_grad);
                    
                    // contribution from diffusion flux
                    flux.zero();
                    for (unsigned int i_dim=0; i_dim<dim; i_dim++)
                    {
                        this->calculate_diffusion_flux(i_dim, p_sol, stress_tensor,
                                                       temp_grad, tmp_vec2_n1);
                        flux.add(face_normals[qp](i_dim), tmp_vec2_n1); // fi ni
                    }
                    B_mat.vector_mult_transpose(tmp_vec1_n2, flux);
                    Fvec.add(JxW[qp], tmp_vec1_n2);
                }
                
            
                if ( request_jacobian && c.get_elem_solution_derivative() )
                {
                    // terms with negative eigenvalues do not contribute to the Jacobian
                    
                    // now calculate the Jacobian for eigenvalues greater than 0,
                    // the characteristics go out of the domain, so that
                    // the flux is evaluated using the local solution
                    tmp_mat_n1n1 = l_eig_vec_inv_tr;
                    for (unsigned int j=0; j<n1; j++)
                        if (eig_val(j, j) > 0)
                            tmp_mat_n1n1.scale_column(j, eig_val(j, j)); // L^-T [omaga]_{+}
                        else
                            tmp_mat_n1n1.scale_column(j, 0.0);
                    tmp_mat_n1n1.right_multiply_transpose(l_eig_vec); // A_{+} = L^-T [omaga]_{+} L^T
                    B_mat.left_multiply(tmp_mat1_n1n2, tmp_mat_n1n1);
                    B_mat.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat1_n1n2); // B^T A_{+} B   (this is flux going out of the solution domain)
                    
                    Kmat.add(-JxW[qp], tmp_mat2_n2n2);
                    
                    if (_if_viscous)
                    {
                        // contribution from diffusion flux
                        for (unsigned int i_dim=0; i_dim<dim; i_dim++)
                            for (unsigned int deriv_dim=0; deriv_dim<dim; deriv_dim++)
                            {
                                this->calculate_diffusion_flux_jacobian(i_dim, deriv_dim, p_sol, tmp_mat_n1n1); // Kij
                                dB_mat[deriv_dim].left_multiply(tmp_mat1_n1n2, tmp_mat_n1n1); // Kij dB/dx_j
                                B_mat.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat1_n1n2); // B^T Kij dB/dx_j
                                Kmat.add(JxW[qp], tmp_mat2_n2n2);
                            }
                    }
                }
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
    
    FEMOperatorMatrix B_mat;
    std::vector<FEMOperatorMatrix> dB_mat(dim);
    std::vector<DenseMatrix<Real> > Ai_advection(dim);
    DenseMatrix<Real> LS_mat, LS_sens, Ai_Bi_advection, A_entropy, A_inv_entropy, tmp_mat_n1n2,
    tmp_mat2_n2n2, tmp_mat3;
    DenseVector<Real> flux, tmp_vec1_n1, tmp_vec3_n2, conservative_sol;
    LS_mat.resize(n1, n_dofs); LS_sens.resize(n_dofs, n_dofs); tmp_mat2_n2n2.resize(n_dofs, n_dofs);
    Ai_Bi_advection.resize(dim+2, n_dofs); A_inv_entropy.resize(dim+2, dim+2);
    A_entropy.resize(dim+2, dim+2); tmp_mat_n1n2.resize(dim+2, n_dofs);
    flux.resize(n1); tmp_vec1_n1.resize(n1); tmp_vec3_n2.resize(n_dofs); conservative_sol.resize(dim+2);
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
        
    for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
        // first update the variables at the current quadrature point
        this->update_solution_at_quadrature_point(vars, qp, c, true, c.elem_fixed_solution,
                                                  conservative_sol, primitive_sol,
                                                  B_mat, dB_mat);

        for ( unsigned int i_dim=0; i_dim < dim; i_dim++)
            this->calculate_advection_flux_jacobian ( i_dim, primitive_sol, Ai_advection[i_dim] );
        
        this->calculate_entropy_variable_jacobian ( primitive_sol, A_entropy, A_inv_entropy );
        
        Ai_Bi_advection.zero();
        
        for (unsigned int i_dim=0; i_dim<dim; i_dim++)
        {
            dB_mat[i_dim].left_multiply(tmp_mat_n1n2, Ai_advection[i_dim]);
            Ai_Bi_advection.add(1.0, tmp_mat_n1n2);
        }
        
        if (_if_update_stabilization_per_quadrature_point || (qp == 0))
            this->calculate_differential_operator_matrix(vars, qp, c, c.elem_fixed_solution, primitive_sol, B_mat, dB_mat,
                                                         Ai_advection, Ai_Bi_advection, A_inv_entropy,
                                                         flux_jacobian_sens, LS_mat, LS_sens, diff_val);
        
        // Galerkin contribution to velocity
        B_mat.vector_mult( tmp_vec1_n1, c.elem_solution );
        B_mat.vector_mult_transpose(tmp_vec3_n2, tmp_vec1_n1);
        Fvec.add(JxW[qp], tmp_vec3_n2);
        
        // LS contribution to velocity
        B_mat.vector_mult( tmp_vec1_n1, c.elem_solution );
        LS_mat.vector_mult_transpose(tmp_vec3_n2, tmp_vec1_n1);
        Fvec.add(JxW[qp], tmp_vec3_n2);
        
        
        if (request_jacobian && c.get_elem_solution_derivative())
        {
            // contribution from unsteady term
            // Galerkin contribution of velocity
            B_mat.right_multiply_transpose(tmp_mat2_n2n2, B_mat);
            Kmat.add(JxW[qp], tmp_mat2_n2n2);
            
            // LS contribution of velocity
            LS_mat.get_transpose(tmp_mat3);
            B_mat.left_multiply(tmp_mat2_n2n2, tmp_mat3); // LS^T tau Bmat
            Kmat.add(JxW[qp], tmp_mat2_n2n2);
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
    for (unsigned int i=0; i<this->qoi.size(); i++)
        if (i == 0)
            libMesh::out << "Lift: " << std::setw(25) << std::setprecision(14) << qoi[i] << std::endl;
        else if (i == 1)
            libMesh::out << "Drag: " << std::setw(25) << std::setprecision(14) << qoi[i] << std::endl;
        else if (i == 2)
            libMesh::out << "Total V: " << std::setw(25) << std::setprecision(14) << qoi[i] << std::endl;
        else if (i == 3)
            libMesh::out << "Entropy Error: " << std::setw(25) << std::setprecision(14) << sqrt(qoi[i]/qoi[i-1]) << std::endl;
}





void EulerSystem::evaluate_recalculate_dc_flag()
{
    Real norm = this->calculate_norm(*(this->solution), this->vars[0], L2),
    relative_change0 = fabs(_rho_norm_curr - _rho_norm_old)/_rho_norm_curr,
    relative_change1 = fabs(norm - _rho_norm_curr)/norm;
    
    
    libMesh::out << "Rho L2-norm delta: " << relative_change0 << " , " << relative_change1 << std::endl;

    // if the relative change in the density norm is less than the threshold, then hold
    // the dc_coeffs constant
    if ((relative_change0 < this->dc_recalculate_tolerance) &&
        (relative_change1 < this->dc_recalculate_tolerance))
    {
        this->if_use_stored_dc_coeff = true;
        libMesh::out << "Using stored dc coeff" << std::endl;
    }
    else
    {
        if_use_stored_dc_coeff = false;
        libMesh::out << "Recalculating dc coeff" << std::endl;
    }
    
    _rho_norm_old = _rho_norm_curr;
    _rho_norm_curr = norm;
}





Real get_var_val(const std::string& var_name, const PrimitiveSolution& p_sol, Real p0, Real q0)
{
    if (var_name == "u")
        return p_sol.u1;
    else if (var_name == "v")
        return p_sol.u2;
    else if (var_name == "w")
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



class PrimitiveFEMFunction : public FEMFunctionBase<Number>
{
public:
    // Constructor
    PrimitiveFEMFunction(MeshFunction& fluid_func, std::vector<std::string>& vars,
                        Real cp, Real cv, Real p0, Real q0):
    FEMFunctionBase<Number>(),
    _fluid_function(fluid_func),
    _vars(vars), _cp(cp), _cv(cv), _p0(p0), _q0(q0)
    {}
    
    // Destructor
    virtual ~PrimitiveFEMFunction () {}
    
    virtual AutoPtr<FEMFunctionBase<Number> > clone () const
    {return AutoPtr<FEMFunctionBase<Number> >( new PrimitiveFEMFunction
                                              (_fluid_function, _vars, _cp, _cv,
                                               _p0, _q0) ); }
    
    virtual void operator() (const FEMContext& c, const Point& p,
                             const Real t, DenseVector<Number>& val)
    {
        DenseVector<Number> fluid_sol;
        _fluid_function(p, t, fluid_sol);
        PrimitiveSolution p_sol;
        
        p_sol.init(c.dim, fluid_sol, _cp, _cv, false);
        
        for (unsigned int i=0; i<_vars.size(); i++)
            val(i) = get_var_val(_vars[i], p_sol, _p0, _q0);
    }
    
    
    virtual Number component(const FEMContext& c, unsigned int i_comp,
                             const Point& p, Real t=0.)
    {
        DenseVector<Number> fluid_sol;
        _fluid_function(p, t, fluid_sol);
        PrimitiveSolution p_sol;
        
        p_sol.init(c.dim, fluid_sol, _cp, _cv, false);
        
        return get_var_val(_vars[i_comp], p_sol, _p0, _q0);
    }
    
    
    virtual Number operator() (const FEMContext&, const Point& p,
                               const Real time = 0.)
    {libmesh_error();}
    
private:
    
    MeshFunction& _fluid_function;
    std::vector<std::string>& _vars;
    Real _cp, _cv, _p0, _q0;
};




void FluidPostProcessSystem::init_data()
{
    const unsigned int dim = this->get_mesh().mesh_dimension();
    
    GetPot infile("euler.in");
    
    unsigned int o = infile("fe_order", 1);
    std::string fe_family = infile("fe_family", std::string("LAGRANGE"));
    FEFamily fefamily = Utility::string_to_enum<FEFamily>(fe_family);
    
    
    const Order order = static_cast<Order>(fmin(2, o));
    
    u = this->add_variable("u", order, fefamily);
    if (dim > 1)
        v = this->add_variable("v", order, fefamily);
    if (dim > 2)
        w = this->add_variable("w", order, fefamily);
    T = this->add_variable("T", order, fefamily);
    s = this->add_variable("s", order, fefamily);
    p = this->add_variable("p", order, fefamily);
    cp = this->add_variable("cp", order, fefamily);
    a = this->add_variable("a", order, fefamily);
    M = this->add_variable("M", order, fefamily);
    
    System::init_data();
}





void FluidPostProcessSystem::postprocess()
{
    // get the solution vector from
    const EulerSystem& fluid = this->get_equation_systems().get_system<EulerSystem>("EulerSystem");
    
    std::vector<unsigned int> fluid_vars(fluid.n_vars());
    std::vector<std::string> post_process_var_names(this->n_vars());
    
    for (unsigned int i=0; i<fluid_vars.size(); i++)
        fluid_vars[i] = i;
    for (unsigned int i=0; i<this->n_vars(); i++)
        post_process_var_names[i] = this->variable_name(i);
    
    
    AutoPtr<NumericVector<Number> > soln = NumericVector<Number>::build(this->get_equation_systems().comm());
    std::vector<Number> global_soln;
    fluid.update_global_solution(global_soln);
    soln->init(fluid.solution->size(), true, SERIAL);
    (*soln) = global_soln;
    
    AutoPtr<MeshFunction> mesh_function
    (new MeshFunction(this->get_equation_systems(), *soln,
                      fluid.get_dof_map(), fluid_vars));
    mesh_function->init();
    
    AutoPtr<FEMFunctionBase<Number> > post_process_function
    (new PrimitiveFEMFunction(*mesh_function, post_process_var_names, fluid.cp, fluid.cv,
                              fluid.p_inf, fluid.q0_inf));
    
    this->project_solution(post_process_function.get());
    
    this->update();
}




#endif // LIBMESH_USE_COMPLEX_NUMBERS

