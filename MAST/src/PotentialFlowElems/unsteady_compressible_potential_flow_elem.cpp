//
//  unsteady_compressible_potential_flow_elem.cpp
//  RealSolver
//
//  Created by Manav Bhatia on 10/7/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

// MAST includes
#include "PotentialFlowElems/unsteady_compressible_potential_flow_elem.h"

// libMesh includes
#include "libmesh/string_to_enum.h"
#include "libmesh/fem_context.h"
#include "libmesh/quadrature.h"
#include "libmesh/dense_vector.h"


Real unsteady_compressible_potential_solution_value(const Point& p,
                                                    const Parameters& parameters,
                                                    const std::string& sys_name,
                                                    const std::string& var_name)
{
    libmesh_assert_equal_to (sys_name, "FluidSystem");
    
    // since we are initializing the solution, variable values at all points is same
    
    if (var_name == "rho")
        
    {
        return parameters.get<Real> ("rho_inf");
    }
    else if (var_name == "phi")
    {
        Point dphi;
        dphi(0) = parameters.get<Real>("u1_inf"),
        dphi(1) = parameters.get<Real>("u2_inf"),
        dphi(2) = parameters.get<Real>("u3_inf");
        // assuming phi=0 at (0,0,0), and using the gradients, calculate the
        // value of phi at point p
        
        Real phi=0.;
        for (unsigned int i=0; i<3; i++)
            phi += dphi(i) * p(i);
        
        return phi;
    }
    else
        libmesh_assert(false);
}



void init_compressible_potential_variables(EquationSystems& es,
                                           const std::string& system_name)
{
    // It is a good idea to make sure we are initializing
    // the proper system.
    libmesh_assert_equal_to (system_name, "UnsteadyCompressiblePotentialSystem");
    
    // Get a reference to the Convection-Diffusion system object.
    UnsteadyCompressiblePotentialFlow & system =
    es.get_system<UnsteadyCompressiblePotentialFlow>("UnsteadyCompressiblePotentialSystem");
    
    // Project initial conditions at time 0
    libmesh_assert_equal_to(system.time, 0.);
    
    system.project_solution(unsteady_compressible_potential_solution_value,
                            NULL, es.parameters);
}


void UnsteadyCompressiblePotentialFlow::init_data ()
{
    this->use_fixed_solution = true;
    
    dim = this->get_mesh().mesh_dimension();
    
    vars.resize(2);
    
    // initialize the fluid values
    PotentialFlowElemBase::init_data();
    
    // Check the input file for Reynolds number, application type,
    // approximation type
    
    // set parameter values
    Parameters& params = this->get_equation_systems().parameters;
    
    unsigned int o = _infile("fe_order", 1);
    std::string fe_family = _infile("fe_family", std::string("LAGRANGE"));
    FEFamily fefamily = Utility::string_to_enum<FEFamily>(fe_family);

    // velocity potential
    vars[0]  = this->add_variable ( "phi", static_cast<Order>(o), fefamily);
    this->time_evolving(vars[0]);
    params.set<Real> ("u1_inf") =
    flight_condition->velocity_magnitude*flight_condition->drag_normal(0);
    params.set<Real> ("u2_inf") =
    flight_condition->velocity_magnitude*flight_condition->drag_normal(1);
    params.set<Real> ("u3_inf") =
    flight_condition->velocity_magnitude*flight_condition->drag_normal(2);
    
    // density
    vars[1] = this->add_variable ("rho", static_cast<Order>(o), fefamily);
    this->time_evolving(vars[1]);
    params.set<Real> ("rhou_inf") = flight_condition->rho();
    
    // Useful debugging options
    this->verify_analytic_jacobians = _infile("verify_analytic_jacobians", 0.);
    this->print_jacobians = _infile("print_jacobians", false);
    this->print_element_jacobians = _infile("print_element_jacobians", false);
    
    // Do the parent's initialization after variables and boundary constraints are defined
    FEMSystem::init_data();
}



void UnsteadyCompressiblePotentialFlow::init_context(DiffContext &context)
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
        c.get_side_fe( vars[i], elem_side_fe[i]);
        elem_side_fe[i]->get_JxW();
        elem_side_fe[i]->get_phi();
        elem_side_fe[i]->get_dphi();
        elem_side_fe[i]->get_xyz();
    }
}


bool UnsteadyCompressiblePotentialFlow::element_time_derivative (bool request_jacobian,
                                                                 DiffContext &context)
{
    FEMContext &c = libmesh_cast_ref<FEMContext&>(context);
    
    FEBase* elem_fe;
    
    c.get_element_fe(vars[0], elem_fe);
    
    // Element Jacobian * quadrature weights for interior integration
    const std::vector<Real> &JxW = elem_fe->get_JxW();
    
    // The number of local degrees of freedom in each variable
    unsigned int n_dofs = 0;
    for (unsigned int i=0; i<2; i++)
        n_dofs += c.get_dof_indices( vars[i] ).size();
    
    libmesh_assert_equal_to (n_dofs, 2*c.get_dof_indices( vars[0] ).size());
    
    // The subvectors and submatrices we need to fill:
    DenseMatrix<Number>& Kmat = c.get_elem_jacobian();
    DenseVector<Number>& Fvec = c.get_elem_residual();
    
    // Now we will build the element Jacobian and residual.
    // Constructing the residual requires the solution and its
    // gradient from the previous timestep.  This must be
    // calculated at each quadrature point by summing the
    // solution degree-of-freedom values by the appropriate
    // weight functions.
    const unsigned int n_qpoints = c.get_element_qrule().n_points(), n1 = 2;
    
    FEMOperatorMatrix B_mat;
    std::vector<FEMOperatorMatrix> dB_mat(dim);
    DenseMatrix<Real> LS_mat, tmp_mat1_n1n1, tmp_mat2_n1n2, tmp_mat3_n2n2, tmp_mat4;
    DenseVector<Real> tmp_vec1_n1, tmp_vec2_n2;
    
    LS_mat.resize(n1, n_dofs); tmp_mat1_n1n1.resize(n1, n1);
    tmp_mat2_n1n2.resize(n1, n_dofs); tmp_mat3_n2n2.resize(n_dofs, n_dofs);
    tmp_vec1_n1.resize(n1); tmp_vec2_n2.resize(n_dofs);
    
    Point uvec;
    Real rho;
    const Real uinf = flight_condition->velocity_magnitude,
    gamma = flight_condition->gas_property.gamma,
    ainf  = flight_condition->gas_property.a,
    rhoinf = flight_condition->gas_property.rho;
    
    const std::vector<std::vector<Real> >& phi = elem_fe->get_phi(); // assuming that all variables have the same interpolation
    const unsigned int n_phi = (unsigned int)phi.size();
    
    for (unsigned int qp=0; qp != n_qpoints; qp++) {
        // first update the variables at the current quadrature point
        this->update_solution_at_quadrature_point(vars, qp, c,
                                                  true, c.get_elem_solution(),
                                                  tmp_vec1_n1, uvec,
                                                  B_mat, dB_mat,
                                                  LS_mat);
        rho = tmp_vec1_n1(1); // get value from the interpolated sol
        
        for (unsigned int i_dim=0; i_dim<dim; i_dim++) {
            // convective: Galerkin
            dB_mat[i_dim].vector_mult(tmp_vec1_n1, c.get_elem_solution()); // dB/dx_i V
            tmp_vec1_n1(0) *= 0.5*uvec(i_dim);  // [u_i 0] dphi/dx_i  : potential
            tmp_vec1_n1(1)  = 0.;               // 0                  : density
            B_mat.vector_mult_transpose(tmp_vec2_n2, tmp_vec1_n1);  // B [C1] dB/dx_i V
            Fvec.add(-JxW[qp], tmp_vec2_n2);
            
            // convective: Least-Squares
            LS_mat.vector_mult_transpose(tmp_vec2_n2, tmp_vec1_n1);  // L*^T tau [C1] dB/dx_i V
            Fvec.add(-JxW[qp], tmp_vec2_n2);
            
            
            // diffusive : Galerkin
            dB_mat[i_dim].vector_mult(tmp_vec1_n1, c.get_elem_solution()); // dB/dx_i V
            tmp_vec1_n1(1) = rho*tmp_vec1_n1(0);               // [rho  0] dphi/dx_i  : density
            tmp_vec1_n1(0) *= 0.;                              // [0    0] 0          : potential
            dB_mat[i_dim].vector_mult_transpose(tmp_vec2_n2, tmp_vec1_n1);  // dB/dx_i [C2] dB/dx_i V
            Fvec.add(JxW[qp], tmp_vec2_n2);

            // diffusive: Least-Squares term neglected
            // TODO: this requires a 2nd order differential of the flux

        }

        
        // source term
        tmp_vec1_n1.zero();
        tmp_vec1_n1(0) = 0.5*pow(uinf,2) -
        pow(ainf,2)/(gamma-1.) * (pow(rho/rhoinf, gamma-1.) - 1.);
        B_mat.vector_mult_transpose(tmp_vec2_n2, tmp_vec1_n1);
        Fvec.add(JxW[qp], tmp_vec2_n2);
        
        
        if (request_jacobian && c.elem_solution_derivative)
        {
            for (unsigned int i_dim=0; i_dim<dim; i_dim++)
            {
                // convective: Galerkin
                tmp_mat1_n1n1.zero();
                tmp_mat1_n1n1(0,0) = 0.5*uvec(i_dim);
                dB_mat[i_dim].left_multiply(tmp_mat2_n1n2, tmp_mat1_n1n1);
                B_mat.right_multiply_transpose(tmp_mat3_n2n2,
                                               tmp_mat2_n1n2);// B^T C1 dB/dx_i
                Kmat.add(-JxW[qp], tmp_mat3_n2n2);
                
                // convective: Least-Squares
                LS_mat.get_transpose(tmp_mat4);
                tmp_mat4.right_multiply_transpose(tmp_mat2_n1n2);// LS^T C1 dB/dx_i
                Kmat.add(-JxW[qp], tmp_mat4);
                
                // diffusive: Galerkin
                tmp_mat1_n1n1.zero();
                tmp_mat1_n1n1(1,0) = rho;
                dB_mat[i_dim].left_multiply(tmp_mat2_n1n2, tmp_mat1_n1n1);
                dB_mat[i_dim].right_multiply_transpose(tmp_mat3_n2n2,
                                                       tmp_mat2_n1n2);// dB/dx_i^T C2 dB/dx_i
                Kmat.add(JxW[qp], tmp_mat3_n2n2);
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





bool UnsteadyCompressiblePotentialFlow::side_time_derivative (bool request_jacobian,
                                                              DiffContext &context)
{
    FEMContext &c = libmesh_cast_ref<FEMContext&>(context);
    
    // check for the boundary tags to check if the boundary condition needs to be applied to this element
    
    FluidBoundaryConditionType mechanical_bc_type;
    
    unsigned int n_mechanical_bc = 0;
    
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
                    
                case FAR_FIELD:
                    mechanical_bc_type = FAR_FIELD;
                    n_mechanical_bc++;
                    break;
            }
        }
    }
    
    // return if no boundary condition is applied
    if (n_mechanical_bc == 0)
        return request_jacobian;
    
    
    libmesh_assert_equal_to(n_mechanical_bc, 1); // make sure that only one is active
    
    const unsigned int n1 = 2,
    spatial_dim = this->get_mesh().spatial_dimension();
    
    // The number of local degrees of freedom for this element
    unsigned int n_dofs = 0;
    for (unsigned int i=0; i<2; i++)
        n_dofs += c.get_dof_indices( vars[i] ).size();
    
    libmesh_assert_equal_to (n_dofs, (2)*c.get_dof_indices( vars[0] ).size());
    
    FEBase * side_fe;
    c.get_side_fe(vars[0], side_fe); // assuming all variables have the same FE
    
    // Element Jacobian * quadrature weights for interior integration
    const std::vector<Real> &JxW = side_fe->get_JxW();
    
    // Physical location of the quadrature points
    const std::vector<Point>& qpoint = side_fe->get_xyz();
    
    // boundary normals
    const std::vector<Point>& face_normals = side_fe->get_normals();
    
    DenseMatrix<Number>& Kmat = c.get_elem_jacobian();
    DenseVector<Number>& Fvec = c.get_elem_residual();
    DenseMatrix<Real> LS_mat;
    DenseVector<Real> tmp_vec1_n1, tmp_vec2_n2;
    
    tmp_vec1_n1.resize(n1); tmp_vec2_n2.resize(n_dofs);
    LS_mat.resize(n1, n_dofs);
    
    FEMOperatorMatrix B_mat;
    std::vector<FEMOperatorMatrix> dB_mat(dim);

    Point uvec;
    Real rho;
    
    switch (mechanical_bc_type)
    // adiabatic and isothermal are handled in the no-slip wall BC
    {
        case SYMMETRY_WALL: // inviscid boundary condition
                            // vi ni = 0       (slip wall)
        {
            // since vi ni = 0, and the diffusive flux does not have a rho component,
            // nothing needs to be done here
        }
            break;
            
            
        case SLIP_WALL: // inviscid boundary condition
                        // vi ni = wi_dot (ni + dni) - ui dni   (moving slip wall with deformation)
        {
            Real ui_ni = 0.;
            DenseVector<Real> local_normal, surface_vel, dnormal;
            local_normal.resize(spatial_dim);
            
            for (unsigned int qp=0; qp<qpoint.size(); qp++) {
                // first update the variables at the current quadrature point
                this->update_solution_at_quadrature_point (vars, qp, c,
                                                           false, c.get_elem_solution(),
                                                           tmp_vec1_n1, uvec,
                                                           B_mat, dB_mat,
                                                           LS_mat);
                rho = tmp_vec1_n1(1); // get value from the interpolated sol
                
                // copy the surface normal
                for (unsigned int i_dim=0; i_dim<dim; i_dim++)
                    local_normal(i_dim) = face_normals[qp](i_dim);
                
                // now check if the surface deformation is defined and
                // needs to be applied through transpiration boundary
                // condition
                ui_ni = 0.;
                
                if (surface_motion) { // get the surface motion data
                    surface_motion->surface_velocity_time_domain
                    (this->time, qpoint[qp], face_normals[qp],
                     surface_vel, dnormal);
                    
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
                    for (unsigned int i_dim=0; i_dim<spatial_dim; i_dim++)
                        ui_ni -= uvec(i_dim) * dnormal(i_dim);
                }
                
                tmp_vec1_n1.zero();
                tmp_vec1_n1(1) = rho*ui_ni;
                
                B_mat.vector_mult_transpose(tmp_vec2_n2, tmp_vec1_n1);
                Fvec.add(-JxW[qp], tmp_vec2_n2);
                
                if ( request_jacobian && c.get_elem_solution_derivative() ) {
                    // since the velocity is independent of solution, Jacobian is zero
                }
            }
        }
            break;
            
            
        case FAR_FIELD:
            // conditions enforced are:
            // -- ui_ni = ui_inf ni + (ui_interior ni - ui_inf ni)
        {
            Real ui_ni = 0.;
            
            for (unsigned int qp=0; qp<qpoint.size(); qp++) {
                // first update the variables at the current quadrature point
                this->update_solution_at_quadrature_point (vars, qp, c,
                                                           false, c.get_elem_solution(),
                                                           tmp_vec1_n1, uvec,
                                                           B_mat, dB_mat,
                                                           LS_mat);
                rho = tmp_vec1_n1(1); // get value from the interpolated sol
                
                tmp_vec1_n1.zero();
                tmp_vec1_n1(1) = rho*ui_ni;
                
                B_mat.vector_mult_transpose(tmp_vec2_n2, tmp_vec1_n1);
                Fvec.add(-JxW[qp], tmp_vec2_n2);
                
                if ( request_jacobian && c.get_elem_solution_derivative() ) {
                    // since the velocity is independent of solution, Jacobian is zero
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


bool UnsteadyCompressiblePotentialFlow::mass_residual (bool request_jacobian,
                                                       DiffContext &context)
{
    FEMContext &c = libmesh_cast_ref<FEMContext&>(context);
    
    FEBase* elem_fe;
    
    c.get_element_fe(vars[0], elem_fe);
    
    
    // Element Jacobian * quadrature weights for interior integration
    const std::vector<Real> &JxW = elem_fe->get_JxW();
    
    // The number of local degrees of freedom in each variable
    unsigned int n_dofs = 0;
    for (unsigned int i=0; i<2; i++)
        n_dofs += c.get_dof_indices( vars[i] ).size();
    
    libmesh_assert_equal_to (n_dofs, 2*c.get_dof_indices( vars[0] ).size());
    
    // The subvectors and submatrices we need to fill:
    DenseMatrix<Number>& Kmat = c.get_elem_jacobian();
    DenseVector<Number>& Fvec = c.get_elem_residual();
    
    // Now we will build the element Jacobian and residual.
    // Constructing the residual requires the solution and its
    // gradient from the previous timestep.  This must be
    // calculated at each quadrature point by summing the
    // solution degree-of-freedom values by the appropriate
    // weight functions.
    const unsigned int n_qpoints = c.get_element_qrule().n_points(), n1 = 2;
    
    FEMOperatorMatrix B_mat;
    std::vector<FEMOperatorMatrix> dB_mat(dim);
    DenseMatrix<Real> LS_mat, tmp_mat1_n2n2, tmp_mat3;
    DenseVector<Real> tmp_vec1_n1, tmp_vec2_n2;
    LS_mat.resize(n1, n_dofs); tmp_mat1_n2n2.resize(n_dofs, n_dofs);
    tmp_vec1_n1.resize(n1); tmp_vec2_n2.resize(n_dofs);
    
    Point uvec;
    Real rho;
    
    for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
        // first update the variables at the current quadrature point
        this->update_solution_at_quadrature_point (vars, qp, c,
                                                   true, c.get_elem_fixed_solution(),
                                                   tmp_vec1_n1, uvec,
                                                   B_mat, dB_mat,
                                                   LS_mat);
        rho = tmp_vec1_n1(1); // get value from the interpolated sol
        
        
        // Galerkin contribution to velocity
        B_mat.vector_mult( tmp_vec1_n1, c.get_elem_solution() );
        B_mat.vector_mult_transpose(tmp_vec2_n2, tmp_vec1_n1);
        Fvec.add(JxW[qp], tmp_vec2_n2);
        
        // LS contribution to velocity
        B_mat.vector_mult( tmp_vec1_n1, c.get_elem_solution() );
        LS_mat.vector_mult_transpose(tmp_vec2_n2, tmp_vec1_n1);
        Fvec.add(JxW[qp], tmp_vec2_n2);
        
        
        if (request_jacobian && c.get_elem_solution_derivative())
        {
            // contribution from unsteady term
            // Galerkin contribution of velocity
            B_mat.right_multiply_transpose(tmp_mat1_n2n2, B_mat);
            Kmat.add(JxW[qp], tmp_mat1_n2n2);
            
            // LS contribution of velocity
            LS_mat.get_transpose(tmp_mat3);
            B_mat.left_multiply(tmp_mat1_n2n2, tmp_mat3); // LS^T tau Bmat
            Kmat.add(JxW[qp], tmp_mat1_n2n2);
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




void UnsteadyCompressiblePotentialFlow::update_solution_at_quadrature_point
( const std::vector<unsigned int>& vars, const unsigned int qp, FEMContext& c,
 const bool if_elem_domain, const DenseVector<Real>& elem_solution,
 DenseVector<Real>& conservative_sol, Point& uvec,
 FEMOperatorMatrix& B_mat, std::vector<FEMOperatorMatrix>& dB_mat,
 DenseMatrix<Real>& LS_mat)
{
    conservative_sol.zero();

    FEBase* fe;
    
    DenseVector<Real> phi_vals;
    
    unsigned int i_var = 0; // assuming that all FE are same

    if (if_elem_domain)
        c.get_element_fe(vars[i_var], fe);
    else
        c.get_side_fe(vars[i_var], fe);
    
    const std::vector<std::vector<Real> >& phi = fe->get_phi();
    const unsigned int n_phi = (unsigned int)phi.size();
    
    phi_vals.resize(n_phi);
    
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
        phi_vals(i_nd) = phi[i_nd][qp];
    
    B_mat.reinit(c.n_vars(), phi_vals); // initialize the operator matrix
    
    // since the velocity needs to be evaluated even for boundaries, the
    // shape function derivatives are always calculated
    const std::vector<std::vector<RealVectorValue> >& dphi =
    fe->get_dphi();
    
    for ( unsigned int i_dim=0; i_dim<dim; i_dim++ )
    {
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vals(i_nd) = dphi[i_nd][qp](i_dim);
        dB_mat[i_dim].reinit(c.n_vars(), phi_vals);
        dB_mat[i_dim].vector_mult(conservative_sol, elem_solution);
        uvec(i_dim) = conservative_sol(0); // u_i = dpotential/dx_i
    }
    
    B_mat.vector_mult( conservative_sol, elem_solution );
    
    DenseMatrix<Real> mat1, mat2;
    mat1.resize(2,2); mat2.resize(B_mat.m(), B_mat.n());
    LS_mat.zero();
    mat1(0,0) = 1.; mat1(1,1) = 1.;
    B_mat.left_multiply(mat2, mat1);
    LS_mat.add(1., mat2);
    
    mat1.zero();
    for (unsigned int i_dim=0; i_dim<dim; i_dim++) {
        mat1(0,0) = uvec(i_dim);
        dB_mat[i_dim].left_multiply(mat2, mat1);
        LS_mat.add(1., mat2);
    }
    
    // currently the diffusive term is neglected

    // calculate the tau matrix
    // calculate the dot product of velocity times gradient of shape function
    DenseVector<Real> dN; dN.resize(dim);
    Real tau = 0., h = 0, u_val = uvec.size();
    DenseVector<Real> u;
    for (unsigned int i_dim=0; i_dim<dim; i_dim++)
        u(i_dim) = uvec(i_dim);
    u.scale(1.0/u_val);
    
    for (unsigned int i_nodes=0; i_nodes<dphi.size(); i_nodes++)
    {
        // set value of shape function gradient
        for (unsigned int i_dim=0; i_dim<dim; i_dim++)
            dN(i_dim) = dphi[i_nodes][qp](i_dim);
        
        h += fabs(dN.dot(u));
    }
    
    h = 2.0/h;
    tau = h/u_val;
    
    // now scale with the tau matrix
    LS_mat.scale_row(0, tau);
}





