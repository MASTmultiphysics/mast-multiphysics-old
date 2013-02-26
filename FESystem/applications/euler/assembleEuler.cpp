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
#include "libmesh/kelly_error_estimator.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/uniform_refinement_estimator.h"
#include "libmesh/diff_solver.h"
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



PrimitiveSolution::PrimitiveSolution()
{
    this->zero();
}


void
PrimitiveSolution::zero()
{
    this->primitive_sol.zero();
    rho = 0.; u1 = 0.; u2 = 0.; u3 = 0.; T = 0.; p = 0.; a = 0.; e_tot = 0.; k = 0.;
    entropy = 0.; mach = 0.;
}


void PrimitiveSolution::init(const unsigned int dim, const DenseVector<Real> &conservative_sol, const Real cp, const Real cv)
{
    const unsigned int n1 = dim+2;
    const Real R = cp-cv, gamma = cp/cv;
    primitive_sol.resize(n1);
    
    rho = conservative_sol(0);
    primitive_sol(0) = rho;
    
    u1 = conservative_sol(1)/rho;
    primitive_sol(1) = u1;
    k += u1*u1;
    
    if (dim > 1)
    {
        u2 = conservative_sol(2)/rho;
        primitive_sol(2) = u2;
        k += u2*u2;
    }
    
    if (dim > 2)
    {
        u3 = conservative_sol(3)/rho;
        primitive_sol(3) = u3;
        k += u3*u3;
    }
    
    k *= 0.5;
    T = (conservative_sol(n1-1)/rho - k)/cv;
    primitive_sol(n1-1) = T;
    
    e_tot = cv*T+k;
    a = sqrt(gamma*R*T);
    p = R*T*rho;
    entropy = log(p/pow(rho,gamma));
    mach = sqrt(2.0*k/(gamma*R*T));
}


Real
PrimitiveSolution::c_pressure(const Real p0, const Real q0)
{
    return (p-p0)/q0;
}


void
PrimitiveSolution::print(std::ostream& out)
{
    out << "Primitive Solution:" << std::endl;
    primitive_sol.print(out);
    std::cout
    << std::setw(15) <<  " rho: " << rho << std::endl
    << std::setw(15) <<  " u1: " << u1 << std::endl
    << std::setw(15) <<  " u2: " << u2 << std::endl
    << std::setw(15) <<  " u3: " << u3 << std::endl
    << std::setw(15) <<  " mach: " << mach << std::endl
    << std::setw(15) <<  " a: " << a << std::endl
    << std::setw(15) <<  " T: " << T << std::endl
    << std::setw(15) <<  " p: " << p << std::endl
    << std::setw(15) <<  " e_tot: " << e_tot << std::endl
    << std::setw(15) <<  " k: " << k << std::endl
    << std::setw(15) <<  " entropy: " << entropy << std::endl << std::endl;
}



void EulerSystem::init_data ()
{
    this->use_fixed_solution = true;
    
    dim = this->get_mesh().mesh_dimension();
    const Real pi = acos(-1.);
    
    vars.resize(dim+2);
    
    // Check the input file for Reynolds number, application type,
    // approximation type
    GetPot infile("navier.in");
    aoa = infile("aoa",0.0);
    rho_inf = infile("rho",1.05);
    mach_inf = infile("mach",0.5);
    temp_inf = infile("temp",300.0);
    cp = infile("cp",1.003e3);
    cv = infile("cv",0.716e3);
    
    gamma = cp/cv;
    R = cp-cv;
    a_inf = sqrt(gamma*R*temp_inf);
    
    u1_inf = mach_inf*a_inf*cos(aoa*pi/180.0);
    u2_inf = mach_inf*a_inf*sin(aoa*pi/180.0);
    u3_inf = 0.0;
    q0_inf = 0.5*rho_inf*(u1_inf*u1_inf+u2_inf*u2_inf+u3_inf*u3_inf);
    p_inf = R*rho_inf*temp_inf;
    
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
    
    vars[dim+2] = this->add_variable ("rhoe", static_cast<Order>(o), fefamily);
    this->time_evolving(vars[dim+2]);
    params.set<Real> ("rhoe_inf") = rho_inf*temp_inf*cv + q0_inf;
    
    // Useful debugging options
    // Set verify_analytic_jacobians to 1e-6 to use
    this->verify_analytic_jacobians = infile("verify_analytic_jacobians", 0.);
    this->print_jacobians = infile("print_jacobians", false);
    this->print_element_jacobians = infile("print_element_jacobians", false);
    
    // Set Dirichlet boundary conditions
    const boundary_id_type top_id = (dim==3) ? 5 : 2;
    
    std::set<boundary_id_type> top_bdys;
    top_bdys.insert(top_id);
    
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
    
    // Physical location of the quadrature points
    const std::vector<Point>& qpoint = elem_fe->get_xyz();
    
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
    Fvec.zero();
    if (request_jacobian && c.elem_solution_derivative)
        Kmat.zero();

    PrimitiveSolution primitive_sol;
    
    for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
        // first update the variables at the current quadrature point
        this->update_solution_at_quadrature_point( qp, c, true, conservative_sol, primitive_sol, B_mat);
        this->update_jacobian_at_quadrature_point( qp, c, primitive_sol, dB_mat, Ai_advection, A_entropy, A_inv_entropy );
        
        Ai_Bi_advection.zero();
        
        for (unsigned int i_dim=0; i_dim<dim; i_dim++)
        {
            tmp_mat = Ai_advection[i_dim];
            tmp_mat.right_multiply(dB_mat[i_dim]);
            Ai_Bi_advection.add(1.0, tmp_mat);
        }
                
        this->calculate_differential_dperator_matrix( qp, c, true, primitive_sol, B_mat, dB_mat, Ai_advection, Ai_Bi_advection, A_inv_entropy, LS_mat, diff_val);
        
//        if (this->if_update_discont_values)
//            (*this->discontinuity_capturing_value)[i] = diff_val;
//        else
//            diff_val = (*this->discontinuity_capturing_value)[i];
        
        for (FESystemUInt i_dim=0; i_dim<dim; i_dim++)
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
            for (FESystemUInt i_dim=0; i_dim<dim; i_dim++)
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
//    Fvec.print(std::cout);
//    if (request_jacobian && c.elem_solution_derivative)
//        Kmat.print(std::cout);
    
    return request_jacobian;
}



bool EulerSystem::side_constraint (bool request_jacobian,
                                   DiffContext &context)
{
    FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

    // check for the boundary tags to check if the boundary condition needs to be applied to this element
    
    bool if_wall_bc = false, if_inf_bc = false;
    
    if ( this->get_mesh().boundary_info->has_boundary_id(c.elem, c.side, 1) ) // wall bc
        if_wall_bc = true;
    if ( this->get_mesh().boundary_info->has_boundary_id(c.elem, c.side, 2) ) // infinite bc
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

    Fvec.zero();
    if ( request_jacobian && c.get_elem_solution_derivative() )
        Kmat.zero();

    
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

        for (FESystemUInt qp=0; qp<qpoint.size(); qp++)
        {
            xini = face_normals[qp] * vel;
            
            // first update the variables at the current quadrature point
            this->update_solution_at_quadrature_point( qp, c, true, conservative_sol, p_sol, B_mat);
            
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
                tmp_mat1 = B_mat;
                tmp_mat1.right_multiply_transpose(tmp_mat2);
                Kmat.add(-JxW[qp], tmp_mat1);
            }
        }
    }

    
    if ( if_inf_bc )
    {
        for (FESystemUInt qp=0; qp<qpoint.size(); qp++)
        {
            // first update the variables at the current quadrature point
            this->update_solution_at_quadrature_point( qp, c, true, conservative_sol, p_sol, B_mat);

            this->calculate_advection_left_eigenvector_and_inverse_for_normal(p_sol, face_normals[qp], eig_val, l_eig_vec, l_eig_vec_inv_tr);
            
            // for all eigenalues that are less than 0, the characteristics are coming into the domain, hence,
            // evaluate them using the given solution.
            tmp_mat1 = l_eig_vec_inv_tr;
            for (FESystemDouble j=0; j<n1; j++)
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
            for (FESystemDouble j=0; j<n1; j++)
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
                for (FESystemDouble j=0; j<n1; j++)
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
//    c.elem->print_info();
//    std::cout << if_inf_bc << "  " << if_wall_bc << std::endl;
//    Fvec.print(std::cout);
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
    
    // Physical location of the quadrature points
    const std::vector<Point>& qpoint = elem_fe->get_xyz();
    
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
    Fvec.zero();
    if (request_jacobian && c.get_elem_solution_derivative())
        Kmat.zero();
    
    PrimitiveSolution primitive_sol;
    
    for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
        // first update the variables at the current quadrature point
        this->update_solution_at_quadrature_point( qp, c, false, conservative_sol, primitive_sol, B_mat);
        this->update_jacobian_at_quadrature_point( qp, c, primitive_sol, dB_mat, Ai_advection, A_entropy, A_inv_entropy );
        
        Ai_Bi_advection.zero();
        
        for (unsigned int i_dim=0; i_dim<dim; i_dim++)
        {
            tmp_mat = Ai_advection[i_dim];
            tmp_mat.right_multiply(dB_mat[i_dim]);
            Ai_Bi_advection.add(1.0, tmp_mat);
        }
        
        this->calculate_differential_dperator_matrix( qp, c, false, primitive_sol, B_mat, dB_mat, Ai_advection, Ai_Bi_advection, A_inv_entropy, LS_mat, diff_val);
        
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
//    c.elem->print_info();
//    Fvec.print(std::cout);
//    if (request_jacobian && c.elem_solution_derivative)
//        Kmat.print(std::cout);

    return request_jacobian;
}




void EulerSystem::postprocess()
{
    
}




void EulerSystem::get_infinity_vars( DenseVector<Real>& vars_inf )
{
    Real k = 0.0;
    vars_inf(0) = rho_inf;
    vars_inf(1) = rho_inf*u1_inf; k += u1_inf*u1_inf;
    if (dim > 1)
    {
        vars_inf(2) = rho_inf*u2_inf;
        k += u2_inf*u2_inf;
    }
    if (dim > 2)
    {
        vars_inf(3) = rho_inf*u3_inf;
        k += u3_inf*u3_inf;
    }
    vars_inf(dim+2-1) = rho_inf*(cv*temp_inf+0.5*k);
}




void EulerSystem::update_solution_at_quadrature_point( const unsigned int qp, FEMContext& c, const bool if_elem_time_derivative,
                                                        DenseVector<Real>& conservative_sol, PrimitiveSolution& primitive_sol,
                                                        DenseMatrix<Real>& B_mat)
{
    conservative_sol.zero();
    B_mat.zero();
    
    FEBase* fe;
    c.get_element_fe(vars[0], fe); // assuming that all variables have same interpolation

    const std::vector<std::vector<Real> >& phi = fe->get_phi();
    const unsigned int n_phi = phi.size();
    
    for ( unsigned int i_var=0; i_var<c.n_vars(); i_var++ )
    {
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            B_mat( i_var, i_var*n_phi+i_nd ) = phi[i_nd][qp];
    }

    
    if ( if_elem_time_derivative ) // forcing function calculation
        B_mat.vector_mult( conservative_sol, c.elem_solution );
    else // time derivative calculation
        B_mat.vector_mult( conservative_sol, c.elem_fixed_solution );
    
    primitive_sol.zero();
    primitive_sol.init(dim, conservative_sol, cp, cv);
}




void EulerSystem::update_jacobian_at_quadrature_point( const unsigned int qp, FEMContext& c, PrimitiveSolution& primitive_sol,
                                                       std::vector<DenseMatrix<Real> >& dB_mat,
                                                       std::vector<DenseMatrix<Real> >& Ai_advection,
                                                       DenseMatrix<Real>& A_entropy, DenseMatrix<Real>& A_inv_entropy )
{
    for ( unsigned int i_dim=0; i_dim<dim; i_dim++ )
        dB_mat[ i_dim ].zero();
    
    FEBase* fe;
    c.get_element_fe(vars[0], fe); // assuming that all variables have same interpolation
    
    const std::vector<std::vector<RealVectorValue> >& dphi = fe->get_dphi();
    const unsigned int n_phi = dphi.size();
    
    for ( unsigned int i_var=0; i_var<c.n_vars(); i_var++ )
        for ( unsigned int i_dim=0; i_dim<dim; i_dim++ )
            for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
                dB_mat[i_dim]( i_var, i_var*n_phi+i_nd ) = dphi[i_nd][qp](i_dim);
    
    for ( unsigned int i_dim=0; i_dim < dim; i_dim++)
        this->calculate_advection_flux_jacobian ( i_dim, primitive_sol, Ai_advection[i_dim] );
    
    this->calculate_entropy_variable_jacobian ( primitive_sol, A_entropy, A_inv_entropy );
}



void
EulerSystem::calculate_advection_flux(const unsigned int calculate_dim,
                                      const PrimitiveSolution& sol,
                                      DenseVector<Real>& flux)
{
    const unsigned int n1 = 2 + dim;
    
    const Real rho = sol.rho,
    u1 = sol.u1,
    u2 = sol.u2,
    u3 = sol.u3,
    p = sol.p,
    e_tot = sol.e_tot;
    
    flux.zero();
    
    // calculate the flux using given flow parameters at this point
    switch (calculate_dim)
    {
        case 0:
        {
            flux(0) =  rho * u1;
            flux(n1-1) =  u1 * (rho * e_tot + p);
            switch (dim)
            {
                case 3:
                    flux(3) =  rho * u1 * u3;
                case 2:
                    flux(2) =  rho * u1 * u2;
                case 1:
                    flux(1) =  rho * u1 * u1 + p;
            }
        }
            break;
            
        case 1:
        {
            flux(0) =  rho * u2;
            flux(n1-1) =  u2 * (rho * e_tot + p);
            switch (dim)
            {
                case 3:
                    flux(3) =  rho * u2 * u3;
                case 2:
                    flux(2) =  rho * u2 * u2 + p;
                case 1:
                    flux(1) =  rho * u2 * u1;
            }
        }
            break;
            
        case 2:
        {
            flux(0) =  rho * u3;
            flux(n1-1) =  u3 * (rho * e_tot + p);
            switch (dim)
            {
                case 3:
                    flux(3) =  rho * u3 * u3 + p;
                case 2:
                    flux(2) =  rho * u3 * u2;
                case 1:
                    flux(1) =  rho * u3 * u1;
            }
        }
            break;
    }
}




void
EulerSystem::calculate_advection_flux_jacobian(const unsigned int calculate_dim,
                                               const PrimitiveSolution& sol,
                                               DenseMatrix<Real>& mat)
{
    // calculate Ai = d F_adv / d x_i, where F_adv is the Euler advection flux vector
    
    const unsigned int n1 = 2 + dim;
    //    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    //
    //    FESystemAssert4((s.first == n1) && (s.second == n1), FESystem::Numerics::MatrixSizeMismatch, s.first, s.second, n1, n1);
    //    FESystemAssert0(div_coord < dim, FESystem::Exception::InvalidValue);
    
    mat.zero();
    
    const Real rho = sol.rho,
    u1 = sol.u1,
    u2 = sol.u2,
    u3 = sol.u3,
    k = sol.k,
    e_tot = sol.e_tot,
    T = sol.T;
    
    switch (calculate_dim)
    {
        case 0:
        {
            switch (dim)
            {
                case 3:
                {
                    mat(1, 3) = -u3*R/cv;
                    
                    mat(3, 0) = -u1*u3;
                    mat(3, 1) = u3;
                    mat(3, 3) = u1;
                    
                    mat(n1-1, 3) = -u1*u3*R/cv;
                }
                    
                case 2:
                {
                    mat(1, 2) = -u2*R/cv;
                    
                    mat(2, 0) = -u1*u2;
                    mat(2, 1) = u2;
                    mat(2, 2) = u1;
                    
                    mat(n1-1, 2) = -u1*u2*R/cv;
                }
                    
                case 1:
                {
                    mat(0, 1) = 1.0; // d U / d (rho u1)
                    
                    mat(1, 0) = -u1*u1+R*k/cv;
                    mat(1, 1) = u1*(2.0-R/cv);
                    mat(1, n1-1) = R/cv;
                    
                    mat(n1-1, 0) = u1*(R*(-e_tot+2.0*k)-e_tot*cv)/cv;
                    mat(n1-1, 1) = e_tot+R*(T-u1*u1/cv);
                    mat(n1-1, n1-1) = u1*gamma;
                }
                    break;
            }
        }
            break;
            
        case 1:
        {
            switch (dim)
            {
                case 3:
                {
                    mat(2, 3) = -u3*R/cv;
                    
                    mat(3, 0) = -u2*u3;
                    mat(3, 2) = u3;
                    mat(3, 3) = u2;
                    
                    mat(n1-1, 3) = -u2*u3*R/cv;
                }
                    
                case 2:
                {
                    mat(0, 2) = 1.0; // d U / d (rho u2)
                    
                    mat(1, 0) = -u1*u2;
                    mat(1, 1) = u2;
                    mat(1, 2) = u1;
                    
                    mat(2, 0) = -u2*u2+R*k/cv;
                    mat(2, 1) = -u1*R/cv;
                    mat(2, 2) = u2*(2.0-R/cv);
                    mat(2, n1-1) = R/cv;
                    
                    mat(n1-1, 0) = u2*(R*(-e_tot+2.0*k)-e_tot*cv)/cv;
                    mat(n1-1, 1) = -u1*u2*R/cv;
                    mat(n1-1, 2) = e_tot+R*(T-u2*u2/cv);
                    mat(n1-1, n1-1) = u2*gamma;
                }
                    break;
                    
                case 1:
                    // if second coordinate divergence is being asked for, then the element is atleast 2D
                    FESystemAssert0(false, FESystem::Exception::InvalidValue);
                    break;
            }
        }
            break;
            
        case 2:
        {
            mat(0, 3) = 1.0; // d U / d (rho u3)
            
            mat(1, 0) = -u1*u3;
            mat(1, 1) = u3;
            mat(1, 3) = u1;
            
            mat(2, 0) = -u2*u3;
            mat(2, 2) = u3;
            mat(2, 3) = u2;
            
            mat(3, 0) = -u3*u3+R*k/cv;
            mat(3, 1) = -u1*R/cv;
            mat(3, 2) = -u2*R/cv;
            mat(3, 3) = u3*(2.0-R/cv);
            mat(3, n1-1) = R/cv;
            
            mat(n1-1, 0) = u3*(R*(-e_tot+2.0*k)-e_tot*cv)/cv;
            mat(n1-1, 1) = -u1*u3*R/cv;
            mat(n1-1, 2) = -u2*u3*R/cv;
            mat(n1-1, 3) = e_tot+R*(T-u3*u3/cv);
            mat(n1-1, n1-1) = u3*gamma;
        }
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
    }
}


void
EulerSystem::calculate_advection_left_eigenvector_and_inverse_for_normal(const PrimitiveSolution& sol,
                                                                         const Point& normal, DenseMatrix<Real>& eig_vals,
                                                                         DenseMatrix<Real>& l_eig_mat, DenseMatrix<Real>& l_eig_mat_inv_tr)
{
    const unsigned int n1 = 2 + dim;
    
    //    const std::pair<FESystemUInt, FESystemUInt> s_eig_val = eig_vals.getSize(), s_l_eig_mat = l_eig_mat.getSize(), s_l_eig_mat_inv_tr = l_eig_mat_inv_tr.getSize();
    //
    //    FESystemAssert4((s_eig_val.first == n1) && (s_eig_val.first == n1), FESystem::Numerics::MatrixSizeMismatch, s_eig_val.first, s_eig_val.second, n1, n1);
    //    FESystemAssert4((s_l_eig_mat.first == n1) && (s_l_eig_mat.first == n1), FESystem::Numerics::MatrixSizeMismatch, s_l_eig_mat.first, s_l_eig_mat.second, n1, n1);
    //    FESystemAssert4((s_l_eig_mat_inv_tr.first == n1) && (s_l_eig_mat_inv_tr.first == n1), FESystem::Numerics::MatrixSizeMismatch, s_l_eig_mat_inv_tr.first, s_l_eig_mat_inv_tr.second, n1, n1);
    
    eig_vals.zero(); l_eig_mat.zero(); l_eig_mat_inv_tr.zero();
    
    Real nx=0., ny=0., nz=0., u=0.;
    unsigned int dim_for_eig_vec=100; // initializing with arbitrarily high value
    
    const Real rho = sol.rho,
    u1 = sol.u1,
    u2 = sol.u2,
    u3 = sol.u3,
    k = sol.k,
    T = sol.T,
    a = sol.a;
    
    // initialize the values
    switch (dim)
    {
        case 3:
        {
            nz = normal(2);
            u += u3*nz;
        }
            
        case 2:
        {
            ny = normal(1);
            u += u2*ny;
        }
            
        case 1:
        {
            nx = normal(0);
            u += u1*nx;
        }
    }
    
    // select the largest value of surface normal component, so that the appropriate section can be chosen
    if ((fabs(nx)>=fabs(ny)) && (fabs(nx)>=fabs(nz)))
        dim_for_eig_vec = 1;
    else if ((fabs(ny)>=fabs(nx)) && (fabs(ny)>=fabs(nz)))
        dim_for_eig_vec = 2;
    else if ((fabs(nz)>=fabs(nx)) && (fabs(nz)>=fabs(ny)))
        dim_for_eig_vec = 3;
    
    // set eigenvalues
    switch (dim)
    {
        case 3:
            eig_vals(2, 2) = u;
        case 2:
            eig_vals(1, 1) = u;
        case 1:
        {
            eig_vals(0, 0) = u;
            eig_vals(n1-2, n1-2) = u-a;
            eig_vals(n1-1, n1-1) = u+a;
        }
    }
    
    // set last two columns of the eigenvector matrices, and all columns of the eigenvector inverse.
    // Note that the dim_for_eig_vec column of the inverse matrix will be overwritten by the appropriate matrix
    switch (dim)
    {
        case 3:
        {
            // for u-a
            l_eig_mat(3, n1-2) = -u3-nz*a/(gamma-1.0);
            
            // for u+a
            l_eig_mat(3, n1-1) = -u3+nz*a/(gamma-1.0);
            
            // for u
            l_eig_mat_inv_tr(3, 0) = (gamma-1.0)*u1*u3/(a*a)-nx*nz;
            
            // for u
            l_eig_mat_inv_tr(3, 1) = (gamma-1.0)*u2*u3/(a*a)-ny*nz;
            
            // for u
            l_eig_mat_inv_tr(0, 2) = (gamma-1.0)*u3/(a*a);
            l_eig_mat_inv_tr(1, 2) = (gamma-1.0)*u3*u1/(a*a)-nz*nx;
            l_eig_mat_inv_tr(2, 2) = (gamma-1.0)*u3*u2/(a*a)-nz*ny;
            l_eig_mat_inv_tr(3, 2) = (gamma-1.0)*u3*u3/(a*a)+1.0-nz*nz;
            l_eig_mat_inv_tr(n1-1, 2) = (gamma-1.0)*k*u3/(a*a)+u3-nz*u;
            
            // overwrite for u
            l_eig_mat_inv_tr(3, dim_for_eig_vec-1) = u3;
            
            // for u-a
            l_eig_mat_inv_tr(3, n1-2) = 0.5*(gamma-1.0)*(-nz+u3/a)/a;
            
            // for u+a
            l_eig_mat_inv_tr(3, n1-1) = 0.5*(gamma-1.0)*(nz+u3/a)/a;
            
        }
            
        case 2:
        {
            // for u-a
            l_eig_mat(2, n1-2) = -u2-ny*a/(gamma-1.0);
            
            // for u+a
            l_eig_mat(2, n1-1) = -u2+ny*a/(gamma-1.0);
            
            // for u
            l_eig_mat_inv_tr(2, 0) = (gamma-1.0)*u1*u2/(a*a)-nx*ny;
            
            // for u
            l_eig_mat_inv_tr(0, 1) = (gamma-1.0)*u2/(a*a);
            l_eig_mat_inv_tr(1, 1) = (gamma-1.0)*u2*u1/(a*a)-ny*nx;
            l_eig_mat_inv_tr(2, 1) = (gamma-1.0)*u2*u2/(a*a)+1.0-ny*ny;
            l_eig_mat_inv_tr(n1-1, 1) = (gamma-1.0)*k*u2/(a*a)+u2-ny*u;
            
            // overwrite for u
            l_eig_mat_inv_tr(2, dim_for_eig_vec-1) = u2;
            
            // for u-a
            l_eig_mat_inv_tr(2, n1-2) = 0.5*(gamma-1.0)*(-ny+u2/a)/a;
            
            // for u+a
            l_eig_mat_inv_tr(2, n1-1) = 0.5*(gamma-1.0)*(ny+u2/a)/a;
        }
            
        case 1:
        {
            // for u-a
            l_eig_mat(0, n1-2) = k+u*a/(gamma-1.0);
            l_eig_mat(1, n1-2) = -u1-nx*a/(gamma-1.0);
            l_eig_mat(n1-1, n1-2) = 1.0;
            
            // for u+a
            l_eig_mat(0, n1-1) = k-u*a/(gamma-1.0);
            l_eig_mat(1, n1-1) = -u1+nx*a/(gamma-1.0);
            l_eig_mat(n1-1, n1-1) = 1.0;
            
            // for u
            l_eig_mat_inv_tr(0, 0) = (gamma-1.0)*u1/(a*a);
            l_eig_mat_inv_tr(1, 0) = (gamma-1.0)*u1*u1/(a*a)+1.0-nx*nx;
            l_eig_mat_inv_tr(n1-1, 0) = (gamma-1.0)*k*u1/(a*a)+u1-nx*u;
            
            // overwrite for u
            l_eig_mat_inv_tr(0, dim_for_eig_vec-1) = 1.0;
            l_eig_mat_inv_tr(1, dim_for_eig_vec-1) = u1;
            l_eig_mat_inv_tr(n1-1, dim_for_eig_vec-1) = k;
            l_eig_mat_inv_tr.scale_column(dim_for_eig_vec-1, -1.0/(cv*T*gamma));
            
            // for u-a
            l_eig_mat_inv_tr(0, n1-2) = 1.0/(2.0*cv*T*gamma);
            l_eig_mat_inv_tr(1, n1-2) = 0.5*(gamma-1.0)*(-nx+u1/a)/a;
            l_eig_mat_inv_tr(n1-1, n1-2) = 0.5*(1.0+(gamma-1.0)*(-u+k/a)/a);
            
            // for u+a
            l_eig_mat_inv_tr(0, n1-1) = 1.0/(2.0*cv*T*gamma);
            l_eig_mat_inv_tr(1, n1-1) = 0.5*(gamma-1.0)*(nx+u1/a)/a;
            l_eig_mat_inv_tr(n1-1, n1-1) = 0.5*(1.0+(gamma-1.0)*(u+k/a)/a);
            
        }
    }
    
    
    switch (dim_for_eig_vec)
    {
        case 1:
        {
            // set values in the left eigenvector matrix and eigenvalue matrix
            switch (dim)
            {
                case 3:
                {
                    // for u
                    l_eig_mat(0, 2) = nz*u1/nx-u3;
                    l_eig_mat(1, 2) = -nz/nx;
                    l_eig_mat(3, 2) = 1.0;
                }
                    
                case 2:
                {
                    // for u
                    l_eig_mat(0, 1) = ny*u1/nx-u2;
                    l_eig_mat(1, 1) = -ny/nx;
                    l_eig_mat(2, 1) = 1.0;
                }
                    
                case 1:
                {
                    // for u
                    l_eig_mat(0, 0) = u*u1/nx-(cv*T*gamma+k);
                    l_eig_mat(1, 0) = -u/nx;
                    l_eig_mat(n1-1, 0) = 1.0;
                }
            }
        }
            break;
            
        case 2:
        {
            // set values in the left eigenvector matrix and eigenvalue matrix
            switch (dim)
            {
                case 3:
                {
                    // for u
                    l_eig_mat(0, 2) = nz*u2/ny-u3;
                    l_eig_mat(2, 2) = -nz/ny;
                    l_eig_mat(3, 2) = 1.0;
                }
                    
                case 2:
                case 1:
                {
                    // for u
                    l_eig_mat(0, 0) = nx*u2/ny-u1;
                    l_eig_mat(1, 0) = 1.0;
                    l_eig_mat(2, 0) = -nx/ny;
                    
                    // for u
                    l_eig_mat(0, 1) = u*u2/ny-(cv*T*gamma+k);
                    l_eig_mat(2, 1) = -u/ny;
                    l_eig_mat(n1-1, 1) = 1.0;
                }
            }
        }
            break;
            
        case 3:
        {
            // set values in the left eigenvector matrix and eigenvalue matrix
            
            // for u
            l_eig_mat(0, 0) = nx*u3/nz-u1;
            l_eig_mat(1, 0) = 1.0;
            l_eig_mat(3, 0) = -nx/nz;
            
            // for u
            l_eig_mat(0, 1) = ny*u3/nz-u2;
            l_eig_mat(2, 1) = 1.0;
            l_eig_mat(3, 1) = -ny/nz;
            
            // for u
            l_eig_mat(0, 2) = u*u3/nz-(cv*T*gamma+k);
            l_eig_mat(3, 2) = -u/nz;
            l_eig_mat(n1-1, 2) = 1.0;
        }
            break;
    }
}



void
EulerSystem::calculate_advection_flux_jacobian_for_moving_solid_wall_boundary(const PrimitiveSolution& sol,
                                                                              const Real xi_ni, const Point& nvec, DenseMatrix<Real>& mat)
{
    // calculate Ai = d F_adv / d x_i, where F_adv is the Euler advection flux vector
    
    const unsigned int n1 = 2 + dim;
    //    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    //
    //    FESystemAssert4((s.first == n1) && (s.second == n1), FESystem::Numerics::MatrixSizeMismatch, s.first, s.second, n1, n1);
    
    mat.zero();
    const Real u1 = sol.u1,
    u2 = sol.u2,
    u3 = sol.u3,
    k = sol.k;
    
    switch (dim)
    {
        case 3:
        {
            mat(1, 3) = -R/cv*nvec(0)*u3;
            
            mat(2, 3) = -R/cv*nvec(1)*u3;
            
            mat(3, 0) = nvec(2)*R/cv*k;
            mat(3, 1) = -R/cv*nvec(2)*u1;
            mat(3, 2) = -R/cv*nvec(2)*u2;
            mat(3, 3) = xi_ni-R/cv*nvec(2)*u3;
            mat(3, n1-1) = R/cv*nvec(2);
            
            mat(n1-1, 3) = xi_ni*R/cv*u3;
        }
            
        case 2:
        {
            mat(1, 2) = -R/cv*nvec(0)*u2;
            
            mat(2, 0) = nvec(1)*R/cv*k;
            mat(2, 1) = -R/cv*nvec(1)*u1;
            mat(2, 2) = xi_ni-R/cv*nvec(1)*u2;
            mat(2, n1-1) = R/cv*nvec(1);
            
            mat(n1-1, 2) = xi_ni*R/cv*u2;
        }
            
        case 1:
        {
            mat(0, 1) = xi_ni; // d U / d (rho u1)
            
            mat(1, 0) = nvec(0)*R/cv*k;
            mat(1, 1) = xi_ni-R/cv*nvec(0)*u1;
            mat(1, n1-1) = R/cv*nvec(0);
            
            mat(n1-1, 0) = xi_ni*R/cv*k;
            mat(n1-1, 1) = xi_ni*R/cv*u1;
            mat(n1-1, n1-1) = xi_ni*(R+cv)/cv;
        }
            break;
    }
}




void
EulerSystem::calculate_entropy_variable_jacobian(const PrimitiveSolution& sol,
                                                 DenseMatrix<Real>& dUdV, DenseMatrix<Real>& dVdU)
{
    // calculates dU/dV where V is the Entropy variable vector
    
    // calculate A0 = d U / d Y , where U = conservation variable vector, and Y = unknown variable vector
    // note that for conservation variables as the unknown, this is an identity matrix
    
    const unsigned int n1 = 2 + dim;
    const Real rho = sol.rho,
    u1 = sol.u1,
    u2 = sol.u2,
    u3 = sol.u3,
    k = sol.k,
    T = sol.T,
    e_tot = sol.e_tot;
    
    dUdV.zero(); dVdU.zero();
    
    // du/dv
    switch (dim)
    {
        case 3:
        {
            dUdV(0, 3) = u3;
            
            dUdV(1, 3) = u1*u3;
            
            dUdV(2, 3) = u2*u3;
            
            dUdV(3, 0) = dUdV(0, 3);
            dUdV(3, 1) = dUdV(1, 3);
            dUdV(3, 2) = dUdV(2, 3);
            dUdV(3, 3) = u3*u3+cv*T*(gamma-1.0);
            dUdV(3, n1-1) = u3*(cv*T*gamma+k);
            
            dUdV(n1-1, 3) = dUdV(3, n1-1);
        }
            
        case 2:
        {
            dUdV(0, 2) = u2;
            
            dUdV(1, 2) = u1*u2;
            
            dUdV(2, 0) = dUdV(0, 2);
            dUdV(2, 1) = dUdV(1, 2);
            dUdV(2, 2) = u2*u2+cv*T*(gamma-1.0);
            dUdV(2, n1-1) = u2*(cv*T*gamma+k);
            
            dUdV(n1-1, 2) = dUdV(2, n1-1);
        }
            
        case 1:
        {
            dUdV(0, 0) = 1.0;
            dUdV(0, 1) = u1;
            dUdV(0, n1-1) = e_tot;
            
            dUdV(1, 0) = dUdV(0, 1);
            dUdV(1, 1) = u1*u1+cv*T*(gamma-1.0);
            dUdV(1, n1-1) = u1*(cv*T*gamma+k);
            
            dUdV(n1-1, 0) = dUdV(0, n1-1);
            dUdV(n1-1, 1) = dUdV(1, n1-1);
            dUdV(n1-1, n1-1) = k*k+gamma*cv*T*(cv*T+2*k);
            
        }
            break;
            
        default:
            break;
    }
    
    dUdV.scale(rho/(gamma-1.0));
    
    
    // dv/du
    switch (dim)
    {
        case 3:
        {
            dVdU(0, 3) = -u3*k;
            
            dVdU(1, 3) = u1*u3;
            
            dVdU(2, 3) = u2*u3;
            
            dVdU(3, 0) = dVdU(0, 3);
            dVdU(3, 1) = dVdU(1, 3);
            dVdU(3, 2) = dVdU(2, 3);
            dVdU(3, 3) = u3*u3+cv*T;
            dVdU(3, n1-1) = -u3;
            
            dVdU(n1-1, 3) = dVdU(3, n1-1);
        }
            
        case 2:
        {
            dVdU(0, 2) = -u2*k;
            
            dVdU(1, 2) = u1*u2;
            
            dVdU(2, 0) = dVdU(0, 2);
            dVdU(2, 1) = dVdU(1, 2);
            dVdU(2, 2) = u2*u2+cv*T;
            dVdU(2, n1-1) = -u2;
            
            dVdU(n1-1, 2) = dVdU(2, n1-1);
        }
            
        case 1:
        {
            dVdU(0, 0) = k*k+cv*cv*T*T*gamma;
            dVdU(0, 1) = -u1*k;
            dVdU(0, n1-1) = -e_tot+2.0*k;
            
            dVdU(1, 0) = dVdU(0, 1);
            dVdU(1, 1) = u1*u1+cv*T;
            dVdU(1, n1-1) = -u1;
            
            dVdU(n1-1, 0) = dVdU(0, n1-1);
            dVdU(n1-1, 1) = dVdU(1, n1-1);
            dVdU(n1-1, n1-1) = 1.0;
        }
            break;
            
        default:
            break;
    }
    
    dVdU.scale(1.0/(rho*cv*cv*T*T));
}




void
EulerSystem::calculate_artificial_diffusion_operator(const unsigned int qp, FEMContext& c,
                                                     const PrimitiveSolution& sol,
                                                     DenseMatrix<Real>& streamline_operator)
{
    const unsigned int n1 = 2 + dim;
    
    FEBase* fe;
    c.get_element_fe(vars[0], fe);
    std::vector<std::vector<RealVectorValue> > dphi = fe->get_dphi(); // assuming that all variables have the same interpolation
    
    DenseVector<Real> u, dN;
    u.resize(dim); dN.resize(dim);
    
    const Real u1 = sol.u1,
    u2 = sol.u2,
    u3 = sol.u3,
    a = sol.a,
    dt = c.get_deltat_value();
    
    streamline_operator.zero();
    
    // calculate the gradients
    switch (dim)
    {
        case 3:
            u(2) = u3;
            
        case 2:
            u(1) = u2;
            
        case 1:
            u(0) = u1;
            break;
            
        default:
            break;
    }
    
    // calculate the dot product of velocity times gradient of shape function
    Real h = 0, u_val = u.l2_norm(), tau_rho, tau_m, tau_e;
    u.scale(1.0/u_val);
    
    for (FESystemUInt i_nodes=0; i_nodes<dphi.size(); i_nodes++)
    {
        // set value of shape function gradient
        for (FESystemUInt i_dim=0; i_dim<dim; i_dim++)
            dN(i_dim) = dphi[i_nodes][qp](i_dim);
        
        h += fabs(dN.dot(u));
    }
    
    h = 2.0/h;
    
    // now set the value of streamwise dissipation
    tau_rho = 1.0/sqrt(pow(2.0/dt, 2)+ pow(2.0/h*(u_val+a), 2));
    tau_m = 1.0/sqrt(pow(2.0/dt, 2)+ pow(2.0/h*(u_val+a), 2));
    tau_e = 1.0/sqrt(pow(2.0/dt, 2)+ pow(2.0/h*(u_val+a), 2));
    
    streamline_operator(0, 0) = tau_rho;
    for (FESystemUInt i_dim=0; i_dim<dim; i_dim++)
        streamline_operator(1+i_dim, 1+i_dim) = tau_m;
    streamline_operator(n1-1, n1-1) = tau_e;
}


void
EulerSystem::calculate_dxidX (const unsigned int qp, FEMContext& c,
                              DenseMatrix<Real>& dxi_dX)
{
    // initialize dxi_dX
    dxi_dX.zero();
    Real val=0.;
    FEBase* fe;
    c.get_element_fe(vars[0], fe); // assuming that all elements have the same interpolation fe
    
    for (unsigned int i_dim=0; i_dim<dim; i_dim++)
        for (unsigned int j_dim=0; j_dim<dim; j_dim++)
        {
            switch (i_dim)
            {
                case 0:
                {
                    switch (j_dim)
                    {
                        case 0:
                            val = fe->get_dxidx()[qp];
                            break;
                        case 1:
                            val = fe->get_dxidy()[qp];
                            break;
                        case 2:
                            val = fe->get_dxidz()[qp];
                            break;
                    }
                }
                    break;
                case 1:
                {
                    switch (j_dim)
                    {
                        case 0:
                            val = fe->get_detadx()[qp];
                            break;
                        case 1:
                            val = fe->get_detady()[qp];
                            break;
                        case 2:
                            val = fe->get_detadz()[qp];
                            break;
                    }
                }
                    break;
                case 2:
                {
                    switch (j_dim)
                    {
                        case 0:
                            val = fe->get_dzetadx()[qp];
                            break;
                        case 1:
                            val = fe->get_dzetady()[qp];
                            break;
                        case 2:
                            val = fe->get_dzetadz()[qp];
                            break;
                    }
                }
                    break;
            }
            dxi_dX(i_dim, j_dim) = val;
        }
}


void EulerSystem::calculate_differential_dperator_matrix( const unsigned int qp, FEMContext& c, const bool if_elem_time_derivative, const PrimitiveSolution& sol,
                                                          const DenseMatrix<Real>& B_mat, const std::vector<DenseMatrix<Real> >& dB_mat,
                                                          const std::vector<DenseMatrix<Real> >& Ai_advection, const DenseMatrix<Real>& Ai_Bi_advection,
                                                          const DenseMatrix<Real>& A_inv_entropy,
                                                          DenseMatrix<Real>& LS_operator, Real& discontinuity_val )
{
    const unsigned int n1 = 2 + dim, n2 = B_mat.n();
    
    std::vector<DenseVector<Real> > diff_vec(3);
    DenseMatrix<Real> tmp_mat, tmp_mat_n1n1, diff_operator, dxi_dX;
    DenseVector<Real> vec1, vec2;
    tmp_mat.resize(n1, n2); tmp_mat_n1n1.resize(n1, n1); diff_operator.resize(n1, n2);
    dxi_dX.resize(dim, dim);
    vec1.resize(n1); vec2.resize(n1);
    for (FESystemUInt i=0; i<dim; i++) diff_vec[i].resize(n1);
    
    this->calculate_dxidX (qp, c, dxi_dX);
    
    // contribution of unsteady term
    LS_operator.zero();
    
    diff_operator = B_mat;
    diff_operator += Ai_Bi_advection;
    
    FESystemDouble val1 = 0.0;
    
    vec2.zero();
    
    // contribution of advection flux term
    for (FESystemUInt i=0; i<dim; i++)
    {
        Ai_advection[i].get_transpose(tmp_mat);
        tmp_mat.right_multiply(dB_mat[i]);
        LS_operator.add(1.0, tmp_mat);  // A_i^T B_i
        
        if (if_elem_time_derivative)
            dB_mat[i].vector_mult(diff_vec[i], c.elem_solution);
        else
            dB_mat[i].vector_mult(diff_vec[i], c.elem_fixed_solution);
        
        Ai_advection[i].vector_mult(vec1, diff_vec[i]);
        
        vec2.add(1.0, vec1); // sum A_i dU/dx_i
    }
    
    // add the velocity and calculate the numerator of the discontinuity capturing term coefficient
    //vec2 += c.elem_solution; // add velocity TODO: how to get the velocity for all calls to this method
    A_inv_entropy.vector_mult(vec1, vec2);
    discontinuity_val = vec1.dot(vec2); // this is the numerator term
    
    // now evaluate the dissipation factor for the discontinuity capturing term
    // this is the denominator term
    
    val1 = 0.0;
    for (FESystemUInt i=0; i<dim; i++)
    {
        vec1.zero();
        tmp_mat.zero();
        
        for (FESystemUInt j=0; j<dim; j++)
            vec1.add(dxi_dX(i, j), diff_vec[j]);
        
        // calculate the value of denominator
        A_inv_entropy.vector_mult(vec2, vec1);
        val1 += vec1.dot(vec2);
    }
    
    //    // now calculate the discontinuity capturing operator
    if ((fabs(val1) > 0.0) &&  (fabs(discontinuity_val) > 0.0))
        discontinuity_val = sqrt(discontinuity_val/val1);
    else
        discontinuity_val = 0.0;
    
    // scale the LS matrix with the correct factor
    this->calculate_artificial_diffusion_operator(qp, c, sol, tmp_mat_n1n1);
    LS_operator.left_multiply(tmp_mat_n1n1);
}


// The main program.
int main (int argc, char** argv)
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
    const unsigned int coarsegridsize    = infile("coarsegridsize", 1);
    const unsigned int coarserefinements = infile("coarserefinements", 0);
    const unsigned int max_adaptivesteps = infile("max_adaptivesteps", 10);
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
    
    // Use the MeshTools::Generation mesh generator to create a uniform
    // grid on the square [-1,1]^D.  We instruct the mesh generator
    // to build a mesh of 8x8 \p Quad9 elements in 2D, or \p Hex27
    // elements in 3D.  Building these higher-order elements allows
    // us to use higher-order approximation, as in example 3.
    if (dim == 2)
        MeshTools::Generation::build_square (mesh,
                                             coarsegridsize,
                                             coarsegridsize,
                                             0., 1.,
                                             0., 1.,
                                             QUAD4);
    else if (dim == 3)
        MeshTools::Generation::build_cube (mesh,
                                           coarsegridsize,
                                           coarsegridsize,
                                           coarsegridsize,
                                           0., 1.,
                                           0., 1.,
                                           0., 1.,
                                           HEX8);
    
    mesh_refinement.uniformly_refine(coarserefinements);
    
    // Print information about the mesh to the screen.
    mesh.print_info();
    
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
    
    // Set the time stepping options
    system.deltat = deltat;
    
    // And the nonlinear solver options
    DiffSolver &solver = *(system.time_solver->diff_solver().get());
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
    
    // And the linear solver options
    solver.max_linear_iterations =
    infile("max_linear_iterations", 50000);
    solver.initial_linear_tolerance =
    infile("initial_linear_tolerance", 1.e-3);
    
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
            std::vector<Real> weights(2,1.0);  // u, v
            if (dim == 3)
                weights.push_back(1.0);          // w
            weights.push_back(0.0);            // p
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
