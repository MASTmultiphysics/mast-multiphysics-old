//
//  shock_tube_fluid_elem.cpp
//  MAST
//
//  Created by Manav Bhatia on 5/16/14.
//  Copyright (c) 2014 Manav Bhatia. All rights reserved.
//

// C++ includes
#include <math.h>

// MAST includes
#include "FluidElems/shock_tube_fluid_elem.h"

// libMesh includes
#include "libmesh/quadrature.h"

//#ifndef LIBMESH_USE_COMPLEX_NUMBERS


MAST::ShockTubeFluidElem::ShockTubeFluidElem(libMesh::EquationSystems& es,
                                             const std::string& name_in,
                                             const unsigned int number_in):
FluidSystem(es, name_in, number_in) {
    
}


MAST::ShockTubeFluidElem::~ShockTubeFluidElem() {
    
}

bool
MAST::ShockTubeFluidElem::element_time_derivative(bool request_jacobian,
                                                  libMesh::DiffContext& context) {
    
    libmesh_assert(dim == 1);
    
    FEMContext &c = libmesh_cast_ref<FEMContext&>(context);
    
    libMesh::FEBase* elem_fe;
    
    c.get_element_fe(vars[0], elem_fe);
    
    // Element Jacobian * quadrature weights for interior integration
    const std::vector<libMesh::Real> &JxW = elem_fe->get_JxW();
    
    // The subvectors and submatrices we need to fill:
    libMesh::DenseMatrix<libMesh::Real>& Kmat = c.get_elem_jacobian();
    libMesh::DenseVector<libMesh::Real>& Fvec = c.get_elem_residual();
    
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
    std::vector<libMesh::DenseMatrix<libMesh::Real> >  Ai_advection(dim);
    libMesh::DenseMatrix<libMesh::Real> LS_mat, LS_sens, Ai_Bi_advection, tmp_mat_n1n1,
    tmp_mat_n1n2, tmp_mat2_n2n2, tmp_mat3, A_sens, stress_tensor,
    dprim_dcons, dcons_dprim;
    
    libMesh::DenseVector<libMesh::Real> flux, tmp_vec1_n1, tmp_vec2_n1, tmp_vec3_n2,
    conservative_sol, delta_vals, diff_val, temp_grad;
    
    LS_mat.resize(n1, n_dofs); LS_sens.resize(n_dofs, n_dofs);
    Ai_Bi_advection.resize(dim+2, n_dofs);
    tmp_mat_n1n2.resize(dim+2, n_dofs); A_sens.resize(n1, n_dofs);
    stress_tensor.resize(dim, dim); dprim_dcons.resize(dim+2, dim+2);
    dcons_dprim.resize(dim+2, dim+2);
    tmp_mat2_n2n2.resize(n_dofs, n_dofs);
    
    flux.resize(n1); tmp_vec1_n1.resize(n1); tmp_vec2_n1.resize(n1);
    tmp_vec3_n2.resize(n_dofs); conservative_sol.resize(dim+2);
    delta_vals.resize(n_qpoints); temp_grad.resize(dim); diff_val.resize(dim);
    
    for (unsigned int i=0; i<dim; i++)
        Ai_advection[i].resize(dim+2, dim+2);
    
    std::vector<std::vector<libMesh::DenseMatrix<libMesh::Real> > > flux_jacobian_sens;
    flux_jacobian_sens.resize(dim);
    for (unsigned int i_dim=0; i_dim<dim; i_dim++)
    {
        flux_jacobian_sens[i_dim].resize(n1); // number of variables for sensitivity
        for (unsigned int i_cvar=0; i_cvar<n1; i_cvar++)
            flux_jacobian_sens[i_dim][i_cvar].resize(n1, n1);
    }
    
    
    
    PrimitiveSolution primitive_sol;
    const std::vector<std::vector<libMesh::Real> >& phi = elem_fe->get_phi(); // assuming that all variables have the same interpolation
    const unsigned int n_phi = phi.size();
    libMesh::Real area, darea_dx, x;
    
    
    for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
        x =  elem_fe->get_xyz()[qp](0);
        area     = 1.398 + .347 * tanh(0.8 * x - 4.);
        darea_dx = .347 * (1. - pow(tanh(0.8 * x - 4.),2)) * 0.8;
     
        // first update the variables at the current quadrature point
        this->update_solution_at_quadrature_point
        (vars, qp, c,
         true, c.get_elem_solution(),
         conservative_sol, primitive_sol,
         B_mat, dB_mat);
        
        for ( unsigned int i_dim=0; i_dim < dim; i_dim++)
            this->calculate_advection_flux_jacobian ( i_dim, primitive_sol,
                                                     Ai_advection[i_dim] );
        
        Ai_Bi_advection.zero();
        
        for (unsigned int i_dim=0; i_dim<dim; i_dim++)
        {
            dB_mat[i_dim].left_multiply(tmp_mat_n1n2, Ai_advection[i_dim]);
            Ai_Bi_advection.add(1.0, tmp_mat_n1n2);
        }
        
        for (unsigned int i_dim=0; i_dim<dim; i_dim++)
            this->calculate_advection_flux_jacobian_sensitivity_for_conservative_variable
            (i_dim, primitive_sol, flux_jacobian_sens[i_dim]);
        
        this->calculate_differential_operator_matrix
        (vars, qp, c,
         c.get_elem_solution(), primitive_sol,
         B_mat, dB_mat, Ai_advection,
         Ai_Bi_advection, flux_jacobian_sens,
         LS_mat, LS_sens);
        
        this->calculate_aliabadi_discontinuity_operator(vars, qp, c,
                                                        primitive_sol,
                                                        c.get_elem_solution(),
                                                        dB_mat,
                                                        Ai_Bi_advection,
                                                        diff_val);
        
        for (unsigned int i_dim=0; i_dim<dim; i_dim++)
        {
            // Galerkin contribution from the advection flux terms
            this->calculate_advection_flux(i_dim, primitive_sol, flux); // F^adv_i
            dB_mat[i_dim].vector_mult_transpose(tmp_vec3_n2, flux); // dBw/dx_i F^adv_i
            Fvec.add(JxW[qp], tmp_vec3_n2);
            
            // discontinuity capturing operator
            dB_mat[i_dim].vector_mult(flux, c.get_elem_solution());
            dB_mat[i_dim].vector_mult_transpose(tmp_vec3_n2, flux);
            Fvec.add(-JxW[qp]*diff_val(i_dim), tmp_vec3_n2);
        }
        
        // Least square contribution from divergence of advection flux
        Ai_Bi_advection.vector_mult(flux, c.get_elem_solution()); // d F^adv_i / dxi
        LS_mat.vector_mult_transpose(tmp_vec3_n2, flux); // LS^T tau F^adv_i
        Fvec.add(-JxW[qp], tmp_vec3_n2);
        
        
        // Galerkin contribution from the source term
        this->calculate_source_flux(primitive_sol, flux); // S
        B_mat.vector_mult_transpose(tmp_vec3_n2, flux);   // Bw S
        Fvec.add(darea_dx/area*JxW[qp], tmp_vec3_n2);
        
        // Least square contribution from source term
        LS_mat.vector_mult_transpose(tmp_vec3_n2, flux); // LS^T tau S
        Fvec.add(darea_dx/area*JxW[qp], tmp_vec3_n2);
        
        
        if (request_jacobian && c.elem_solution_derivative)
        {
            A_sens.zero();
            
            for (unsigned int i_dim=0; i_dim<dim; i_dim++)
            {
                // Galerkin contribution from the advection flux terms
                B_mat.left_multiply(tmp_mat_n1n2, Ai_advection[i_dim]);
                dB_mat[i_dim].right_multiply_transpose(tmp_mat2_n2n2,
                                                       tmp_mat_n1n2);// dBw/dx_i^T  dF^adv_i/ dU
                Kmat.add(JxW[qp], tmp_mat2_n2n2);
                
                // discontinuity capturing term
                dB_mat[i_dim].right_multiply_transpose(tmp_mat2_n2n2,
                                                       dB_mat[i_dim]);
                Kmat.add(-JxW[qp]*diff_val(i_dim), tmp_mat2_n2n2);
                
                // sensitivity of Ai_Bi with respect to U:   [dAi/dUj.Bi.U  ...  dAi/dUn.Bi.U]
                dB_mat[i_dim].vector_mult(tmp_vec2_n1,
                                          c.get_elem_solution());
                for (unsigned int i_cvar=0; i_cvar<n1; i_cvar++)
                {
                    flux_jacobian_sens[i_dim][i_cvar].vector_mult(tmp_vec1_n1,
                                                                  tmp_vec2_n1);
                    for (unsigned int i_phi=0; i_phi<n_phi; i_phi++)
                        A_sens.add_column((n_phi*i_cvar)+i_phi,
                                          phi[i_phi][qp],
                                          tmp_vec1_n1); // assuming that all variables have same n_phi
                }
            }
            
            // Least square contribution of flux gradient
            LS_mat.get_transpose(tmp_mat3);
            tmp_mat3.right_multiply(Ai_Bi_advection); // LS^T tau d^2F^adv_i / dx dU   (Ai constant)
            Kmat.add(-JxW[qp], tmp_mat3);
            
            LS_mat.get_transpose(tmp_mat3);
            tmp_mat3.right_multiply(A_sens); // LS^T tau d^2F^adv_i / dx dU  (Ai sensitivity)
            Kmat.add(-JxW[qp], tmp_mat3);
            
            // contribution sensitivity of the LS.tau matrix
            Kmat.add(-JxW[qp], LS_sens);
            
            
            // Galerkin contribution from source term
            this->calculate_source_flux_jacobian(primitive_sol, tmp_mat_n1n2); // dS/dU
            B_mat.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat_n1n2);       // B^T dS/dU
            Kmat.add(darea_dx/area*JxW[qp], tmp_mat2_n2n2);
            
            // Least square contribution from source term
            LS_mat.get_transpose(tmp_mat3);
            tmp_mat3.right_multiply(tmp_mat_n1n2);       // LS^T dS/dU
            Kmat.add(darea_dx/area*JxW[qp], tmp_mat3);
            
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



void
MAST::ShockTubeFluidElem::calculate_source_flux(const PrimitiveSolution& sol,
                                                DenseVector<Real>& vec) {
    calculate_advection_flux(0, sol, vec);
    vec(1) -= sol.p;
    vec.scale(-1.);
}



void
MAST::ShockTubeFluidElem::calculate_source_flux_jacobian(const PrimitiveSolution& sol,
                                                         DenseMatrix<Real>& mat) {
    libMesh::DenseVector<Real> dpress_dp, vec;
    libMesh::DenseMatrix<Real> dcdp, dpdc;
    dpress_dp.resize(3); vec.resize(3);
    dcdp.resize(3,3); dpdc.resize(3,3);
    
    
    calculate_advection_flux_jacobian(0, sol, mat);
    
    // this calculates dp/dU =  dp/dV * dV/dU , where V and U are the
    // primitive and conservative variables
    calculate_conservative_variable_jacobian(sol, dcdp, dpdc);
    dpress_dp(0) = (sol.cp - sol.cv)*sol.T; // R T
    dpress_dp(2) = (sol.cp - sol.cv)*sol.rho; // R rho
    dpdc.vector_mult_transpose(vec, dpress_dp);
    for (unsigned int i=0; i<3; i++)
        mat(1,i) -= vec(i);
    mat.scale(-1.);
}

//#endif
