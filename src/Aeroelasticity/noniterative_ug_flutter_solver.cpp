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
#include "Aeroelasticity/noniterative_ug_flutter_solver.h"



void
MAST::NoniterativeUGFlutterSolver::scan_for_roots() {
    
    // if the initial scanning has not been done, then do it now
    if (!_flutter_solutions.size()) {
        // the outer loop consists of the reference velocity at which
        // aerodynamics is calculated
        Real current_v_ref = v_ref_range.second,
        delta_v_ref = (v_ref_range.second-v_ref_range.first)/n_v_ref_divs;
        
        std::vector<Real> v_ref_vals(n_v_ref_divs+1);
        for (unsigned int i=0; i<n_v_ref_divs+1; i++) {
            v_ref_vals[i] = current_v_ref;
            current_v_ref -= delta_v_ref;
        }
        v_ref_vals[n_v_ref_divs] = v_ref_range.first; // to get around finite-precision arithmetic

        for (unsigned int j=0; j<n_v_ref_divs+1; j++) {
            
            current_v_ref = v_ref_vals[j];
            
            // march from the upper limit to the lower to find the roots
            Real current_k_red = k_red_range.second,
            delta_k_red = (k_red_range.second-k_red_range.first)/n_k_red_divs;
            
            std::vector<Real> k_vals(n_k_red_divs+1);
            for (unsigned int i=0; i<n_k_red_divs+1; i++) {
                k_vals[i] = current_k_red;
                current_k_red -= delta_k_red;
            }
            k_vals[n_k_red_divs] = k_red_range.first; // to get around finite-precision arithmetic

            MAST::FlutterSolutionBase* prev_sol = NULL;
            
            for (unsigned int i=0; i<n_k_red_divs+1; i++) {
                current_k_red = k_vals[i];
                std::auto_ptr<MAST::FlutterSolutionBase> sol =
                analyze(current_k_red,
                        current_v_ref,
                        prev_sol);
                
                
                sol->print(_output);

                // add the solution to this solver
                _insert_new_solution(current_v_ref, sol.release());
                
                // now get a pointer to the previous solution
                // get the solution from the database for this reduced frequency
                std::map<Real, MAST::FlutterSolutionBase*>::iterator it =
                _flutter_solutions.find(current_k_red);
                
                libmesh_assert(it != _flutter_solutions.end());
                prev_sol = it->second;
            }
            
        }
        _identify_crossover_points();
    }
}



std::pair<bool, MAST::FlutterSolutionBase*>
MAST::NoniterativeUGFlutterSolver::bisection_search(const std::pair<MAST::FlutterSolutionBase*,
                                                    MAST::FlutterSolutionBase*>& ref_sol_range,
                                                    const unsigned int root_num,
                                                    const Real g_tol,
                                                    const unsigned int max_iters) {
    
    // assumes that the upper k_val has +ve g val and lower k_val has -ve
    // k_val
    Real
    lower_ref_val   = ref_sol_range.first->get_root(root_num).k_red_ref,
    lower_v_ref_val = ref_sol_range.first->get_root(root_num).V,
    lower_g         = ref_sol_range.first->get_root(root_num).g,
    upper_ref_val   = ref_sol_range.second->get_root(root_num).k_red_ref,
    upper_v_ref_val = ref_sol_range.second->get_root(root_num).V,
    upper_g         = ref_sol_range.second->get_root(root_num).g,
    new_k = 0.,
    new_v_ref = lower_v_ref_val +
    (upper_v_ref_val-lower_v_ref_val)/(upper_g-lower_g)*(0.-lower_g); // linear interpolation
    unsigned int n_iters = 0;
    
    std::auto_ptr<MAST::FlutterSolutionBase> new_sol;
    std::pair<bool, MAST::FlutterSolutionBase*> rval(false, NULL);
    
    while (n_iters < max_iters) {
        
        new_k = lower_ref_val +
        (upper_ref_val-lower_ref_val)/(upper_g-lower_g)*(0.-lower_g); // linear interpolation

        new_sol.reset(analyze(new_k,
                          new_v_ref,
                          ref_sol_range.first).release());
        
        new_sol->print(_output);

        // add the solution to this solver
        _insert_new_solution(new_v_ref, new_sol.release());

        // get the solution from the database for this reduced frequency
        std::map<Real, MAST::FlutterSolutionBase*>::iterator it =
        _flutter_solutions.find(new_k);
        
        libmesh_assert(it != _flutter_solutions.end());
        rval.second = it->second;
        const MAST::FlutterRootBase& root = rval.second->get_root(root_num);

        // use the estimated flutter velocity to get the next
        // V_ref value for aerodynamic matrices.
        new_v_ref = root.V;
        
        // check if the new damping value
        if (fabs(root.g) <= g_tol) {
            rval.first = true;
            return  rval;
        }
        
        // update the k_val
        if (root.g < 0.) {
            lower_ref_val = new_k;
            lower_g = root.g;
        }
        else {
            upper_ref_val = new_k;
            upper_g = root.g;
        }
        
        n_iters++;
    }
    
    // return false, along with the latest sol
    rval.first = false;
    
    return rval;
}



std::pair<bool, MAST::FlutterSolutionBase*>
MAST::NoniterativeUGFlutterSolver::newton_search(const MAST::FlutterSolutionBase& init_sol,
                                                 const unsigned int root_num,
                                                 const Real tol,
                                                 const unsigned int max_iters) {
    
    std::pair<bool, MAST::FlutterSolutionBase*> rval(false, NULL);
    // assumes that the upper k_val has +ve g val and lower k_val has -ve
    // k_val
    Real k_red, v_ref;
    unsigned int n_iters = 0;
    
    std::auto_ptr<MAST::FlutterSolutionBase> new_sol;

    RealVectorX res, sol, dsol;
    RealMatrixX jac;
    ComplexMatrixX mat_A, mat_B, mat_A_sens, mat_B_sens;
    ComplexVectorX v;

    res.resize(2); sol.resize(2); dsol.resize(2);
    jac.resize(2,2);
    
    // initialize the solution to the values of init_sol
    k_red  = init_sol.get_root(root_num).k_red_ref;
    v_ref  = init_sol.get_root(root_num).V_ref;
    sol(0) = k_red;
    sol(1) = v_ref;

    const MAST::FlutterSolutionBase* prev_sol = &init_sol;
    
    bool if_continue = true;
    
    while (if_continue) {

        // evaluate the residual and Jacobians
        std::auto_ptr<MAST::FlutterSolutionBase> ug_sol =
        this->analyze(k_red, v_ref, prev_sol);

        ug_sol->print(_output);

        // add the solution to this solver
        _insert_new_solution(v_ref, ug_sol.release());
        
        // now get a pointer to the previous solution
        // get the solution from the database for this reduced frequency
        std::map<Real, MAST::FlutterSolutionBase*>::iterator it =
        _flutter_solutions.find(k_red);
        
        libmesh_assert(it != _flutter_solutions.end());
        rval.second = it->second;
        prev_sol = it->second;

        // solve the Newton update problem
        const MAST::FlutterRootBase& root = prev_sol->get_root(root_num);
        Complex eig = root.root, eig_k_red_sens = 0., den = 0., eig_V_ref_sens = 0.;
        
        // initialize the baseline matrices
        initialize_matrices(k_red, v_ref, mat_A, mat_B);
        
        // solve the sensitivity problem
        // first with respect to k_red
        // next we need the sensitivity of k_red before we can calculate
        // the sensitivity of flutter eigenvalue
        initialize_matrix_sensitivity_for_reduced_freq(k_red,
                                                       v_ref,
                                                       mat_A_sens,
                                                       mat_B_sens);
        
        // now calculate the quotient for sensitivity wrt k_red
        // calculate numerator
        mat_B_sens *= -eig;
        mat_B_sens += mat_A_sens;
        v = mat_B_sens*root.eig_vec_right;
        den = root.eig_vec_left.dot(mat_B*root.eig_vec_right);
        eig_k_red_sens = root.eig_vec_left.dot(v) / den;
        
        /*std::auto_ptr<MAST::FlutterSolutionBase> ug_dsol =
        this->analyze(k_red+.001, v_ref, prev_sol);
        eig -= ug_dsol->get_root(root_num).root;
        eig /= -.001;
        
        std::cout << eig << std::endl;
        std::cout << eig_k_red_sens << std::endl;*/
        
        // next, sensitivity wrt V_ref
        initialize_matrix_sensitivity_for_V_ref(k_red,
                                                v_ref,
                                                mat_A_sens,
                                                mat_B_sens);

        // now calculate the quotient for sensitivity wrt V_ref
        // calculate numerator
        mat_B_sens *= -eig;
        mat_B_sens += mat_A_sens;
        v = mat_B_sens*root.eig_vec_right;
        den = root.eig_vec_left.dot(mat_B*root.eig_vec_right);
        eig_V_ref_sens = root.eig_vec_left.dot(v) / den;

        /*std::auto_ptr<MAST::FlutterSolutionBase> ug_dsol =
         this->analyze(k_red, v_ref+.001, prev_sol);
         eig -= ug_dsol->get_root(root_num).root;
         eig /= -.001;
         
         std::cout << eig << std::endl;
         std::cout << eig_V_ref_sens << std::endl;*/

        // residual
        res(0) = eig.imag()/eig.real();
        res(1) = v_ref*std::sqrt(eig.real()) - 1.;
        
        // Jacobian
        jac(0,0) =
        eig_k_red_sens.imag()/eig.real() -
        eig.imag()*pow(eig.real(),-2)*eig_k_red_sens.real();
        jac(0,1) = eig_V_ref_sens.imag()/eig.real() -
        eig.imag()*pow(eig.real(),-2)*eig_V_ref_sens.real();
        jac(1,0) = 0.5*v_ref*pow(eig.real(), -0.5)*eig_k_red_sens.real();
        jac(1,1) = std::sqrt(eig.real()) + 0.5*v_ref*pow(eig.real(), -0.5)*eig_V_ref_sens.real();

        // now calculate the updates
        //     r0 + J *dx = 0
        // =>  dx = - inv(J) * r0
        // =>  x1 = x0 + dx
        
        dsol = -jac.inverse()*res;
        sol += dsol;
        
        // get the updated parameter values
        k_red = sol(0);
        v_ref = sol(1);
        
        // increment the iteration counter
        n_iters++;
        
        // set the flag
        if (fabs(root.g) < tol) {
            rval.first = true;
            return rval;
        }
        
        if (n_iters >= max_iters)
            if_continue = false;
    }
    
    // return false, along with the latest sol
    rval.first = false;
    
    return rval;
}



void
MAST::NoniterativeUGFlutterSolver::_insert_new_solution(const Real v_ref,
                                                        MAST::FlutterSolutionBase* sol) {
    
    const Real k_red = sol->ref_val();
    
    // if the value does not already exist, then insert it in the map
    std::map<Real, MAST::FlutterSolutionBase*>::iterator it =
    _flutter_solutions.find(k_red);
    
    if ( it == _flutter_solutions.end()) {
        bool if_success =
        _flutter_solutions.insert(std::pair<Real, MAST::FlutterSolutionBase*>
                                  (k_red, sol)).second;
        
        libmesh_assert(if_success);
        
        // do not delete the solution since it is stored in the solver
    }
    else {
        // assuming that the roots are sorted, compare them with the existing
        // root and swap the ones that are closer to v_ref

        // first make sure that the two have the same number of roots
        libmesh_assert_equal_to(it->second->n_roots(), sol->n_roots());

        Real dV_old = 0., v_ref_old = 0., dV_new = 0.;
        MAST::FlutterRootBase *old_root, *new_root;
        
        for (unsigned int i=0; i<sol->n_roots(); i++) {
            old_root  = &(it->second->get_root(i));
            new_root  = &(sol->get_root(i));
            v_ref_old = old_root->V_ref;
            dV_old    = fabs(old_root->V - v_ref_old);
            dV_new    = fabs(new_root->V - v_ref);
            
            if (dV_new < dV_old)
                old_root->copy_root(*new_root);
        }
        
        // delete the solution since its relevants roots have been swapped
        // out
        delete sol;
    }
}



void
MAST::NoniterativeUGFlutterSolver::calculate_sensitivity(MAST::FlutterRootBase& root,
                                                         const libMesh::ParameterVector& params,
                                                         const unsigned int i) {
    // make sure that the aero_structural_model is a valid pointer
    libmesh_assert(aero_structural_model);
    
    libMesh::out
    << " ====================================================" << std::endl
    << "UG Sensitivity Solution" << std::endl
    << "   k_red = " << std::setw(10) << root.k_red << std::endl
    << "   V_ref = " << std::setw(10) << root.V << std::endl;
    
    Complex eig = root.root, sens = 0., k_sens = 0., vref_sens = 0., den = 0.;
    
    // matrix to store the sensitivity of constraint wrt kref and vref
    // first row is constraint  f1(k,v) = g = im(lambda)/re(lambda) = 0
    // second row is constraint f2(k,v) = vref*sqrt(re(lambda))-1.  = 0
    RealMatrixX jac;
    
    // residual vector for the sensitivity problem, and total derivative vec
    // first row is for constraint f1(k,v)
    // first row is for constraint f2(k,v)
    RealVectorX res, dsens;
    jac.setZero(2,2);
    res.setZero(2);
    
    // get the sensitivity of the matrices
    ComplexMatrixX mat_A, mat_B, mat_A_sens, mat_B_sens;
    ComplexVectorX v;
    
    // initialize the baseline matrices
    initialize_matrices(root.k_red_ref, root.V_ref, mat_A, mat_B);
    
    // calculate the eigenproblem sensitivity
    initialize_matrix_sensitivity_for_param(params, i,
                                            root.k_red_ref,
                                            root.V_ref,
                                            mat_A_sens,
                                            mat_B_sens);
    
    // the eigenproblem is     A x - lambda B x = 0
    // therefore, the denominator is obtained from the inner product of
    // x^T B x
    // sensitivity is
    //   -dlambda/dp x^T B x = - x^T (dA/dp - lambda dB/dp)
    // or
    //   dlambda/dp = [x^T (dA/dp - lambda dB/dp)]/(x^T B x)
    
    // now calculate the quotient for sensitivity
    // numerator =  ( dA/dp - lambda dB/dp)
    mat_B_sens *= -eig;
    mat_B_sens += mat_A_sens;
    v = mat_B_sens*root.eig_vec_right;
    den = root.eig_vec_left.dot(mat_B*root.eig_vec_right);
    sens = root.eig_vec_left.dot(v)/den;
    
    // sensitivity of f1(k,v) = im(lambda)/re(lambda) = 0
    res(0) =
    sens.imag()/eig.real() - eig.imag()/pow(eig.real(),2) * sens.real();
    // sensitivity of f2(k,v) = vref*sqrt(real(lambda)) - 1 = 0
    res(1) = root.V_ref*0.5*pow(eig.real(),-0.5)*sens.real();
    
    
    // next we need the sensitivity of k_red before we can calculate
    // the sensitivity of flutter eigenvalue
    initialize_matrix_sensitivity_for_reduced_freq(root.k_red_ref,
                                                   root.V_ref,
                                                   mat_A_sens,
                                                   mat_B_sens);
    
    // now calculate the quotient for sensitivity wrt k_red
    // calculate numerator
    mat_B_sens *= -eig;
    mat_B_sens += mat_A_sens;
    v = mat_B_sens*root.eig_vec_right;
    k_sens = root.eig_vec_left.dot(v) / den;
    
    // use this to calculate the partial derivative of f1 wrt k_red
    jac(0,0) =
    k_sens.imag()/eig.real() - eig.imag()/pow(eig.real(),2)*k_sens.real();
    jac(1,0) = root.V_ref*0.5*pow(eig.real(),-0.5)*k_sens.real();
    
    
    // next we need the sensitivity of Vref constraint before we can calculate
    // the sensitivity of flutter eigenvalue
    initialize_matrix_sensitivity_for_V_ref(root.k_red_ref,
                                            root.V_ref,
                                            mat_A_sens,
                                            mat_B_sens);
    
    // now calculate the quotient for sensitivity wrt k_red
    // calculate numerator
    mat_B_sens *= -eig;
    mat_B_sens += mat_A_sens;
    v = mat_B_sens*root.eig_vec_right;
    vref_sens = root.eig_vec_left.dot(v) / den;
    
    // use this to calculate the partial derivative of f2 wrt vref
    jac(0,1) =
    vref_sens.imag()/eig.real() - eig.imag()/pow(eig.real(),2)*vref_sens.real();
    jac(1,1) =
    sqrt(eig.real()) + root.V_ref*0.5*pow(eig.real(),-0.5)*vref_sens.real();

    // now invert the Jacobian to calculate the sensitivity
    dsens = -jac.inverse()*res;
    
    // finally add the correction to the flutter sensitivity
    sens += k_sens * dsens(0) + vref_sens * dsens(1);
    
    // set value in the return root
    root.has_sensitivity_data = true;
    root.root_sens  = sens;
    root.V_sens     = -.5*sens.real()/pow(eig.real(), 1.5);
    
    libMesh::out
    << "Finished UG Sensitivity Solution" << std::endl
    << " ====================================================" << std::endl;
    
}



