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
                
                
                sol->print(_output, _mode_output);

                // add the solution to this solver
                bool if_delete = _insert_new_solution(current_v_ref, *sol);
                
                // if this was not a new reduced frequency, then this sol
                // contains the lower quality approximations with respect to
                // v_ref, so delete it
                if (if_delete)
                    sol.reset();
                else
                    sol.release();

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
        
        new_sol->print(_output, _mode_output);

        // add the solution to this solver
        bool if_delete = _insert_new_solution(new_v_ref, *new_sol);

        // if this was not a new reduced frequency, then this sol contains
        // the lower quality approximations with respect to v_ref, so delete
        // it
        if (if_delete)
            new_sol.reset();
        else
            new_sol.release();
        
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


bool
MAST::NoniterativeUGFlutterSolver::_insert_new_solution(const Real v_ref,
                                                        MAST::FlutterSolutionBase& sol) {
    
    const Real k_red = sol.ref_val();
    
    // if the value does not already exist, then insert it in the map
    std::map<Real, MAST::FlutterSolutionBase*>::iterator it =
    _flutter_solutions.find(k_red);
    
    if ( it == _flutter_solutions.end()) {
        bool if_success =
        _flutter_solutions.insert(std::pair<Real, MAST::FlutterSolutionBase*>
                                  (k_red, &sol)).second;
        
        libmesh_assert(if_success);
        
        // return false so that the calling function knows not to delete it
        return false;
    }
    else {
        // assuming that the roots are sorted, compare them with the existing
        // root and swap the ones that are closer to v_ref

        // first make sure that the two have the same number of roots
        libmesh_assert_equal_to(it->second->n_roots(), sol.n_roots());

        Real dV_old = 0., v_ref_old = 0., dV_new = 0.;
        const MAST::FlutterRootBase *old_root, *new_root;
        
        for (unsigned int i=0; i<sol.n_roots(); i++) {
            old_root  = &(it->second->get_root(i));
            new_root  = &(sol.get_root(i));
            v_ref_old = old_root->V_ref;
            dV_old    = fabs(it->second->get_root(i).V - v_ref_old);
            dV_new    = fabs(sol.get_root(i).V - v_ref);
            
            if (dV_new < dV_old)
                it->second->swap_root(sol, i);
        }
        
        // tell the calling function to delete the sol since it is no longer
        // of use
        return true;
    }
}

