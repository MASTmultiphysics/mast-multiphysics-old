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
#include "Aeroelasticity/flutter_solver_base.h"



MAST::FlutterSolverBase::~FlutterSolverBase() {
    this->clear_solutions();
}




void
MAST::FlutterSolverBase::clear_solutions() {

    std::map<Real, MAST::FlutterSolutionBase*>::iterator it =
    _flutter_solutions.begin();
    
    for ( ; it != _flutter_solutions.end(); it++)
        delete it->second;
    
    std::map<Real, MAST::FlutterRootCrossoverBase*>::iterator cross_it =
    _flutter_crossovers.begin();
    
    for ( ; cross_it != _flutter_crossovers.end(); cross_it++)
        delete cross_it->second;
    
    _flutter_solutions.clear();
    _flutter_crossovers.clear();
}


unsigned int
MAST::FlutterSolverBase::n_roots_found() const {
    
    std::map<Real, MAST::FlutterRootCrossoverBase*>::const_iterator
    it = _flutter_crossovers.begin(),
    end = _flutter_crossovers.end();
    
    unsigned int n = 0;
    for ( ; it!=end; it++)
        if (it->second->root) // a valid root pointer has been assigned
            n++;
    
    return n;
}



const MAST::FlutterRootBase&
MAST::FlutterSolverBase::get_root(const unsigned int n) const {
    
    libmesh_assert(n < n_roots_found());
    
    std::map<Real, MAST::FlutterRootCrossoverBase*>::const_iterator
    it = _flutter_crossovers.begin(),
    end = _flutter_crossovers.end();
    
    unsigned int root_num = 0;
    for ( ; it!=end; it++) {
        if (it->second->root) { // a valid root pointer has been assigned
            if (root_num == n)  // root num matches the one being requested
                return *(it->second->root);
            else
                root_num++;
        }
    }
    
    libmesh_assert(false); // should not get here
}




void MAST::FlutterSolverBase::scan_for_roots()
{
    // if the initial scanning has not been done, then do it now
    if (!_flutter_solutions.size()) {
        // march from the upper limit to the lower to find the roots
        Real current_k_red = ref_val_range.second,
        delta_k_red = (ref_val_range.second-ref_val_range.first)/n_ref_val_divs;
        
        std::vector<Real> k_vals(n_ref_val_divs+1);
        for (unsigned int i=0; i<n_ref_val_divs+1; i++) {
            k_vals[i] = current_k_red;
            current_k_red -= delta_k_red;
        }
        k_vals[n_ref_val_divs] = ref_val_range.first; // to get around finite-precision arithmetic
        
        MAST::FlutterSolutionBase* prev_sol = NULL;
        for (unsigned int i=0; i<n_ref_val_divs+1; i++) {
            current_k_red = k_vals[i];
            prev_sol = analyze(current_k_red, prev_sol);
        }
        
        _identify_crossover_points();
    }
}




std::pair<bool, const MAST::FlutterRootBase*>
MAST::FlutterSolverBase::find_next_root(const Real g_tol,
                                        const unsigned int n_bisection_iters)
{
    // iterate over the cross-over points and calculate the next that has
    // not been evaluated
    std::map<Real, MAST::FlutterRootCrossoverBase*>::iterator
    it = _flutter_crossovers.begin(), end = _flutter_crossovers.end();
    while ( it != end)
    {
        MAST::FlutterRootCrossoverBase* cross = it->second;
        
        if (!cross->root) {
            const unsigned int root_num = cross->root_num;
            std::pair<bool, MAST::FlutterSolutionBase*> bisection_sol =
            bisection_search(cross->crossover_solutions,
                             root_num, g_tol, n_bisection_iters);
            cross->root = &(bisection_sol.second->get_root(root_num));
            
            // now, remove this entry from the _flutter_crossover points and
            // reinsert it with the actual critical velocity
            _flutter_crossovers.erase(it);
            std::pair<Real, MAST::FlutterRootCrossoverBase*>
            val(cross->root->V, cross);
            _flutter_crossovers.insert(val);
            return std::pair<bool, const MAST::FlutterRootBase*> (true, cross->root);
        }
        
        it++;
    }
    
    // if it gets here, no new root was found
    return std::pair<bool, MAST::FlutterRootBase*> (false, NULL);
}



std::pair<bool, const MAST::FlutterRootBase*>
MAST::FlutterSolverBase::find_critical_root(const Real g_tol,
                                            const unsigned int n_bisection_iters)
{
    // iterate over the cross-over points and calculate the next that has
    // not been evaluated
    std::map<Real, MAST::FlutterRootCrossoverBase*>::iterator
    it = _flutter_crossovers.begin(), end = _flutter_crossovers.end();
    
    if (it == end) // no potential cross-over points were identified
        return std::pair<bool, MAST::FlutterRootBase*> (false, NULL);
    
    // it is possible that once the root has been found, its velocity end up
    // putting it at a higher velocity in the map, so we need to check if
    // the critical root has changed
    while (!it->second->root)
    {
        MAST::FlutterRootCrossoverBase* cross = it->second;
        
        if (!cross->root) {
            const unsigned int root_num = cross->root_num;
            std::pair<bool, MAST::FlutterSolutionBase*> bisection_sol =
            bisection_search(cross->crossover_solutions,
                             root_num, g_tol, n_bisection_iters);
            cross->root = &(bisection_sol.second->get_root(root_num));
            
            // now, remove this entry from the _flutter_crossover points and
            // reinsert it with the actual critical velocity
            _flutter_crossovers.erase(it);
            std::pair<Real, MAST::FlutterRootCrossoverBase*>
            val(cross->root->V, cross);
            _flutter_crossovers.insert(val);
        }
        
        // update the iterator to make sure that this is updated
        it = _flutter_crossovers.begin(); end = _flutter_crossovers.end();
    }
    
    // if it gets here, then the root was successfully found
    return std::pair<bool, const MAST::FlutterRootBase*> (true, it->second->root);
}





std::pair<bool, MAST::FlutterSolutionBase*>
MAST::FlutterSolverBase::bisection_search(const std::pair<MAST::FlutterSolutionBase*,
                                          MAST::FlutterSolutionBase*>& ref_sol_range,
                                          const unsigned int root_num,
                                          const Real g_tol,
                                          const unsigned int max_iters)
{
    // assumes that the upper k_val has +ve g val and lower k_val has -ve
    // k_val
    Real lower_ref_val = ref_sol_range.first->ref_val(),
    lower_g = ref_sol_range.first->get_root(root_num).g,
    upper_ref_val = ref_sol_range.second->ref_val(),
    upper_g = ref_sol_range.second->get_root(root_num).g,
    new_k = 0.;
    unsigned int n_iters = 0;
    
    MAST::FlutterSolutionBase* new_sol = NULL;
    
    while (n_iters < max_iters) {
        
        new_k = lower_ref_val +
        (upper_ref_val-lower_ref_val)/(upper_g-lower_g)*(0.-lower_g); // linear interpolation
        new_sol = analyze(new_k, ref_sol_range.first);
        const MAST::FlutterRootBase& root = new_sol->get_root(root_num);
        
        // check if the new damping value
        if (fabs(root.g) <= g_tol)
            return std::pair<bool, MAST::FlutterSolutionBase*>
            (true, new_sol);
        
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
    
    return std::pair<bool, MAST::FlutterSolutionBase*>
    (false, new_sol); // return false, along with the latest sol
}



void
MAST::FlutterSolverBase::print_sorted_roots(std::ostream* output)
{
    if (!output)
        output = &_output;
    
    std::map<Real, MAST::FlutterSolutionBase*>::const_iterator
    sol_it = _flutter_solutions.begin(),
    sol_end = _flutter_solutions.end();
    libmesh_assert(sol_it != sol_end); // solutions should have been evaluated
    
    unsigned int nvals = sol_it->second->n_roots();
    libmesh_assert(nvals); // should not be zero
    
    // each root is written separately one after another
    for (unsigned int i=0; i<nvals; i++)
    {
        // print the headers
        *output
        << "** Root # "
        << std::setw(5) << i << " **" << std::endl
        << std::setw(15) << "k"
        << std::setw(15) << "Re"
        << std::setw(15) << "Im"
        << std::setw(15) << "g"
        << std::setw(15) << "V"
        << std::setw(15) << "omega" << std::endl;
        
        // update the iterator for this analysis
        sol_it = _flutter_solutions.begin();
        
        // write the data from all solutions
        for ( ; sol_it != sol_end; sol_it++)
        {
            const MAST::FlutterRootBase& root =
            sol_it->second->get_root(i);
            
            *output
            << std::setw(15) << root.k_red
            << std::setw(15) << std::real(root.root)
            << std::setw(15) << std::imag(root.root)
            << std::setw(15) << root.g;
            
            // check if the root might be the odd artifact of the UG method
            if (!root.if_nonphysical_root)
                *output << std::setw(15) << root.V;
            else
                *output << "** " << std::setw(12) << root.V;
            *output
            << std::setw(15) << root.omega << std::endl;
        }
        *output << std::endl << std::endl;
    }
    
    
    // write the roots identified using iterative search technique
    std::streamsize prec = output->precision();
    
    unsigned int nroots = this->n_roots_found();
    *output << std::endl
    << "n critical roots identified: " << nroots << std::endl;
    for (unsigned int i=0; i<nroots; i++)
    {
        const MAST::FlutterRootBase& root = this->get_root(i);
        *output
        << "** Root : " << std::setw(5) << i << " **" << std::endl
        << "g      = " << std::setw(15) << root.g << std::endl
        << "V      = " << std::setw(35) << std::setprecision(15) << root.V << std::endl
        << "omega  = " << std::setw(35) << std::setprecision(15) << root.omega << std::endl
        << std::setprecision(prec) // set the precision to the default value
        << "k_red  = " << std::setw(15) << root.k_red << std::endl
        << "Modal Participation : " << std::endl ;
        for (unsigned int j=0; j<nvals; j++)
            *output
            << "(" << std::setw(5) << j << "): "
            << std::setw(10) << root.modal_participation(j)
            << std::setw(3)  << " ";
        *output << std::endl << std::endl;
    }
    
    
}


void
MAST::FlutterSolverBase::print_crossover_points(std::ostream* output)
{
    if (!output)
        output = &_output;
    
    *output << "n crossover points found: "
    << std::setw(5) << _flutter_crossovers.size() << std::endl;
    
    std::multimap<Real, MAST::FlutterRootCrossoverBase*>::const_iterator
    it = _flutter_crossovers.begin(), end = _flutter_crossovers.end();
    
    unsigned int i=0;
    
    for ( ; it != end; it++) {
        *output << "** libMesh::Point : " << std::setw(5) << i << " **" << std::endl;
        it->second->print(*output);
        *output << std::endl;
        i++;
    }
}



