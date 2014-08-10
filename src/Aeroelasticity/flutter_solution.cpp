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
#include "Aeroelasticity/flutter_solution.h"
#include "Aeroelasticity/flutter_solver_base.h"




void
MAST::FlutterSolutionBase::sort(const MAST::FlutterSolutionBase& sol)
{
    const unsigned int nvals = this->n_roots();
    libmesh_assert_equal_to(nvals, sol.n_roots());
    
    // two roots with highest modal_participation dot product are sorted
    // in the same serial order
    for (unsigned int i=0; i<nvals-1; i++)
    {
        const MAST::FlutterRootBase& r = sol.get_root(i);
        Real max_val = 0.;
        Complex val = 0.;
        unsigned int max_val_root = nvals+1;
        for (unsigned int j=i; j<nvals; j++) {
            val = r.eig_vec_left.dot(_Bmat*_roots[j]->eig_vec_right);
            //_roots[j]->modal_participation.dot(r.modal_participation);
            // scale by the eigenvalue separation with the assumption that
            // the roots will be closer to each other than any other
            // root at two consecutive eigenvalues. In other words,
            // we are penalizing the dot product with the eigenvalue
            // distance
            val /= abs(r.root-_roots[j]->root);
            if (abs(val) > max_val) {
                max_val = abs(val);
                max_val_root = j;
            }
        }
        
        // now we should have the one with highest dot product
        std::swap(_roots[i], _roots[max_val_root]);
    }
}



void
MAST::FlutterSolutionBase::swap_root(MAST::FlutterSolutionBase& sol,
                                     unsigned int root_num) {
    libmesh_assert(root_num < _roots.size());
    libmesh_assert(root_num < sol._roots.size());
    
    std::swap(_roots[root_num], sol._roots[root_num]);
}



void
MAST::FlutterSolutionBase::mark_inconsistent_roots_as_nonphysical(const Real ref_val,
                                                                  const Real tol) {
    // iterate over the roots and identify the roots where the reference value
    // differs from the calculated value by more than the tolerance
    
    for (unsigned int i=0; i<_roots.size(); i++)
        if (!_roots[i]->if_nonphysical_root &&
            fabs(_roots[i]->V - ref_val) > tol)
            _roots[i]->if_nonphysical_root = true;
}







void
MAST::FlutterSolutionBase::print(std::ostream &output)
{
    // make sure that some roots exist for this solution
    libmesh_assert(this->n_roots() > 0);
    
    const unsigned int nvals = this->n_roots();
    unsigned int n_participation_vals =
    (unsigned int)this->get_root(0).modal_participation.size();
    libmesh_assert(nvals);

    // first write the reference values of the root
    output << " Flutter Root " << std::endl;

    // now write the root
    output
    << std::setw(5) << "#"
    << std::setw(15) << "k_ref"
    << std::setw(15) << "V_ref"
    << std::setw(15) << "g"
    << std::setw(15) << "V"
    << std::setw(15) << "k"
    << std::setw(15) << "omega"
    << std::setw(15) << "Re"
    << std::setw(15) << "Im"
    << std::setw(3)  << " | ";
    
    // output the headers for flutter mode participation
    for (unsigned int i=0; i<n_participation_vals; i++)
        output
        << std::setw(2) << " "
        << std::setw(5) << "Mode "
        << std::setw(5) << i
        << std::setw(2) << " ";
    
    // output the headers for flutter mode
    for (unsigned int i=0; i<nvals; i++)
        output
        << std::setw(10) << "|         "
        << std::setw(5) << "Mode "
        << std::setw(5) << i
        << std::setw(10) << "         |";
    output << std::endl;
    
    for (unsigned int i=0; i<nvals; i++)
    {
        const MAST::FlutterRootBase& root = this->get_root(i);
        
        if (root.if_nonphysical_root)
        {
            std::stringstream oss;
            oss << "**" << i;
            output << std::setw(5) << oss.str();
        }
        else
            output << std::setw(5) << i;
        
        // flutter root details
        output
        << std::setw(15) << root.k_red_ref
        << std::setw(15) << root.V_ref
        << std::setw(15) << root.g
        << std::setw(15) << root.V
        << std::setw(15) << root.k_red
        << std::setw(15) << root.omega
        << std::setw(15) << std::real(root.root)
        << std::setw(15) << std::imag(root.root)
        << std::setw(3)  << " | ";
        
        // now write the modal participation
        for (unsigned int j=0; j<n_participation_vals; j++)
            output
            << std::setw(12) << root.modal_participation(j)
            << std::setw(2) << " ";
        
        // now write the flutter mode
        for (unsigned int j=0; j<nvals; j++)
        {
            output
            << std::setw(2) << "| "
            << std::setw(12) << std::real(root.eig_vec_right(j))
            << std::setw(2) << " "
            << std::setw(12) << std::imag(root.eig_vec_right(j))
            << std::setw(2) << " |";
        }
        output << std::endl;
    }
}



void
MAST::FlutterRootCrossoverBase::print(std::ostream &output) const
{
    const MAST::FlutterRootBase
    &lower = crossover_solutions.first->get_root(root_num),
    &upper = crossover_solutions.second->get_root(root_num);
    
    output
    << " Lower Root: " << std::endl
    << "    k : " << std::setw(15) << lower.k_red << std::endl
    << "    g : " << std::setw(15) << lower.g << std::endl
    << "    V : " << std::setw(15) << lower.V << std::endl
    << "omega : " << std::setw(15) << lower.omega << std::endl
    << " Upper Root: " << std::endl
    << "    k : " << std::setw(15) << upper.k_red << std::endl
    << "    g : " << std::setw(15) << upper.g << std::endl
    << "    V : " << std::setw(15) << upper.V << std::endl
    << "omega : " << std::setw(15) << upper.omega << std::endl;
    
    if (root)
        output
        << "Critical Root: " << std::endl
        << "    k : " << std::setw(15) << root->k_red << std::endl
        << "    g : " << std::setw(15) << root->g << std::endl
        << "    V : " << std::setw(15) << root->V << std::endl
        << "omega : " << std::setw(15) << root->omega << std::endl;
    else
        output << "Critical root not yet calculated." << std::endl;
}
