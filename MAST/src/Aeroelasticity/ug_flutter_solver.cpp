//
//  ug_flutter_solver.cpp
//  MAST
//
//  Created by Manav Bhatia on 9/16/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

// C++ includes
#include <iomanip>

// MAST includes
#include "Aeroelasticity//ug_flutter_solver.h"


void
UGFlutterSolver::UGFlutterRoot::init(const Real k, const Real b_ref,
                                     const Complex num,
                                     const Complex den,
                                     const ComplexMatrixX& Bmat,
                                     const ComplexVectorX& eig_vec)
{
    k_ref = k;

    if (std::abs(den) > 0.)
    {
        root = num/den;
        if (std::real(root) > 0.)
        {
            V     = sqrt(1./std::real(root));
            g     = std::imag(root)/std::real(root);
            omega = k_ref*V/b_ref;
            if_imag_V = false;
        }
        else
        {
            V     = 0.;
            g     = 0.;
            omega = 0.;
            if_imag_V = true;
        }
    }

    // calculate the modal participation vector
    const unsigned int nvals = (int)Bmat.rows();
    mode = eig_vec;
    ComplexVectorX k_q = Bmat * mode;
    modal_participation.resize(nvals, 1);
    for (unsigned int i=0; i<nvals; i++)
        modal_participation(i) =  std::abs(std::conj(mode(i)) * k_q(i));
    modal_participation *= (1./modal_participation.sum());
}




void
UGFlutterSolver::UGFlutterSolution::init (const Real kref, const Real bref,
                                          const LAPACK_ZGGEV& eig_sol)
{
    _k_ref = kref;
    // iterate over the roots and initialize the vector
    const ComplexMatrixX &B = eig_sol.B(), &VR = eig_sol.right_eigenvectors();
    const ComplexVectorX &num = eig_sol.alphas(), &den = eig_sol.betas();
    unsigned int nvals = (int)B.rows();
    
    for (unsigned int i=0; i<nvals; i++)
        _roots[i]->init(kref, bref, num(i), den(i), B, VR.col(i));
}




void
UGFlutterSolver::UGFlutterSolution::sort(const UGFlutterSolver::UGFlutterSolution& sol)
{
    const unsigned int nvals = this->n_roots();
    libmesh_assert_equal_to(nvals, sol.n_roots());
    
    // two roots with highest modal_participation dot product are sorted
    // in the same serial order
    for (unsigned int i=0; i<nvals-1; i++)
    {
        const UGFlutterSolver::UGFlutterRoot& r = sol.get_root(i);
        Real max_val = 0., val = 0.;
        unsigned int max_val_root = nvals+1;
        for (unsigned int j=i; j<nvals; j++) {
            val = _roots[j]->modal_participation.dot(r.modal_participation);
            if (val > max_val) {
                max_val = val;
                max_val_root = j;
            }
        }
        
        // now we should have the one with highest dot product
        std::swap(_roots[i], _roots[max_val_root]);
    }
}



void
UGFlutterSolver::UGFlutterSolution::print(std::ostream &output,
                                          std::ostream &mode_output)
{
    const unsigned int nvals = this->n_roots();
    libmesh_assert(nvals);
    
    output << "k = " << _k_ref << std::endl
    << std::setw(5) << "UG #"
    << std::setw(15) << "Re"
    << std::setw(15) << "Im"
    << std::setw(15) << "g"
    << std::setw(15) << "V"
    << std::setw(15) << "omega";
    
    // output the headers for flutter mode participation
    for (unsigned int i=0; i<nvals; i++)
        output
        << std::setw(2) << " "
        << std::setw(5) << "Mode "
        << std::setw(5) << i
        << std::setw(3) << " ";
    output << std::endl;
    
    // output the headers for flutter mode
    mode_output << "k = " << _k_ref << "  (Right Eigenvectors)  " << std::endl;
    mode_output << std::setw(5) << "UG #";
    for (unsigned int i=0; i<nvals; i++)
        mode_output
        << std::setw(10) << " "
        << std::setw(5) << "Mode "
        << std::setw(5) << i
        << std::setw(10) << " ";
    mode_output << std::endl;
    
    for (unsigned int i=0; i<nvals; i++)
    {
        const UGFlutterSolver::UGFlutterRoot& root = this->get_root(i);
        
        if (root.if_imag_V)
        {
            std::stringstream oss;
            oss << "**" << i;
            output << std::setw(5) << oss.str();
            mode_output << std::setw(5) << oss.str();
        }
        else
        {
            output << std::setw(5) << i;
            mode_output << std::setw(5) << i;
        }
        
        // flutter root details
        output
        << std::setw(15) << std::real(root.root)
        << std::setw(15) << std::imag(root.root)
        << std::setw(15) << root.g
        << std::setw(15) << root.V
        << std::setw(15) << root.omega;
        
        // now write the modal participation
        for (unsigned int j=0; j<nvals; j++)
            output << std::setw(15) << root.modal_participation(j);
        output << std::endl;
        
        // now write the flutter mode
        for (unsigned int j=0; j<nvals; j++)
        {
            mode_output
            << std::setw(13) << std::real(root.mode(j))
            << std::setw(4) << " "
            << std::setw(13) << std::imag(root.mode(j));
        }
        mode_output << std::endl;
    }
}



void
UGFlutterSolver::UGFlutterRootCrossover::print(std::ostream &output) const
{
    const UGFlutterSolver::UGFlutterRoot
    &lower = crossover_solutions.first->get_root(root_num),
    &upper = crossover_solutions.second->get_root(root_num);
    
    output
    << " Lower Root: " << std::endl
    << "    k : " << std::setw(15) << lower.k_ref << std::endl
    << "    g : " << std::setw(15) << lower.g << std::endl
    << "    V : " << std::setw(15) << lower.V << std::endl
    << "omega : " << std::setw(15) << lower.omega << std::endl
    << " Upper Root: " << std::endl
    << "    k : " << std::setw(15) << upper.k_ref << std::endl
    << "    g : " << std::setw(15) << upper.g << std::endl
    << "    V : " << std::setw(15) << upper.V << std::endl
    << "omega : " << std::setw(15) << upper.omega << std::endl;
    
    if (root)
        output
        << "Critical Root: " << std::endl
        << "    k : " << std::setw(15) << root->k_ref << std::endl
        << "    g : " << std::setw(15) << root->g << std::endl
        << "    V : " << std::setw(15) << root->V << std::endl
        << "omega : " << std::setw(15) << root->omega << std::endl;
    else
        output << "Critical root not yet calculated." << std::endl;
}




UGFlutterSolver::~UGFlutterSolver()
{
    std::map<Real, UGFlutterSolver::UGFlutterSolution*>::iterator it =
    _flutter_solutions.begin();
    
    for ( ; it != _flutter_solutions.end(); it++)
        delete it->second;

    std::map<Real, UGFlutterSolver::UGFlutterRootCrossover*>::iterator cross_it =
    _flutter_crossovers.begin();
    
    for ( ; cross_it != _flutter_crossovers.end(); cross_it++)
        delete cross_it->second;
}



unsigned int
UGFlutterSolver::n_roots_found() const
{
    std::map<Real, UGFlutterSolver::UGFlutterRootCrossover*>::const_iterator
    it = _flutter_crossovers.begin(),
    end = _flutter_crossovers.end();
    
    unsigned int n = 0;
    for ( ; it!=end; it++)
        if (it->second->root) // a valid root pointer has been assigned
            n++;
    
    return n;
}



const FlutterRoot&
UGFlutterSolver::get_root(const unsigned int n) const
{
    libmesh_assert(n < n_roots_found());

    std::map<Real, UGFlutterSolver::UGFlutterRootCrossover*>::const_iterator
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





void UGFlutterSolver::scan_for_roots()
{
    // if the initial scanning has not been done, then do it now
    if (!_flutter_solutions.size()) {
        // march from the upper limit to the lower to find the roots
        Real current_k_ref = k_ref_range.second,
        delta_k_ref = (k_ref_range.second-k_ref_range.first)/n_k_divs;
        
        std::vector<Real> k_vals(n_k_divs+1);
        for (unsigned int i=0; i<n_k_divs+1; i++)
        {
            k_vals[i] = current_k_ref;
            current_k_ref -= delta_k_ref;
        }
        k_vals[n_k_divs] = k_ref_range.first; // to get around finite-precision arithmetic
        
        UGFlutterSolver::UGFlutterSolution* prev_sol = NULL;
        for (unsigned int i=0; i<n_k_divs+1; i++)
        {
            current_k_ref = k_vals[i];
            prev_sol = analyze(current_k_ref, prev_sol);
        }
        
        const Real tol = 1.0e-5, max_allowable_g = 0.75;
        
        //
        // for some cases some roots trail along the g=0 axis
        // and should not be considered as flutter. These are simply
        // modes where aerodynamics do not provide any damping.
        // These modes will not have a damping more than tolerance
        //
        const unsigned int nvals = prev_sol->n_roots();
        std::vector<bool> modes_to_neglect(nvals);
        std::fill(modes_to_neglect.begin(),
                  modes_to_neglect.end(), false);

        // look for the max g val for a mode, which will indicate if the
        // mode is undamped or not
        for (unsigned int i=0; i<nvals; i++) {
            std::map<Real, UGFlutterSolver::UGFlutterSolution*>::const_iterator
            sol_it    = _flutter_solutions.begin(),
            sol_end   = _flutter_solutions.end();
            Real max_g_val = 0., val = 0.;
            for ( ; sol_it!=sol_end; sol_it++) {
                val = fabs(sol_it->second->get_root(i).g);
                if (val > max_g_val)
                    max_g_val = val;
            }
            // check the maximum damping seen for this mode
            if (max_g_val < tol)
                modes_to_neglect[i] = true;
        }
        
        // identify the flutter cross-overs. For this, move from
        // higher k to lower k, and handle the k=0 cases separately
        {
            std::map<Real, UGFlutterSolver::UGFlutterSolution*>::const_iterator
            sol_it    = _flutter_solutions.begin(), // first of the pair
            sol_end   = _flutter_solutions.end();
            if (sol_it == sol_end)
                return;
            
            // check if k=0 exists, and identify
            if (fabs(sol_it->second->k_ref()) < tol) { // k = 0
                
                // k=0 makes sense only for divergence roots. Do not use them
                // for crossover points if a finite damping was seen. Hence,
                // do nothing here, and move to the next iterator
                for (unsigned int i=0; i<nvals; i++) {
                    if (!sol_it->second->get_root(i).if_imag_V &&
                        fabs(sol_it->second->get_root(i).g) < tol) {
                        UGFlutterSolver::UGFlutterRootCrossover* cross =
                        new UGFlutterSolver::UGFlutterRootCrossover;
                        cross->crossover_solutions.first = sol_it->second;
                        cross->crossover_solutions.second = sol_it->second;
                        cross->root = &sol_it->second->get_root(i);
                        cross->root_num = i;
                        std::pair<Real, UGFlutterSolver::UGFlutterRootCrossover*>
                        val( cross->root->V, cross);
                        _flutter_crossovers.insert(val);
                    }
                }
            }
        }
        
        // now look for oscillatory roots crossover points in decreasing
        // order of k_ref
        for (unsigned int i=0; i<nvals; i++) {
            std::map<Real, UGFlutterSolver::UGFlutterSolution*>::const_reverse_iterator
            sol_rit    = _flutter_solutions.rbegin(), // first of the pair
            sol_ritp1    = _flutter_solutions.rbegin(), // first of the pair
            sol_rend   = _flutter_solutions.rend();
            if (sol_rit == sol_rend)
                return;

            sol_ritp1++; // increment for the next pair of results
            while (sol_ritp1 != sol_rend) {
                // do not use k_ref = 0, or if the root is invalid
                if (sol_rit->second->get_root(i).if_imag_V ||
                    sol_ritp1->second->get_root(i).if_imag_V ||
                    fabs(sol_rit->second->k_ref()) < tol ||
                    fabs(sol_ritp1->second->k_ref()) < tol ||
                    fabs(sol_rit->second->get_root(i).g) > max_allowable_g ||
                    fabs(sol_ritp1->second->get_root(i).g) > max_allowable_g) {
                    // do nothing
                }
                else if (!modes_to_neglect[i]) { // look for the flutter roots
                    UGFlutterSolver::UGFlutterSolution *lower = sol_rit->second,
                    *upper = sol_ritp1->second;
                    
                    if ((lower->get_root(i).g <= 0.) &&
                        (upper->get_root(i).g > 0.)) {
                        UGFlutterSolver::UGFlutterRootCrossover* cross =
                        new UGFlutterSolver::UGFlutterRootCrossover;
                        cross->crossover_solutions.first  = lower; // -ve g
                        cross->crossover_solutions.second = upper; // +ve g
                        cross->root_num = i;
                        std::pair<Real, UGFlutterSolver::UGFlutterRootCrossover*>
                        val( lower->get_root(i).V, cross);
                        _flutter_crossovers.insert(val);
                    }
                    else if ((lower->get_root(i).g > 0.) &&
                             (upper->get_root(i).g <= 0.)) {
                        UGFlutterSolver::UGFlutterRootCrossover* cross =
                        new UGFlutterSolver::UGFlutterRootCrossover;
                        cross->crossover_solutions.first  = upper; // -ve g
                        cross->crossover_solutions.second = lower; // +ve g
                        cross->root_num = i;
                        std::pair<Real, UGFlutterSolver::UGFlutterRootCrossover*>
                        val( upper->get_root(i).V, cross);
                        _flutter_crossovers.insert(val);
                    }
                }
                
                // increment the pointers for next pair of roots
                sol_rit++;
                sol_ritp1++;
            }
        }
    }
}




std::pair<bool, const FlutterRoot*>
UGFlutterSolver::find_next_root()
{
    // iterate over the cross-over points and calculate the next that has
    // not been evaluated
    std::map<Real, UGFlutterSolver::UGFlutterRootCrossover*>::iterator
    it = _flutter_crossovers.begin(), end = _flutter_crossovers.end();
    while ( it != end)
    {
        UGFlutterSolver::UGFlutterRootCrossover* cross = it->second;
        
        if (!cross->root) {
            const unsigned int root_num = cross->root_num;
            std::pair<bool, UGFlutterSolver::UGFlutterSolution*> bisection_sol =
            bisection_search(cross->crossover_solutions,
                             root_num, 1.0e-5, 10);
            cross->root = &(bisection_sol.second->get_root(root_num));
            
            // now, remove this entry from the _flutter_crossover points and
            // reinsert it with the actual critical velocity
            _flutter_crossovers.erase(it);
            std::pair<Real, UGFlutterSolver::UGFlutterRootCrossover*>
            val(cross->root->V, cross);
            _flutter_crossovers.insert(val);
            return std::pair<bool, const FlutterRoot*> (true, cross->root);
        }
        
        it++;
    }
    
    // if it gets here, no new root was found
    return std::pair<bool, FlutterRoot*> (false, NULL);
}



std::pair<bool, const FlutterRoot*>
UGFlutterSolver::find_critical_root()
{
    // iterate over the cross-over points and calculate the next that has
    // not been evaluated
    std::map<Real, UGFlutterSolver::UGFlutterRootCrossover*>::iterator
    it = _flutter_crossovers.begin(), end = _flutter_crossovers.end();
    
    if (it == end) // no potential cross-over points were identified
        return std::pair<bool, FlutterRoot*> (false, NULL);
    
    // it is possible that once the root has been found, its velocity end up
    // putting it at a higher velocity in the map, so we need to check if
    // the critical root has changed
    while (!it->second->root)
    {
        UGFlutterSolver::UGFlutterRootCrossover* cross = it->second;
        
        if (!cross->root) {
            const unsigned int root_num = cross->root_num;
            std::pair<bool, UGFlutterSolver::UGFlutterSolution*> bisection_sol =
            bisection_search(cross->crossover_solutions,
                             root_num, 1.0e-5, 10);
            cross->root = &(bisection_sol.second->get_root(root_num));
            
            // now, remove this entry from the _flutter_crossover points and
            // reinsert it with the actual critical velocity
            _flutter_crossovers.erase(it);
            std::pair<Real, UGFlutterSolver::UGFlutterRootCrossover*>
            val(cross->root->V, cross);
            _flutter_crossovers.insert(val);
        }
        
        // update the iterator to make sure that this is updated
        it = _flutter_crossovers.begin(); end = _flutter_crossovers.end();
    }

    // if it gets here, then the root was successfully found
    return std::pair<bool, const FlutterRoot*> (true, it->second->root);
}



UGFlutterSolver::UGFlutterSolution*
UGFlutterSolver::analyze(const Real k_ref,
                         const UGFlutterSolver::UGFlutterSolution* prev_sol)
{
    ComplexMatrixX m, k;

    libMesh::out
    << " ====================================================" << std::endl
    << "UG Solution for k_red = "
    << std::setw(10) << k_ref << std::endl;
    
    initialize_matrices(k_ref, m, k);
    LAPACK_ZGGEV ges;
    ges.compute(m, k);
    ges.scale_eigenvectors_to_identity_innerproduct();
    
    // now insert the root
    std::pair<Real, UGFlutterSolver::UGFlutterSolution*>
    val(k_ref, new UGFlutterSolver::UGFlutterSolution((int)m.rows()));
    bool success = _flutter_solutions.insert(val).second;
    libmesh_assert (success); // make sure that it was successfully added
    val.second->init(k_ref, flight_condition->ref_chord, ges);
    val.second->print(_output, _mode_output);
    
    if (prev_sol)
        val.second->sort(*prev_sol);

    libMesh::out
    << "Finished UG Solution" << std::endl
    << " ====================================================" << std::endl;

    
    return val.second;
}


std::pair<bool, UGFlutterSolver::UGFlutterSolution*>
UGFlutterSolver::bisection_search(const std::pair<UGFlutterSolver::UGFlutterSolution*,
                                  UGFlutterSolver::UGFlutterSolution*>& k_ref_sol_range,
                                  const unsigned int root_num,
                                  const Real g_tol,
                                  const unsigned int max_iters)
{
    // assumes that the upper k_val has +ve g val and lower k_val has -ve
    // k_val
    bool found = false;
    Real lower_k = k_ref_sol_range.first->k_ref(),
    lower_g = k_ref_sol_range.first->get_root(root_num).g,
    upper_k = k_ref_sol_range.second->k_ref(),
    upper_g = k_ref_sol_range.second->get_root(root_num).g,
    new_g = 0., new_k = 0.;
    unsigned int n_iters = 0;
    
    UGFlutterSolver::UGFlutterSolution* new_sol = NULL;
    
    while (n_iters < max_iters) {

        new_k = lower_k + (upper_k-lower_k)/(upper_g-lower_g)*(0.-lower_g); // linear interpolation
        new_sol = analyze(new_k, k_ref_sol_range.first);
        const UGFlutterSolver::UGFlutterRoot& root = new_sol->get_root(root_num);
        
        // check if the new damping value
        if (fabs(root.g) <= g_tol)
            return std::pair<bool, UGFlutterSolver::UGFlutterSolution*>
            (true, new_sol);
        
        // update the k_val
        if (root.g < 0.) {
            lower_k = new_k;
            lower_g = root.g;
        }
        else {
            upper_k = new_k;
            upper_g = root.g;
        }
        
        n_iters++;
    }

    return std::pair<bool, UGFlutterSolver::UGFlutterSolution*>
    (false, new_sol); // return false, along with the latest sol
}


void UGFlutterSolver::initialize_matrices(Real k_ref,
                                          ComplexMatrixX& m, // mass & aero
                                          ComplexMatrixX& k) // stiffness
{
    bool has_matrix = false;
    RealMatrixX mat_r;
    
    // stiffness matrix forms the rhs of the eigenvalue problem
    has_matrix = aero_structural_model->get_structural_stiffness_matrix(mat_r);
    k = mat_r.cast<Complex>();
    libmesh_assert(has_matrix);
    
    // combination of mass and aero matrix forms lhs of the eigenvalue problem
    has_matrix = aero_structural_model->get_structural_mass_matrix(mat_r);
    libmesh_assert(has_matrix);
    
    
    has_matrix = aero_structural_model->get_aero_operator_matrix(k_ref, m);
    libmesh_assert(has_matrix);
    
    m *= 0.5 * flight_condition->gas_property.rho;
    mat_r *= pow(k_ref/flight_condition->ref_chord, 2);
    m += mat_r.cast<Complex>();
}


void
UGFlutterSolver::print_sorted_roots(std::ostream* output)
{
    if (!output)
        output = &_output;
    
    std::map<Real, UGFlutterSolver::UGFlutterSolution*>::const_iterator
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
        << "** UG Root # "
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
            const UGFlutterSolver::UGFlutterRoot& root =
            sol_it->second->get_root(i);

            *output
            << std::setw(15) << root.k_ref
            << std::setw(15) << std::real(root.root)
            << std::setw(15) << std::imag(root.root)
            << std::setw(15) << root.g;
            
            // check if the root might be the odd artifact of the UG method
            if (!root.if_imag_V)
                *output << std::setw(15) << root.V;
            else
                *output << "** " << std::setw(12) << root.V;
            *output
            << std::setw(15) << root.omega << std::endl;
        }
        *output << std::endl << std::endl;
    }
    
    
    // write the roots identified using iterative search technique
    unsigned int nroots = this->n_roots_found();
    *output << std::endl
    << "n critical roots identified: " << nroots << std::endl;
    for (unsigned int i=0; i<nroots; i++)
    {
        const FlutterRoot& root = this->get_root(i);
        *output
        << "** Root : " << std::setw(5) << i << " **" << std::endl
        << "g      = " << std::setw(15) << root.g << std::endl
        << "V      = " << std::setw(15) << root.V << std::endl
        << "omega  = " << std::setw(15) << root.omega << std::endl
        << "k_ref  = " << std::setw(15) << root.k_ref << std::endl
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
UGFlutterSolver::print_crossover_points(std::ostream* output)
{
    if (!output)
        output = &_output;

    *output << "n crossover points found: "
    << std::setw(5) << _flutter_crossovers.size() << std::endl;
    
    std::multimap<Real, UGFlutterSolver::UGFlutterRootCrossover*>::const_iterator
    it = _flutter_crossovers.begin(), end = _flutter_crossovers.end();
    
    unsigned int i=0;
    
    for ( ; it != end; it++) {
        *output << "** Point : " << std::setw(5) << i << " **" << std::endl;
        it->second->print(*output);
        *output << std::endl;
        i++;
    }
}



