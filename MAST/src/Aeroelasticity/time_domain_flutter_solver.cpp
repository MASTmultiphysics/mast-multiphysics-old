//
//  time_domain_flutter_solver.cpp
//  RealSolver
//
//  Created by Manav Bhatia on 4/15/14.
//  Copyright (c) 2014 Manav Bhatia. All rights reserved.
//

// C++ includes
#include <iomanip>

// MAST includes
#include "Aeroelasticity/time_domain_flutter_solver.h"



MAST::TimeDomainFlutterSolver::~TimeDomainFlutterSolver()
{ }





void
MAST::TimeDomainFlutterSolver::_identify_crossover_points()
{
    // if the initial scanning has not been done, then do it now
    const Real tol = 1.0e-5, max_allowable_g = 0.75;
    
    const unsigned int nvals = _flutter_solutions.begin()->second->n_roots();
    // make sure that the solution has been generated
    libmesh_assert(nvals);
    
    //
    // for some cases some roots trail along the g=0 axis
    // and should not be considered as flutter. These are simply
    // modes where aerodynamics do not provide any damping.
    // These modes will not have a damping more than tolerance
    //
    
    std::vector<bool> modes_to_neglect(nvals);
    std::fill(modes_to_neglect.begin(),
              modes_to_neglect.end(), false);
    
    // look for the max g val for a mode, which will indicate if the
    // mode is undamped or not
    for (unsigned int i=0; i<nvals; i++) {
        std::map<Real, MAST::FlutterSolutionBase*>::const_iterator
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
        std::map<Real, MAST::FlutterSolutionBase*>::const_iterator
        sol_it    = _flutter_solutions.begin(), // first of the pair
        sol_end   = _flutter_solutions.end();
        if (sol_it == sol_end)
            return;
        
        // check if k=0 exists, and identify
        if (fabs(sol_it->second->ref_val()) < tol) { // k = 0
            
            // k=0 makes sense only for divergence roots. Do not use them
            // for crossover points if a finite damping was seen. Hence,
            // do nothing here, and move to the next iterator
            for (unsigned int i=0; i<nvals; i++) {
                if (!sol_it->second->get_root(i).if_nonphysical_root &&
                    fabs(sol_it->second->get_root(i).g) < tol) {
                    MAST::FlutterRootCrossoverBase* cross =
                    new MAST::FlutterRootCrossoverBase;
                    cross->crossover_solutions.first = sol_it->second;
                    cross->crossover_solutions.second = sol_it->second;
                    cross->root = &sol_it->second->get_root(i);
                    cross->root_num = i;
                    std::pair<Real, MAST::FlutterRootCrossoverBase*>
                    val( cross->root->V, cross);
                    _flutter_crossovers.insert(val);
                }
            }
        }
    }
    
    // now look for oscillatory roots crossover points in decreasing
    // order of k_ref
    for (unsigned int i=0; i<nvals; i++) {
        std::map<Real, MAST::FlutterSolutionBase*>::const_reverse_iterator
        sol_rit    = _flutter_solutions.rbegin(), // first of the pair
        sol_ritp1    = _flutter_solutions.rbegin(), // first of the pair
        sol_rend   = _flutter_solutions.rend();
        if (sol_rit == sol_rend)
            return;
        
        sol_ritp1++; // increment for the next pair of results
        while (sol_ritp1 != sol_rend) {
            // do not use k_ref = 0, or if the root is invalid
            if (sol_rit->second->get_root(i).if_nonphysical_root ||
                sol_ritp1->second->get_root(i).if_nonphysical_root ||
                fabs(sol_rit->second->ref_val()) < tol ||
                fabs(sol_ritp1->second->ref_val()) < tol ||
                fabs(sol_rit->second->get_root(i).g) > max_allowable_g ||
                fabs(sol_ritp1->second->get_root(i).g) > max_allowable_g) {
                // do nothing
            }
            else if (!modes_to_neglect[i]) { // look for the flutter roots
                MAST::FlutterSolutionBase *lower = sol_rit->second,
                *upper = sol_ritp1->second;
                
                if ((lower->get_root(i).g <= 0.) &&
                    (upper->get_root(i).g > 0.)) {
                    MAST::FlutterRootCrossoverBase* cross =
                    new MAST::FlutterRootCrossoverBase;
                    cross->crossover_solutions.first  = lower; // -ve g
                    cross->crossover_solutions.second = upper; // +ve g
                    cross->root_num = i;
                    std::pair<Real, MAST::FlutterRootCrossoverBase*>
                    val( lower->get_root(i).V, cross);
                    _flutter_crossovers.insert(val);
                }
                else if ((lower->get_root(i).g > 0.) &&
                         (upper->get_root(i).g <= 0.)) {
                    MAST::FlutterRootCrossoverBase* cross =
                    new MAST::FlutterRootCrossoverBase;
                    cross->crossover_solutions.first  = upper; // -ve g
                    cross->crossover_solutions.second = lower; // +ve g
                    cross->root_num = i;
                    std::pair<Real, MAST::FlutterRootCrossoverBase*>
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



MAST::FlutterSolutionBase*
MAST::TimeDomainFlutterSolver::analyze(const Real ref_val,
                                       const MAST::FlutterSolutionBase* prev_sol) {
    RealMatrixX a, b;
    
    libMesh::out
    << " ====================================================" << std::endl
    << "Eigensolution for V = "
    << std::setw(10) << ref_val << std::endl;
    
    initialize_matrices(ref_val, a, b);
    LAPACK_DGGEV ges;
    ges.compute(a, b);
    ges.scale_eigenvectors_to_identity_innerproduct();
    
    // now insert the root
    std::pair<Real, MAST::FlutterSolutionBase*>
    val(ref_val, new MAST::TimeDomainFlutterSolution(*this));
    bool success = _flutter_solutions.insert(val).second;
    libmesh_assert (success); // make sure that it was successfully added
    dynamic_cast<MAST::TimeDomainFlutterSolution*>(val.second)->init
    (ref_val, flight_condition->ref_chord, ges);
    val.second->print(_output, _mode_output);
    
    if (prev_sol)
        val.second->sort(*prev_sol);
    
    libMesh::out
    << "Finished Eigensolution" << std::endl
    << " ====================================================" << std::endl;
    
    
    return val.second;
}



void MAST::TimeDomainFlutterSolver::initialize_matrices(Real ref_val,
                                                        RealMatrixX& a, // LHS
                                                        RealMatrixX& b) // RHS
{
    bool has_matrix = false;
    RealMatrixX mat1, mat2;
    Real q0 = flight_condition->q0();
    
    // mass matrix forms the LHS of the eigenvalue problem
    has_matrix = aero_structural_model->get_structural_mass_matrix(mat1);
    const unsigned int n = (int)mat1.rows();
    a.block(n,n,n,n) = mat1;
    libmesh_assert(has_matrix);


    // structural and aerodynamic stiffness matrices
    has_matrix = aero_structural_model->get_structural_stiffness_matrix(mat1);
    libmesh_assert(has_matrix);
    has_matrix = aero_structural_model->get_aero_stiffness_matrix(mat2);
    libmesh_assert(has_matrix);
    b.block(n, 0, n, n) = (q0*mat2-mat1);
    
    // structural and aerodynamic stiffness matrices
    has_matrix = aero_structural_model->get_structural_damping_matrix(mat1);
    libmesh_assert(has_matrix);
    has_matrix = aero_structural_model->get_aero_damping_matrix(mat2);
    libmesh_assert(has_matrix);
    b.block(n, 0, n, n) = (q0*mat2-mat1);

    // finally set the identity blocks in the matrix
    for (unsigned int i=0; i<n; i++) {
        a(i,  i) = 1.;
        b(i,i+n) = 1.;
    }
}



