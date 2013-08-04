//
//  ug_flutter_solver.h
//  MAST
//
//  Created by Manav Bhatia on 7/27/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef MAST_ug_flutter_solver_h
#define MAST_ug_flutter_solver_h

// MAST includes
#include "Aeroelasticity/flutter_solver_base.h"


class UGFlutterSolver: public FlutterSolverBase
{
public:
    UGFlutterSolver():
    FlutterSolverBase(),
    k_ref_range(std::pair<Real, Real>(0.1, 1.)),
    _previous_k_ref(0.)
    {}
    
    
    virtual ~UGFlutterSolver()
    {}
    
    /*!
     *   range of reduced frequencies within which to find flutter roots
     *   Default range is 0.1 to 1.0
     */
    std::pair<Real, Real> k_ref_range;
    
    /*!
     *   Finds and appends the root to the vector of flutter roots. Returns
     *   \p true if a root was successfully found, \p false otherwise.
     */
    virtual bool find_next_root();
    
protected:
    
    /*!
     *   reduced frequency for the previous evaluation
     */
    Real _previous_k_ref;
    
    /*!
     *    bisection method search
     */
    bool bisection_search(const std::pair<Real, Real>& active_k_ref_range);
    
    /*!
     *    initializes the matrices for the specified k_ref. UG does not account
     *    for structural damping.
     */
    void initialize_matrices(Real k_ref,
                             ComplexMatrixX& m, // mass & aero
                             ComplexMatrixX& k); // aero operator

};


bool UGFlutterSolver::find_next_root()
{
    // active range in which to search
    std::pair<Real, Real> active_k_ref_range(0.,0.);
    
    // if this is the first solve, then the active range is same as the
    // k_ref range
    if (n_roots_found() == 0)
        active_k_ref_range = k_ref_range;
    else
        active_k_ref_range = std::pair<Real, Real>
        (k_ref_range.first, _previous_k_ref);
    
    // now march from the upper limit to the lower to find the roots
    Real current_k_ref = active_k_ref_range.second,
    delta_k_ref = (active_k_ref_range.second-active_k_ref_range.first)/20;
    bool continue_finding = true, found_root = false;
    
    unsigned int n_positive_damping = 0;
    
    while (continue_finding)
    {
        // to be programmed
        libmesh_assert(false);
        
        if (n_positive_damping == 0)
        {
            // if a root was not found, then move to the next k_ref
            active_k_ref_range.second -= delta_k_ref;
            current_k_ref -= delta_k_ref;
        }
        else
        {
            // if a root was identified, then use bisection method for accurate
            // root estimation
            std::pair<Real, Real> bisection_k_ref_range
            (active_k_ref_range.second,
             active_k_ref_range.second + delta_k_ref);
            
            // now identify the root
            found_root = bisection_search(bisection_k_ref_range);
            continue_finding = false;
        }
    }
    
}



bool UGFlutterSolver::bisection_search(const std::pair<Real, Real>& active_k_ref_range)
{
    // to be programmed
    libmesh_assert(false);
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
    has_matrix = aero_structural_model->get_structural_stiffness_matrix(mat_r);
    libmesh_assert(has_matrix);
    
    
    has_matrix = aero_structural_model->get_aero_operator_matrix(k_ref, m);
    libmesh_assert(has_matrix);

    m *= 0.5 * flight_condition.gas_property.rho;
    mat_r *= pow(k_ref/flight_condition.ref_chord, 2);
    m += mat_r.cast<Complex>();
}

#endif
