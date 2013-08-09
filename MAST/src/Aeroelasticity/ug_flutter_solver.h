//
//  ug_flutter_solver.h
//  MAST
//
//  Created by Manav Bhatia on 7/27/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef MAST_ug_flutter_solver_h
#define MAST_ug_flutter_solver_h

// C++
#include <iostream>
#include <iomanip>

// MAST includes
#include "Aeroelasticity/flutter_solver_base.h"
#include "Solvers/lapack_interface.h"


class UGFlutterSolver: public FlutterSolverBase
{
public:
    UGFlutterSolver():
    FlutterSolverBase(),
    k_ref_range(std::pair<Real, Real>(0.01, 0.75)),
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

    void evaluate_roots(const Real k_ref,
                        const ComplexVectorX& num, const ComplexVectorX& den);
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
    
    ComplexMatrixX m, k;
    
    while (current_k_ref >= active_k_ref_range.first)
    {
        initialize_matrices(current_k_ref, m, k);
        LAPACK_ZGGEV ges;
        ges.compute(m, k);
        evaluate_roots(current_k_ref, ges.alphas(), ges.betas());
        
        current_k_ref -= delta_k_ref;
    }
    
    return found_root;
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
    has_matrix = aero_structural_model->get_structural_mass_matrix(mat_r);
    libmesh_assert(has_matrix);
    
    
    has_matrix = aero_structural_model->get_aero_operator_matrix(k_ref, m);
    libmesh_assert(has_matrix);

    m *= 0.5 * flight_condition->gas_property.rho;
    mat_r *= pow(k_ref/flight_condition->ref_chord, 2);
    m += mat_r.cast<Complex>();
}


void
UGFlutterSolver::evaluate_roots(const Real k_ref,
                                const ComplexVectorX & num,
                                const ComplexVectorX& den)
{
    libmesh_assert(num.rows() == den.rows());
    libmesh_assert(num.cols() == den.cols());
    
    unsigned int nvals = std::max(num.rows(), num.cols());
    
    Real g, V, omega;
    Complex eig;
    bool if_imag_v = false;
    
    std::cout << "k = " << k_ref << std::endl
    << std::setw(5) << "#"
    << std::setw(15) << "Re"
    << std::setw(15) << "Im"
    << std::setw(15) << "g"
    << std::setw(15) << "V"
    << std::setw(15) << "omega" << std::endl;
    
    for (unsigned int i=0; i<nvals; i++)
    {
        g = 0.; V= 0., omega = 0.;
        eig = 0.;
        if (std::abs(den(i)) > 0.)
        {
            eig = num(i)/den(i);
            if (std::real(eig) > 0.)
            {
                V     = sqrt(1./std::real(eig));
                g     = std::imag(eig)/std::real(eig);
                omega = k_ref*V/flight_condition->ref_chord;
            }
            else
                if_imag_v = true;
        }

        std::cout
        << std::setw(5) << i
        << std::setw(15) << std::real(eig)
        << std::setw(15) << std::imag(eig)
        << std::setw(15) << g
        << std::setw(15) << V
        << std::setw(15) << omega;
        if (if_imag_v)
            std::cout << std::setw(20) << "[Imag. Velocity]";
        std::cout << std::endl;
        
        if_imag_v = false;
    }
    
}

#endif
