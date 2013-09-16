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
#include <map>

// MAST includes
#include "Aeroelasticity/flutter_solver_base.h"
#include "Solvers/lapack_interface.h"



class UGFlutterSolver: public FlutterSolverBase
{
public:
    UGFlutterSolver():
    FlutterSolverBase(),
    k_ref_range(std::pair<Real, Real>(0., 0.35)),
    n_k_divs(10),
    _previous_k_ref(0.)
    {}
    
    
    virtual ~UGFlutterSolver();
    
    /*!
     *   range of reduced frequencies within which to find flutter roots
     *   Default range is 0.1 to 1.0
     */
    std::pair<Real, Real> k_ref_range;

    /*!
     *    number of division in the reduced frequency range for initial scanning
     */
    unsigned int n_k_divs;
    
    /*!
     *   Finds and appends the root to the vector of flutter roots. Returns
     *   \p true if a root was successfully found, \p false otherwise.
     */
    virtual bool find_next_root();

    
    class UGFlutterRoot: public FlutterRoot
    {
    public:
        UGFlutterRoot(): FlutterRoot(), if_imag_V(false) { }
        void init(const Real k, const Real b_ref,
                  const Complex num, const Complex den,
                  const ComplexMatrixX& Bmat,
                  const ComplexVectorX& eig_vec);
        bool if_imag_V;
    };
    
    
    class UGFlutterSolution
    {
    public:
        UGFlutterSolution(unsigned int n):
        _k_ref(0.) {
            _roots.resize(n);
            std::fill(_roots.begin(), _roots.end(),
                      new UGFlutterSolver::UGFlutterRoot() );
        }

        /*!
         *    delete the flutter root objects
         */
        ~UGFlutterSolution() {
            std::vector<UGFlutterSolver::UGFlutterRoot*>::iterator
            it = _roots.begin();
            for ( ; it != _roots.end(); it++)
                delete *it;
        }

        /*!
         *   number of roots in this solution
         */
        unsigned int n_roots() const{
            return _roots.size();
        }
        
        /*!
         *    returns the root
         */
        const UGFlutterRoot& get_root(const unsigned int i) const {
            libmesh_assert_less(i, _roots.size());
            return *_roots[i];
        }
        
        /*!
         *   initializes the UG flutter solution for an eigensolution
         */
        void init (const Real kref, const Real bref,
                   const LAPACK_ZGGEV& eig_sol);
        
        /*!
         *    sort this root with respect to the given solution from a previous
         *    eigen solution. This method relies on the modal participation.
         *    Flutter roots from previous and current solutions with highest
         *    dot product of modal participation vector are considered to be
         *    similar.
         */
        void sort(const UGFlutterSolver::UGFlutterSolution& sol);
        
    protected:
        Real _k_ref;
        std::vector<UGFlutterSolver::UGFlutterRoot*> _roots;
    };

protected:
    
    /*!
     *   reduced frequency for the previous evaluation
     */
    Real _previous_k_ref;

    /*!
     *   map of reduced frequency vs flutter solutions
     */
    std::map<Real, UGFlutterSolver::UGFlutterSolution*> _flutter_solutions;
    
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
                        const LAPACK_ZGGEV& ges);
};

#endif
