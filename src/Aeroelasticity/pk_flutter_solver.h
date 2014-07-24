//
//  pk_flutter_solver.h
//  MAST
//
//  Created by Manav Bhatia on 9/18/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST__pk_flutter_solver__
#define __MAST__pk_flutter_solver__



// C++
#include <iostream>
#include <map>

// MAST includes
#include "Aeroelasticity/flutter_solver_base.h"
#include "Solvers/lapack_interface.h"

class PKFlutterSolver: public FlutterSolverBase
{
public:
    PKFlutterSolver():
    FlutterSolverBase(),
    k_ref_range(std::pair<Real, Real>(0., 0.35)),
    n_k_divs(10),
    _previous_k_ref(0.)
    {}
    
    
    virtual ~PKFlutterSolver();
    
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
     *   Scans for flutter roots in the range specified, and identified the
     *   divergence (if k_ref = 0. is specified) and flutter crossover points.
     *   The roots are organized in terms of increasing velocity.
     */
    virtual void scan_for_roots();
    
    
    /*!
     *    finds the number of critical points already identified in the
     *    procedure.
     */
    virtual unsigned int n_roots_found() const;
    
    
    /*!
     *   returns the \par n th root in terms of ascending velocity that is
     *   found by the solver
     */
    const FlutterRoot& get_root(const unsigned int n) const;
    
    
    /*!
     *   Looks through the list of flutter cross-over points and iteratively
     *   zooms in to find the cross-over point. This should be called only
     *   after scan_for_roots() has been called. Potential cross-over points
     *   are sorted with increasing velocity, and this method will attempt to
     *   identify the next critical root in the order.
     */
    virtual std::pair<bool, const FlutterRoot*> find_next_root();
    
    /*!
     *   Prints the sorted roots to the \par output
     */
    void print_sorted_roots(std::ostream* output = NULL);
    
    
    class PKFlutterRoot: public FlutterRoot
    {
    public:
        PKFlutterRoot(): FlutterRoot(), if_nonphysical_root(false) { }
        void init(const Real k, const Real b_ref,
                  const Complex num, const Complex den,
                  const ComplexMatrixX& Bmat,
                  const ComplexVectorX& eig_vec);
        bool if_nonphysical_root;
    };
    
    
    class PKFlutterSolution
    {
    public:
        PKFlutterSolution(unsigned int n):
        _k_ref(0.) {
            _roots.resize(n);
            for (unsigned int i=0; i<n; i++)
                _roots[i] = new PKFlutterSolver::PKFlutterRoot();
        }
        
        /*!
         *    delete the flutter root objects
         */
        ~PKFlutterSolution() {
            std::vector<PKFlutterSolver::PKFlutterRoot*>::iterator
            it = _roots.begin();
            for ( ; it != _roots.end(); it++)
                delete *it;
        }
        
        /*!
         *   the reduced frequency for this solution
         */
        Real k_ref() const {
            return _k_ref;
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
        const PKFlutterRoot& get_root(const unsigned int i) const {
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
        void sort(const PKFlutterSolver::PKFlutterSolution& sol);
        
        /*!
         *    prints the data and modes from this solution
         */
        void print(std::ostream& output, std::ostream& mode_output);
        
    protected:
        Real _k_ref;
        std::vector<PKFlutterSolver::PKFlutterRoot*> _roots;
    };
    
    
    class PKFlutterRootCrossover
    {
    public:
        PKFlutterRootCrossover():
        crossover_solutions(NULL, NULL),
        root_num(0),
        root(NULL)
        { }
        
        std::pair<PKFlutterSolver::PKFlutterSolution*, PKFlutterSolver::PKFlutterSolution*>
        crossover_solutions;
        
        unsigned int root_num;
        
        const PKFlutterSolver::PKFlutterRoot* root;
    };
    
protected:
    
    /*!
     *   map of reduced frequency vs flutter solutions
     */
    std::map<Real, PKFlutterSolver::PKFlutterSolution*> _flutter_solutions;
    
    /*!
     *   the map of flutter crossover points versus average velocity of the
     *   two bounding roots
     */
    std::multimap<Real, PKFlutterSolver::PKFlutterRootCrossover*> _flutter_crossovers;
    
    /*!
     *   performs an eigensolution at the specified reduced frequency, and
     *   sort the roots based on the provided solution pointer. If the
     *   pointer is NULL, then no sorting is performed
     */
    PKFlutterSolver::PKFlutterSolution*
    analyze(const Real k_ref,
            const PKFlutterSolver::PKFlutterSolution* prev_sol=NULL);
    
    /*!
     *    bisection method search
     */
    std::pair<bool, PKFlutterSolver::PKFlutterSolution*>
    bisection_search(const std::pair<PKFlutterSolver::PKFlutterSolution*,
                     PKFlutterSolver::PKFlutterSolution*>& k_ref_sol_range,
                     const unsigned int root_num,
                     const Real g_tol,
                     const unsigned int max_iters);
    
    /*!
     *    initializes the matrices for the specified k_ref. UG does not account
     *    for structural damping.
     */
    void initialize_matrices(Real k_ref,
                             ComplexMatrixX& m, // mass & aero
                             ComplexMatrixX& k); // aero operator
};



#endif /* defined(__MAST__pk_flutter_solver__) */
