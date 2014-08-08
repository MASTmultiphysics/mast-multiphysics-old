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

#ifndef __MAST__pk_flutter_solver__
#define __MAST__pk_flutter_solver__



// C++
#include <iostream>
#include <map>

// MAST includes
#include "Aeroelasticity/flutter_solver_base.h"


namespace MAST {
    
    class PKFlutterSolution: public MAST::FrequencyDomainFlutterSolution {
    public:
        
        PKFlutterSolution():
        MAST::FrequencyDomainFlutterSolution()
        { }
        
        
        virtual ~PKFlutterSolution() {}
        
        /*!
         *   initializes the flutter solution from an eigensolution
         */
        virtual void init (const MAST::FlutterSolverBase& solver,
                           const Real k_red,
                           const Real v_ref,
                           const Real bref,
                           const LAPACK_ZGGEV& eig_sol);
        
    };

    
    class PKFlutterRoot: public MAST::FrequencyDomainFlutterRoot {
    public:
        PKFlutterRoot():
        MAST::FrequencyDomainFlutterRoot()
        { }
        
        virtual ~PKFlutterRoot() {}
        
        virtual void init(const Real k_red_ref,
                          const Real V_ref,
                          const Real b_ref,
                          const Complex num, const Complex den,
                          const ComplexMatrixX& Kmat,
                          const ComplexVectorX& evec_right,
                          const ComplexVectorX& evec_left);
    };
    
    
    
    
    class PKFlutterSolver: public MAST::FlutterSolverBase
    {
    public:
        PKFlutterSolver():
        MAST::FlutterSolverBase(),
        v_ref_range(std::pair<Real, Real>(0., 0.)),
        n_v_ref_divs(1),
        k_red_range(std::pair<Real, Real>(0., 0.)),
        n_k_red_divs(1)
        { }
        
        
        virtual ~PKFlutterSolver();
        
        
        /*!
         *   range of reference values within which to find flutter roots
         */
        std::pair<Real, Real> v_ref_range;
        
        /*!
         *    number of division in the reference value range for initial
         *    scanning
         */
        unsigned int n_v_ref_divs;

        /*!
         *   range of reference values within which to find flutter roots
         */
        std::pair<Real, Real> k_red_range;
        
        /*!
         *    number of division in the reference value range for initial
         *    scanning
         */
        unsigned int n_k_red_divs;
        
        
        virtual void scan_for_roots();
        
        /*!
         *    creates a new flutter root and returns pointer to it.
         */
        virtual std::auto_ptr<MAST::FlutterRootBase> build_flutter_root() const;
        
        /*!
         *   Calculate the sensitivity of the flutter root with respect to the
         *   \par i^th parameter in params
         */
        virtual void calculate_sensitivity(MAST::FlutterRootBase& root,
                                           const libMesh::ParameterVector& params,
                                           const unsigned int i);
        
    protected:
        
        /*!
         *   identifies all cross-over and divergence points from analyzed
         *   roots
         */
        virtual void _identify_crossover_points();
        
        
        void _insert_new_solution(const Real k_red_ref,
                                  MAST::FlutterSolutionBase* sol);

        
        virtual std::pair<bool, MAST::FlutterSolutionBase*>
        bisection_search(const std::pair<MAST::FlutterSolutionBase*,
                         MAST::FlutterSolutionBase*>& ref_sol_range,
                         const unsigned int root_num,
                         const Real g_tol,
                         const unsigned int max_iters);
        
        
        
        /*!
         *    Newton method to look for cross-over point method search
         */
        virtual std::pair<bool, MAST::FlutterSolutionBase*>
        newton_search(const MAST::FlutterSolutionBase& init_sol,
                      const unsigned int root_num,
                      const Real tol,
                      const unsigned int max_iters);
        
        /*!
         *   performs an eigensolution at the specified reduced frequency, and
         *   sort the roots based on the provided solution pointer. If the
         *   pointer is NULL, then no sorting is performed
         */
        virtual std::auto_ptr<MAST::FlutterSolutionBase>
        analyze(const Real k_red,
                const Real v_ref,
                const MAST::FlutterSolutionBase* prev_sol=NULL);
        
        
        
        /*!
         *    initializes the matrices for the specified k_red. UG does not account
         *    for structural damping.
         */
        void initialize_matrices(const Real k_red,
                                 const Real v_ref,
                                 ComplexMatrixX& m, // mass & aero
                                 ComplexMatrixX& k); // aero operator
        
        /*!
         *    initializes the matrices for the specified k_red. UG does not account
         *    for structural damping.
         */
        void
        initialize_matrix_sensitivity_for_param(const libMesh::ParameterVector& params,
                                                unsigned int p,
                                                const Real k_red,
                                                const Real v_ref,
                                                ComplexMatrixX& m, // mass & aero
                                                ComplexMatrixX& k); // aero operator
        
        /*!
         *    initializes the matrices for the specified k_red. UG does not account
         *    for structural damping.
         */
        void
        initialize_matrix_sensitivity_for_reduced_freq(const Real k_red,
                                                       const Real v_ref,
                                                       ComplexMatrixX& m, // mass & aero
                                                       ComplexMatrixX& k); // aero operator
        
        /*!
         *    initializes the matrices for the specified k_red. UG does not account
         *    for structural damping.
         */
        void
        initialize_matrix_sensitivity_for_V_ref(const Real k_red,
                                                const Real v_ref,
                                                ComplexMatrixX& m, // mass & aero
                                                ComplexMatrixX& k); // aero operator
        
    };
}


#endif /* defined(__MAST__pk_flutter_solver__) */
