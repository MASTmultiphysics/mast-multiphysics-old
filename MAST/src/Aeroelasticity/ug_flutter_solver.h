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

#ifndef __MAST_ug_flutter_solver_h__
#define __MAST_ug_flutter_solver_h__

// MAST includes
#include "Aeroelasticity/flutter_solver_base.h"

namespace MAST {
    
    
    class UGFlutterRoot: public MAST::FrequencyDomainFlutterRoot {
    public:
        UGFlutterRoot():
        MAST::FrequencyDomainFlutterRoot()
        { }
        
        virtual ~UGFlutterRoot() {}
        
        virtual void init(const Real ref_val, const Real b_ref,
                          const Complex num, const Complex den,
                          const ComplexMatrixX& Bmat,
                          const ComplexVectorX& evec_right,
                          const ComplexVectorX& evec_left);
    };

    
    
    
    class UGFlutterSolver: public MAST::FlutterSolverBase
    {
    public:
        UGFlutterSolver():
        MAST::FlutterSolverBase() {
            this->ref_val_range = std::pair<Real, Real>(0., 0.35);
            this->n_ref_val_divs = 10;
        }
        
        
        virtual ~UGFlutterSolver();
        
        
        /*!
         *    creates a new flutter root and returns pointer to it.
         */
        virtual std::auto_ptr<MAST::FlutterRootBase> build_flutter_root();

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
        
        /*!
         *   performs an eigensolution at the specified reduced frequency, and
         *   sort the roots based on the provided solution pointer. If the
         *   pointer is NULL, then no sorting is performed
         */
        virtual MAST::FlutterSolutionBase*
        analyze(const Real ref_val,
                const MAST::FlutterSolutionBase* prev_sol=NULL);

        
        
        /*!
         *    initializes the matrices for the specified k_red. UG does not account
         *    for structural damping.
         */
        void initialize_matrices(Real k_red,
                                 ComplexMatrixX& m, // mass & aero
                                 ComplexMatrixX& k); // aero operator

        /*!
         *    initializes the matrices for the specified k_red. UG does not account
         *    for structural damping.
         */
        void initialize_matrix_sensitivity_for_param(const libMesh::ParameterVector& params,
                                                     unsigned int p,
                                                     Real k_red,
                                                     ComplexMatrixX& m, // mass & aero
                                                     ComplexMatrixX& k); // aero operator
        
        /*!
         *    initializes the matrices for the specified k_red. UG does not account
         *    for structural damping.
         */
        void initialize_matrix_sensitivity_for_reduced_freq(Real k_red,
                                                            ComplexMatrixX& m, // mass & aero
                                                            ComplexMatrixX& k); // aero operator

    };
}

#endif
