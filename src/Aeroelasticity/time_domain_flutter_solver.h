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

#ifndef __MAST_time_domain_flutter_solver_h__
#define __MAST_time_domain_flutter_solver_h__

// MAST includes
#include "Aeroelasticity/flutter_solver_base.h"

namespace MAST {
    
    
    
    class TimeDomainFlutterRoot: public MAST::FlutterRootBase {
    public:
        TimeDomainFlutterRoot():
        MAST::FlutterRootBase()
        { }
        
        virtual ~TimeDomainFlutterRoot() {}
        
        virtual void init(const Real ref_val, const Real b_ref,
                          const Complex num,
                          const Complex den,
                          const RealMatrixX& Bmat,
                          const ComplexVectorX& evec_right,
                          const ComplexVectorX& evec_left);
    };
    

    
    class TimeDomainFlutterSolution: public MAST::FlutterSolutionBase {
    public:
        
        TimeDomainFlutterSolution():
        MAST::FlutterSolutionBase()
        { }
        
        virtual ~TimeDomainFlutterSolution() {}
        
        /*!
         *   initializes the flutter solution from an eigensolution
         */
        virtual void init (const MAST::FlutterSolverBase& solver,
                           const Real v_ref,
                           const Real bref,
                           const LAPACK_DGGEV& eig_sol);
    };
    
    
    

    class TimeDomainFlutterSolver: public MAST::FlutterSolverBase
    {
    public:
        TimeDomainFlutterSolver():
        MAST::FlutterSolverBase(),
        vel_range(std::pair<Real, Real>(0., 0.)),
        n_vel_divs(1)
        {}
        
        
        virtual ~TimeDomainFlutterSolver();
        
        /*!
         *   range of reference values within which to find flutter roots
         */
        std::pair<Real, Real> vel_range;
        
        /*!
         *    number of division in the reference value range for initial
         *    scanning
         */
        unsigned int n_vel_divs;

        
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
        virtual std::auto_ptr<MAST::FlutterSolutionBase>
        analyze(const Real ref_val,
                const MAST::FlutterSolutionBase* prev_sol=NULL);
        
        /*!
         *    initializes the matrices for the specified velocity.
         */
        void initialize_matrices(Real v_ref,
                                 RealMatrixX& a,  // LHS matrix operator
                                 RealMatrixX& b); // RHS matrix operator
    };
}

#endif // __MAST_time_domain_flutter_solver_h__
