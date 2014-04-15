//
//  time_domain_flutter_solver.h
//  MAST
//
//  Created by Manav Bhatia on 4/15/14.
//  Copyright (c) 2014 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_time_domain_flutter_solver_h__
#define __MAST_time_domain_flutter_solver_h__

// MAST includes
#include "Aeroelasticity/flutter_solver_base.h"

namespace MAST {
    
    class TimeDomainFlutterSolver: public MAST::FlutterSolverBase
    {
    public:
        TimeDomainFlutterSolver():
        MAST::FlutterSolverBase() {
            this->ref_val_range = std::pair<Real, Real>(0., 0.);
            this->n_ref_val_divs = 1;
        }
        
        
        virtual ~TimeDomainFlutterSolver();
        
        
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
         *    initializes the matrices for the specified velocity.
         */
        void initialize_matrices(Real ref_val,
                                 RealMatrixX& a,  // LHS matrix operator
                                 RealMatrixX& b); // RHS matrix operator
    };
}

#endif // __MAST_time_domain_flutter_solver_h__
