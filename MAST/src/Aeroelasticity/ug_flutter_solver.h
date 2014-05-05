//
//  ug_flutter_solver.h
//  MAST
//
//  Created by Manav Bhatia on 7/27/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

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
                          const ComplexVectorX& eig_vec);
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
        virtual MAST::FlutterRootBase* build_flutter_root();

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
         *    initializes the matrices for the specified k_ref. UG does not account
         *    for structural damping.
         */
        void initialize_matrices(Real k_ref,
                                 ComplexMatrixX& m, // mass & aero
                                 ComplexMatrixX& k); // aero operator
    };
}

#endif
