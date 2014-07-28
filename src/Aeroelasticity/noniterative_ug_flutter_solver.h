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

#ifndef __MAST_iterative_ug_flutter_solver_h__
#define __MAST_iterative_ug_flutter_solver_h__

// MAST includes
#include "Aeroelasticity/ug_flutter_solver.h"


namespace MAST {
    
    class NoniterativeUGFlutterSolver:
    public MAST::UGFlutterSolver {
        
    public:
        
        NoniterativeUGFlutterSolver() :
        MAST::UGFlutterSolver(),
        v_ref_range(std::pair<Real, Real>(0., 0.)),
        n_v_ref_divs(1)
        {}
        
        ~NoniterativeUGFlutterSolver()
        {}
        
        /*!
         *   range of reference values within which to find flutter roots
         */
        std::pair<Real, Real> v_ref_range;
        
        /*!
         *    number of division in the reference value range for initial
         *    scanning
         */
        unsigned int n_v_ref_divs;

        
        virtual void scan_for_roots();

        std::pair<bool, MAST::FlutterSolutionBase*>
        bisection_search(const std::pair<MAST::FlutterSolutionBase*,
                         MAST::FlutterSolutionBase*>& ref_sol_range,
                         const unsigned int root_num,
                         const Real g_tol,
                         const unsigned int max_iters);
        
    protected:
        
        bool _insert_new_solution(const Real v_ref,
                                  MAST::FlutterSolutionBase& sol);
        
    };
}



#endif // __MAST_iterative_ug_flutter_solver_h__
