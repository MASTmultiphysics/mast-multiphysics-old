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

#ifndef __MAST__flutter_solver_base__
#define __MAST__flutter_solver_base__

// C++ includes
#include <vector>
#include <memory>
#include <fstream>

// MAST includes
#include "Base/MAST_data_types.h"
#include "Aeroelasticity/flutter_solution.h"
#include "Flight/flight_condition.h"
#include "Aeroelasticity/coupled_aero_structural_model.h"

// libMesh includes
#include "libmesh/parameter_vector.h"


namespace MAST {
        
    
    class FlutterSolverBase
    {
    public:
        
        FlightCondition* flight_condition;
        
        CoupledAeroStructuralModel* aero_structural_model;
        
        /*!
         *    file to which the result will be written
         */
        std::ofstream _output, _mode_output;
        
        /*!
         *    constructor for the flutter solver base object
         */
        FlutterSolverBase()
        {}
        
        virtual ~FlutterSolverBase();
        
        
        virtual void clear_solutions();
        
        /*!
         *    finds the number of critical points already identified in the
         *    procedure.
         */
        virtual unsigned int n_roots_found() const;
        
        
        /*!
         *   returns the \par n th root in terms of ascending velocity that is
         *   found by the solver
         */
        const MAST::FlutterRootBase& get_root(const unsigned int n) const;
        
        

        void set_output_file(std::string& nm)
        {
            _output.close();
            _output.open(nm.c_str(), std::ofstream::out);
            std::ostringstream oss;
            oss << "modes_" << nm;
            _mode_output.open(oss.str().c_str(), std::ofstream::out);
        }
        
        
        /*!
         *   Looks through the list of flutter cross-over points and iteratively
         *   zooms in to find the cross-over point. This should be called only
         *   after scan_for_roots() has been called. Potential cross-over points
         *   are sorted with increasing velocity, and this method will attempt to
         *   identify the next critical root in the order.
         */
        virtual std::pair<bool, const MAST::FlutterRootBase*>
        find_next_root(const Real g_tol,
                       const unsigned int n_bisection_iters);
        
        
        /*!
         *   This method checks if the flutter root corresponding to the
         *   lowest velocity crossover has been calculated. If not, then it
         *   attempts to find that root using an iterative approach
         */
        virtual std::pair<bool, const MAST::FlutterRootBase*>
        find_critical_root(const Real g_tol,
                           const unsigned int n_bisection_iters);

        /*!
         *   Calculate the sensitivity of the flutter root with respect to the
         *   \par i^th parameter in params
         */
        virtual void calculate_sensitivity(MAST::FlutterRootBase& root,
                                           const libMesh::ParameterVector& params,
                                           const unsigned int i) = 0;
        
        /*!
         *   Prints the sorted roots to the \par output
         */
        void print_sorted_roots(std::ostream* output = NULL);
        
        
        /*!
         *   Prints the crossover points output. If no pointer to output is given
         *   then the output defined by set_output_file() is used.
         */
        void print_crossover_points(std::ostream* output = NULL);
        
        /*!
         *   Scans for flutter roots in the range specified, and identified the
         *   divergence (if k_red = 0. is specified) and flutter crossover points.
         *   The roots are organized in terms of increasing velocity.
         */
        virtual void scan_for_roots() = 0;

        
        /*!
         *    creates a new flutter root and returns pointer to it.
         */
        virtual std::auto_ptr<MAST::FlutterRootBase> build_flutter_root() = 0;
        
    protected:
        

        /*!
         *    bisection method search
         */
        std::pair<bool, MAST::FlutterSolutionBase*>
        bisection_search(const std::pair<MAST::FlutterSolutionBase*,
                         MAST::FlutterSolutionBase*>& ref_sol_range,
                         const unsigned int root_num,
                         const Real g_tol,
                         const unsigned int max_iters);
        
        /*!
         *   performs an eigensolution at the specified reference value, and
         *   sort the roots based on the provided solution pointer. If the
         *   pointer is NULL, then no sorting is performed
         */
        virtual MAST::FlutterSolutionBase*
        analyze(const Real k_red,
                const Real v_ref,
                const MAST::FlutterSolutionBase* prev_sol=NULL) = 0;


        /*!
         *   identifies all cross-over and divergence points from analyzed
         *   roots
         */
        virtual void _identify_crossover_points() = 0;
        
        /*!
         *   map of reduced frequency vs flutter solutions
         */
        std::map<Real, MAST::FlutterSolutionBase*> _flutter_solutions;
        
        /*!
         *   the map of flutter crossover points versus average velocity of the
         *   two bounding roots
         */
        std::multimap<Real, MAST::FlutterRootCrossoverBase*> _flutter_crossovers;
        
    };
    
}

#endif /* defined(__MAST__flutter_solver_base__) */
