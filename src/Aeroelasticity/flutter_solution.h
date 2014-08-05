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


#ifndef __MAST_flutter_solution_h__
#define __MAST_flutter_solution_h__

// C++ includes
#include <vector>


// MAST includes
#include "Solvers/lapack_interface.h"


namespace MAST {
    
    // forward declerations
    class FlutterSolverBase;
    
    class FlutterRootBase
    {
    public:
        FlutterRootBase():
        if_nonphysical_root(false),
        has_sensitivity_data(false),
        k_red_ref(0.),
        V_ref(0.),
        V(0.),
        g(0.),
        omega(0.),
        k_red(0.),
        V_sens(0.),
        k_red_sens(0.),
        root(0.),
        root_sens(0.)
        {}
        
        /*!
         *   copy constructor
         */
        FlutterRootBase(const FlutterRootBase& f):
        if_nonphysical_root(f.if_nonphysical_root),
        has_sensitivity_data(f.has_sensitivity_data),
        k_red_ref(f.k_red_ref),
        V_ref(f.V_ref),
        V(f.V),
        g(f.g),
        omega(f.omega),
        k_red(f.k_red),
        V_sens(f.V_sens),
        k_red_sens(f.k_red_sens),
        root(f.root),
        root_sens(f.root_sens),
        eig_vec_right(f.eig_vec_right),
        eig_vec_left(f.eig_vec_left),
        modal_participation(f.modal_participation)
        {}
        
        
        
        virtual ~FlutterRootBase() {}
        
        bool if_nonphysical_root;
        
        bool has_sensitivity_data;
        
        Real k_red_ref, V_ref, V, g, omega, k_red, V_sens;
        
        Complex root, root_sens, k_red_sens;
        
        /*!
         *    right and left eigenvevtors
         */
        ComplexVectorX  eig_vec_right, eig_vec_left;
        
        RealVectorX modal_participation;
        
    };

    
    class FrequencyDomainFlutterRoot: public MAST::FlutterRootBase {
    public:
        FrequencyDomainFlutterRoot():
        MAST::FlutterRootBase()
        { }
        
        virtual ~FrequencyDomainFlutterRoot() {}

        virtual void init(const Real k_red_ref,
                          const Real v_ref,
                          const Real b_ref,
                          const Complex num,
                          const Complex den,
                          const ComplexMatrixX& Bmat,
                          const ComplexVectorX& evec_right,
                          const ComplexVectorX& evec_left) = 0;
    };

    


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
    
    
    
    class FlutterSolutionBase
    {
    public:
        FlutterSolutionBase():
        _ref_val(0.)
        {}
        
        /*!
         *    delete the flutter root objects
         */
        virtual ~FlutterSolutionBase() {
            std::vector<MAST::FlutterRootBase*>::iterator
            it = _roots.begin();
            for ( ; it != _roots.end(); it++)
                delete *it;
        }

        
        /*!
         *   the reduced frequency for this solution
         */
        Real ref_val() const {
            return _ref_val;
        }
        
        /*!
         *   number of roots in this solution
         */
        unsigned int n_roots() const{
            return (unsigned int)_roots.size();
        }
        
        /*!
         *    returns the root
         */
        const MAST::FlutterRootBase& get_root(const unsigned int i) const {
            libmesh_assert_less(i, _roots.size());
            return *_roots[i];
        }
        
        /*!
         *    sort this root with respect to the given solution from a previous
         *    eigen solution. This method relies on the modal participation.
         *    Flutter roots from previous and current solutions with highest
         *    dot product of modal participation vector are considered to be
         *    similar.
         */
        void sort(const MAST::FlutterSolutionBase& sol);
        
        /*!
         *
         */
        void mark_inconsistent_roots_as_nonphysical(const Real ref_val,
                                                    const Real tol);

        /*!
         *
         */
        void swap_root(MAST::FlutterSolutionBase& sol,
                       unsigned int root_num);
        
        /*!
         *    prints the data and modes from this solution
         */
        void print(std::ostream& output, std::ostream& mode_output);
        
    protected:
                
        /*!
         *    Reference value of the sweeping parameter for which this solution 
         *    was obtained. For UG solver, this is k_red, and for time domain
         *    solver this could be velocity. PK solver will need additional 
         *    reference values, provided in the inherited class.
         */
        Real _ref_val;
        
        /*!
         *    Matrix used for scaling of eigenvectors, and sorting of roots
         */
        ComplexMatrixX _Bmat;
        
        std::vector<MAST::FlutterRootBase*> _roots;
    };
    
    
    
    
    class FrequencyDomainFlutterSolution: public MAST::FlutterSolutionBase {
    public:
        
        FrequencyDomainFlutterSolution():
        MAST::FlutterSolutionBase()
        { }
        
        
        virtual ~FrequencyDomainFlutterSolution() {}
        
        /*!
         *   initializes the flutter solution from an eigensolution
         */
        virtual void init (const MAST::FlutterSolverBase& solver,
                           const Real k_red,
                           const Real v_ref,
                           const Real bref,
                           const LAPACK_ZGGEV& eig_sol);
        
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

    
    
    
    class FlutterRootCrossoverBase
    {
    public:
        FlutterRootCrossoverBase():
        crossover_solutions(NULL, NULL),
        root_num(0),
        root(NULL)
        { }
        
        void print(std::ostream& output) const;
        
        std::pair<MAST::FlutterSolutionBase*, MAST::FlutterSolutionBase*>
        crossover_solutions;
        
        unsigned int root_num;
        
        const MAST::FlutterRootBase* root;
    };
    
    
    

}


#endif // __MAST_flutter_solution_h__
