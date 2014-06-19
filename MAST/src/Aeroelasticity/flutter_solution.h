//
//  flutter_solution.h
//  MAST
//
//  Created by Manav Bhatia on 4/15/14.
//  Copyright (c) 2014 Manav Bhatia. All rights reserved.
//


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
        V(0.),
        g(0.),
        omega(0.),
        k_ref(0.),
        root(0.)
        {}
        
        virtual ~FlutterRootBase() {}
        
        bool if_nonphysical_root;
        
        Real V, g, omega, k_ref;
        
        Complex root;
        
        ComplexVectorX  mode;
        
        RealVectorX modal_participation;
        
    };

    
    class FrequencyDomainFlutterRoot: public MAST::FlutterRootBase {
    public:
        FrequencyDomainFlutterRoot():
        MAST::FlutterRootBase()
        { }
        
        virtual ~FrequencyDomainFlutterRoot() {}

        virtual void init(const Real ref_val, const Real b_ref,
                          const Complex num, const Complex den,
                          const ComplexMatrixX& Bmat,
                          const ComplexVectorX& eig_vec) = 0;
    };

    


    class TimeDomainFlutterRoot: public MAST::FlutterRootBase {
    public:
        TimeDomainFlutterRoot():
        MAST::FlutterRootBase()
        { }
        
        virtual ~TimeDomainFlutterRoot() {}
        
        virtual void init(const Real ref_val, const Real b_ref,
                          const Complex num, const Complex den,
                          const RealMatrixX& Bmat,
                          const ComplexVectorX& eig_vec);
    };
    
    
    
    class FlutterSolutionBase
    {
    public:
        FlutterSolutionBase(MAST::FlutterSolverBase& solver):
        _ref_val(0.),
        _solver(solver)
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
         *    prints the data and modes from this solution
         */
        void print(std::ostream& output, std::ostream& mode_output);
        
    protected:
                
        /*!
         *    Reference value of the sweeping parameter for which this solution 
         *    was obtained. For UG solver, this is k_ref, and for time domain
         *    solver this could be velocity. PK solver will need additional 
         *    reference values, provided in the inherited class.
         */
        Real _ref_val;
        
        /*!
         *   solver for which this is initialized
         */
        MAST::FlutterSolverBase& _solver;
        
        std::vector<MAST::FlutterRootBase*> _roots;
    };
    
    
    
    
    class FrequencyDomainFlutterSolution: public MAST::FlutterSolutionBase {
    public:
        
        FrequencyDomainFlutterSolution(MAST::FlutterSolverBase& solver):
        MAST::FlutterSolutionBase(solver)
        { }
        
        
        virtual ~FrequencyDomainFlutterSolution() {}
        
        /*!
         *   initializes the flutter solution from an eigensolution
         */
        virtual void init (const Real ref_val, const Real bref,
                           const LAPACK_ZGGEV& eig_sol);
        
    };

    
    
    
    class TimeDomainFlutterSolution: public MAST::FlutterSolutionBase {
    public:
        
        TimeDomainFlutterSolution(MAST::FlutterSolverBase& solver):
        MAST::FlutterSolutionBase(solver)
        { }
        
        virtual ~TimeDomainFlutterSolution() {}

        /*!
         *   initializes the flutter solution from an eigensolution
         */
        virtual void init (const Real ref_val, const Real bref,
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
