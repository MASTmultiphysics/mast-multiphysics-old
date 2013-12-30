//
//  optimization_interface.h
//  MAST
//
//  Created by Manav Bhatia on 11/12/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_optimization_interface_h__
#define __MAST_optimization_interface_h__

// C++ includes
#include <vector>
#include <iostream>
#include <iomanip>


// MAST includes
#include "Base/MAST_data_types.h"


namespace MAST {

    class FunctionEvaluation {
        
    public:
        
        FunctionEvaluation():
        _n_vars(0),
        _n_eq(0),
        _n_ineq(0),
        _max_iters(0),
        _n_rel_change_iters(5),
        _tol(1.0e-6)
        { }
        
        virtual ~FunctionEvaluation() { }
        
        
        unsigned int n_vars() const {
            return _n_vars;
        }
        
        
        unsigned int n_eq() const {
            return _n_eq;
        }
        
        
        unsigned int n_ineq() const{
            return _n_ineq;
        }
        
        
        unsigned int max_iters() const{
            return _max_iters;
        }
        
        
        unsigned int n_iters_relative_change() const {
            return _n_rel_change_iters;
        }
        
        
        Real tolerance() const{
            return _tol;
        }
        
        
        virtual void init_dvar(std::vector<Real>& x,
                               std::vector<Real>& xmin,
                               std::vector<Real>& xmax) = 0;
        
        /*!
         *   \par grads(k): Derivative of f_i(x) with respect 
         *   to x_j, where k = (j-1)*M + i. 
         */
        virtual void evaluate(const std::vector<Real>& dvars,
                              Real& obj,
                              bool eval_obj_grad,
                              std::vector<Real>& obj_grad,
                              std::vector<Real>& fvals,
                              std::vector<bool>& eval_grads,
                              std::vector<Real>& grads) = 0;
        
        virtual void output(unsigned int iter,
                            const std::vector<Real>& x,
                            Real obj,
                            const std::vector<Real>& fval) const;

        virtual void tabulated_output(std::ostream& output,
                                      unsigned int iter,
                                      const std::vector<Real>& x,
                                      Real obj,
                                      const std::vector<Real>& fval) const;
    protected:
        
        unsigned int _n_vars;
        
        unsigned int _n_eq;
        
        unsigned int _n_ineq;
        
        unsigned int _max_iters;
        
        unsigned int _n_rel_change_iters;
        
        Real _tol;
    };
    

    
    class OptimizationInterface {
    public:
     
        OptimizationInterface():
        _feval(NULL)
        { }
        
        virtual ~OptimizationInterface()
        { }

        
        virtual void optimize() = 0;
        

        void attach_function_evaluation_object
        (MAST::FunctionEvaluation& feval)
        {
            _feval = &feval;
        }
        
    protected:
        
        MAST::FunctionEvaluation* _feval;
    };
}


inline void
MAST::FunctionEvaluation::output(unsigned int iter, const std::vector<Real> &x,
                                 Real obj, const std::vector<Real> &fval) const {

    libmesh_assert_equal_to(x.size(), _n_vars);
    libmesh_assert_equal_to(fval.size(), _n_eq + _n_ineq);
    
    
    libMesh::out
    << " *************************** " << std::endl
    << " *** Optimization Output *** " << std::endl
    << " *************************** " << std::endl
    << std::endl
    << "Iter:            "  << std::setw(10) << iter << std::endl
    << "Nvars:           " << std::setw(10) << x.size() << std::endl
    << "Ncons-Equality:  " << std::setw(10) << _n_eq << std::endl
    << "Ncons-Inquality: " << std::setw(10) << _n_ineq << std::endl
    << std::endl
    << "Obj =                  " << std::setw(20) << obj << std::endl
    << std::endl
    << "Vars:            " << std::endl;
    
    for (unsigned int i=0; i<_n_vars; i++)
        libMesh::out
        << "x     [ " << std::setw(10) << i << " ] = "
        << std::setw(20) << x[i] << std::endl;
    
    if (_n_eq) {
        libMesh::out << std::endl
        << "Equality Constraints: " << std::endl;
        
        for (unsigned int i=0; i<_n_eq; i++)
            libMesh::out
            << "feq [ " << std::setw(10) << i << " ] = "
            << std::setw(20) << fval[i] << std::endl;
    }
    
    if (_n_ineq) {
        libMesh::out << std::endl
        << "Inequality Constraints: " << std::endl;
        
        for (unsigned int i=0; i<_n_ineq; i++)
            libMesh::out
            << "fineq [ " << std::setw(10) << i << " ] = "
            << std::setw(20) << fval[i+_n_eq] << std::endl;
    }
    
    libMesh::out << std::endl
    << " *************************** " << std::endl;
}



inline void
MAST::FunctionEvaluation::tabulated_output(std::ostream& output,
                                           unsigned int iter,
                                           const std::vector<Real>& x,
                                           Real obj,
                                           const std::vector<Real>& fval) const {
    libmesh_assert_equal_to(x.size(), _n_vars);
    libmesh_assert_equal_to(fval.size(), _n_eq + _n_ineq);
    
    // write header for the first iteration
    if (iter == 0) {
        output << std::setw(10) << "Iter";
        for (unsigned int i=0; i < x.size(); i++)
            output << " x [ " << std::setw(10) << i << " ]   ";
        output << std::setw(20) << "Obj";
        for (unsigned int i=0; i<fval.size(); i++)
            output << " f [ " << std::setw(10) << i << " ]   ";
    }
    
    output << std::setw(10) << iter;
    for (unsigned int i=0; i < x.size(); i++)
        output << std::setw(20) << x[i];
    output << std::setw(20) << obj;
    for (unsigned int i=0; i < fval.size(); i++)
        output << std::setw(20) << fval[i];
}


#endif  // __MAST_optimization_interface_h__
