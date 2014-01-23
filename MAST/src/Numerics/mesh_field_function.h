//
//  mesh_field_function.h
//  MAST
//
//  Created by Manav Bhatia on 1/22/14.
//  Copyright (c) 2014 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_mesh_field_function_h__
#define __MAST_mesh_field_function_h__

// C++ includes
#include <vector>
#include <string>

// MAST includes
#include "Numerics/function_base.h"

// libMesh includes
#include "libmesh/mesh_function.h"


namespace MAST {
    
    template <typename ValType>
    class MeshFieldFunction: public MAST::FieldFunction<ValType> {
        
    public:
        MeshFieldFunction(const std::string& nm,
                          System& sys,
                          const std::vector<string>& vars):
        MAST::FieldFunction<ValType>(nm),
        _sys(sys),
        _sol(NULL),
        _sens(NULL),
        _function(NULL),
        _function_total_derivative(NULL) {
            // get the variable numbers
            std::vector<unsigned int> var_id(vars.size());
            for (unsigned int i=0; i<vars.size(); i++)
                var_id[i] = sys.variable_number(vars[i]);
            
            // prepare the localized solution
            _sol = NumericVector<ValType>::build(sys.get_equation_systems().comm()).release();
            _sol->init(sys.solution->size(), true, SERIAL);
            sys.solution->localize(*_sol, sys.get_dof_map().get_send_list());
            
            _function = new MeshFunction(sys,
                                         *_sol,
                                         sys.get_dof_map(),
                                         var_id);
            _function->init();
        }
        

        virtual ~MeshFieldFunction() {
            
            // delete the objects before exiting
            if (_function) {
                delete _function;
                delete _sol;
            }
            
            if (_function_total_derivative) {
                delete _function_total_derivative;
                delete _sens;
            }
        }
        
        
        
        
        /*!
         *    initialize the data structure for sensitivity analysis
         */
        void init_for_sens(const unsigned int pid) {
            // prepare the localized solution
            if (!_function_total_derivative) {
                _sens = NumericVector<ValType>::build(_sys.get_equation_systems().comm()).release();
                _function_total_derivative = new MeshFunction(_sys,
                                                              *_sens,
                                                              _sys.get_dof_map(),
                                                              var_id);
                _function_total_derivative->init();
            }

            _sens->init(_sys.solution->size(), true, SERIAL);
            _sys.get_sensitivity_solution(pid).localize(*_sol, _sys.get_dof_map().get_send_list());
        }
        
        
        
        /*!
         *    Returns the value of this function.
         */
        virtual void operator() (const Point& p, const Real t, ValType& v) const {
            (*_function)(p, t, v);
        }
        
        
        
        /*!
         *    Returns the partial derivative of this function with respect to
         *    the sensitivity parameter \par p
         */
        virtual void partial_derivative (const MAST::SensitivityParameters& par,
                                         const Point& p, const Real t, ValType& v) const {
            // there are no partial derivatives here. A system only calculates
            // total derivatives
            libMesh::err
            << "Error: MeshFieldFunction does not provide partial derivative. "
            << std::endl << "Exiting." << std::endl;
            libmesh_error(false);
        }
        
        
        
        /*!
         *    Returns the total derivative of this function with respect to
         *    the sensitivity parameter \par p
         */
        virtual void total_derivative (const MAST::SensitivityParameters& par,
                                          const Point& p, const Real t, ValType& v) const {
            // make sure that this has been initialized
            libmesh_assert(_function_total_derivative);
            (*_function_total_derivative)(p, t, v);
        }
        
        
    protected:

        /*!
         *   reference to system for which this object provides solution
         */
        System& _sys;
        
        /*!
         *   localized vectors
         */
        NumericVector<ValType> *_sol, *_sens;
        
        
        /*!
         *   mesh function object that provides the interpolation
         */
        MeshFunction *_function;
        
        /*!
         *   total derivative mesh function
         */
        MeshFunction *_function_total_derivative;
    };
}



#endif  // __MAST_mesh_field_function_h__
