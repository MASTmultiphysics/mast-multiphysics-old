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
                          const std::vector<string>& vars) {
            // get the variable numbers
            std::vector<unsigned int> var_id(vars.size());
            for (unsigned int i=0; i<vars.size(); i++)
                var_id[i] = sys.variable_number(vars[i]);
            
            _mesh_function = new MeshFunction(sys,
                                              *sys.solution,
                                              sys.get_dof_map(),
                                              var_id);
            _mesh_function->init();
        }
        

        ~MeshFieldFunction() {
            delete _mesh_function;
        }
        
        
        /*!
         *    initialize the data structure for sensitivity analysis
         */
        void initialize_for_sensitivity(const unsigned int pid,
                                        const MAST::SensitivityParameters& par);
        
        
        /*!
         *    Returns the value of this function.
         */
        virtual void operator() (const Point& p, const Real t, ValType& v) const {
            (*_mesh_function)(p, t, v);
        }
        
        /*!
         *    Returns the partial derivative of this function with respect to
         *    the sensitivity parameter \par p
         */
        virtual ValType partial_derivative (const MAST::SensitivityParameters& par,
                                            const Point& p, const Real t, ValType& v) const {
            // make sure that the object has been initialized for sensitivity analysis
            
        }
        
        /*!
         *    Returns the total derivative of this function with respect to
         *    the sensitivity parameter \par p
         */
        virtual ValType total_derivative (const MAST::SensitivityParameters& par,
                                          const Point& p, const Real t, ValType& v) const;
        
        
    protected:

        /*!
         *   mesh function object that provides the interpolation
         */
        MeshFunction *_function;
        
        /*!
         *   partial derivative mesh function
         */
        MeshFunction *_function_partial_derivative;

        /*!
         *   total derivative mesh function
         */
        MeshFunction *_function_total_derivative;
    };
}



#endif  // __MAST_mesh_field_function_h__
