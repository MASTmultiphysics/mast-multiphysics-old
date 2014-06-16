//
//  unsteady_compressible_potential_flow_elem.h
//  MAST
//
//  Created by Manav Bhatia on 10/7/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST__unsteady_compressible_potential_flow_elem__
#define __MAST__unsteady_compressible_potential_flow_elem__

// MAST includes
#include "PotentialFlowElems/potential_flow_elem_base.h"

// libMesh includes
#include "libmesh/fem_system.h"
#include "libmesh/equation_systems.h"

libMesh::Real unsteady_compressible_potential_solution_value(const libMesh::Point& p,
                                                    const Parameters& parameters,
                                                    const std::string& sys_name,
                                                    const std::string& var_name);



void init_compressible_potential_variables(libMesh::EquationSystems& es,
                                           const std::string& system_name);


/*!
 *    implements the equations for unsteady compressible Potential flow
 */

class UnsteadyCompressiblePotentialFlow :
public FEMSystem, public MAST::PotentialFlowElemBase
{
public:
    UnsteadyCompressiblePotentialFlow(libMesh::EquationSystems& es,
                                      const std::string& name_in,
                                      const unsigned int number_in):
    FEMSystem(es, name_in, number_in),
    PotentialFlowElemBase (*es.parameters.get<GetPot*>("input_file"))
    { }

    
    virtual ~UnsteadyCompressiblePotentialFlow()
    { }

    void init_data();
    
    virtual void init_context(DiffContext &context);

    
    virtual bool element_time_derivative (bool request_jacobian,
                                          DiffContext &context);
    
    virtual bool side_time_derivative (bool request_jacobian,
                                       DiffContext &context);
    
    virtual bool mass_residual (bool request_jacobian,
                                DiffContext& context);
    
protected:
    
    /*!
     *   finite element variables
     */
    std::vector<unsigned int> vars;
  
    
    void update_solution_at_quadrature_point
    ( const std::vector<unsigned int>& vars, const unsigned int qp, FEMContext& c,
     const bool if_elem_domain, const DenseRealVector& elem_solution,
     DenseRealVector& conservative_sol, libMesh::Point& uvec,
     FEMOperatorMatrix& B_mat, std::vector<FEMOperatorMatrix>& dB_mat,
     DenseRealMatrix& LS_mat);
};


#endif /* defined(__MAST__unsteady_compressible_potential_flow_elem__) */
