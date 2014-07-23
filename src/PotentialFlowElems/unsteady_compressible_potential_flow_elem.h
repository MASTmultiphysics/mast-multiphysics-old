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

#ifndef __MAST__unsteady_compressible_potential_flow_elem__
#define __MAST__unsteady_compressible_potential_flow_elem__

// MAST includes
#include "PotentialFlowElems/potential_flow_elem_base.h"

// libMesh includes
#include "libmesh/fem_system.h"
#include "libmesh/equation_systems.h"

Real unsteady_compressible_potential_solution_value(const libMesh::Point& p,
                                                    const libMesh::Parameters& parameters,
                                                    const std::string& sys_name,
                                                    const std::string& var_name);



void init_compressible_potential_variables(libMesh::EquationSystems& es,
                                           const std::string& system_name);


/*!
 *    implements the equations for unsteady compressible Potential flow
 */

class UnsteadyCompressiblePotentialFlow :
public libMesh::FEMSystem, public MAST::PotentialFlowElemBase
{
public:
    UnsteadyCompressiblePotentialFlow(libMesh::EquationSystems& es,
                                      const std::string& name_in,
                                      const unsigned int number_in):
    libMesh::FEMSystem(es, name_in, number_in),
    PotentialFlowElemBase (*es.parameters.get<GetPot*>("input_file"))
    { }

    
    virtual ~UnsteadyCompressiblePotentialFlow()
    { }

    void init_data();
    
    virtual void init_context(libMesh::DiffContext &context);

    
    virtual bool element_time_derivative (bool request_jacobian,
                                          libMesh::DiffContext &context);
    
    virtual bool side_time_derivative (bool request_jacobian,
                                       libMesh::DiffContext &context);
    
    virtual bool mass_residual (bool request_jacobian,
                                libMesh::DiffContext& context);
    
protected:
    
    /*!
     *   finite element variables
     */
    std::vector<unsigned int> vars;
  
    
    void update_solution_at_quadrature_point
    ( const std::vector<unsigned int>& vars, const unsigned int qp, libMesh::FEMContext& c,
     const bool if_elem_domain, const DenseRealVector& elem_solution,
     DenseRealVector& conservative_sol, libMesh::Point& uvec,
     FEMOperatorMatrix& B_mat, std::vector<FEMOperatorMatrix>& dB_mat,
     DenseRealMatrix& LS_mat);
};


#endif /* defined(__MAST__unsteady_compressible_potential_flow_elem__) */
