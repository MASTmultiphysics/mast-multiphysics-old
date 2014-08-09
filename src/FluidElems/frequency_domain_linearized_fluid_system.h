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

#ifndef __MAST__frequency_domain_linearized_euler__
#define __MAST__frequency_domain_linearized_euler__

// libMesh includes
#include "libmesh/libmesh_config.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fem_system.h"
#include "libmesh/auto_ptr.h"
#include "libmesh/mesh_function.h"

// MAST includes
#include "FluidElems/fluid_elem_base.h"




class SurfaceMotionBase;
class FlightCondition;

class FrequencyDomainLinearizedFluidSystem: public libMesh::FEMSystem, public FluidElemBase
{
public:
    FrequencyDomainLinearizedFluidSystem(libMesh::EquationSystems& es,
                                         const std::string& name_in,
                                         const unsigned int number_in):
    libMesh::FEMSystem(es, name_in, number_in),
    FluidElemBase(*es.parameters.get<GetPot*>("input_file")),
    perturbed_surface_motion(NULL),
    if_k_red_sensitivity(false),
    if_Vref_sensitivity(false),
    _if_localized_sol(false)
    { }
    
    void init_data();
    
    virtual void init_context(libMesh::DiffContext &context);
    
    
    void localize_fluid_solution();
    
    
    virtual bool element_time_derivative (bool request_jacobian,
                                          libMesh::DiffContext &context);

    virtual bool element_time_derivative_sens (bool request_jacobian,
                                               libMesh::DiffContext &context);

    virtual bool side_time_derivative (bool request_jacobian,
                                       libMesh::DiffContext &context);

    virtual bool side_time_derivative_sens (bool request_jacobian,
                                            libMesh::DiffContext &context);

    std::vector<unsigned int> vars;

    /*!
     *   this defines the small disturbance surface motion on top of the 
     *   steady surface motion that the body might have seen.
     */
    MAST::SurfaceMotionBase* perturbed_surface_motion;
    
    /*!
     *   flag to tell the system if the quantity being solved for is sensitivity
     *    wrt reduced_frequency.
     */
    bool if_k_red_sensitivity;

    
    /*!
     *   flag to tell the system if the quantity being solved for is sensitivity
     *    wrt reference velocity.
     */
    bool if_Vref_sensitivity;
    
    /*!
     *   sets the base solution for sensitivity analysis
     */
    void set_base_solution(libMesh::NumericVector<Real>& vec);
    
protected:
    
    bool _if_localized_sol;
    
    std::auto_ptr<libMesh::NumericVector<Real> > _local_fluid_solution;
    
    /*!
     *   this is the base solution about which sensitivity is to be calculated
     */
    std::auto_ptr<libMesh::NumericVector<Real> > _base_fluid_solution;
    
};




class FrequencyDomainFluidPostProcessSystem : public libMesh::System
{
public:
    // Constructor
    FrequencyDomainFluidPostProcessSystem(libMesh::EquationSystems& es,
                                          const std::string& name_in,
                                          const unsigned int number_in)
    : libMesh::System(es, name_in, number_in)
    {}
    
    virtual void init_data();
    
    virtual void postprocess();
    
    unsigned int u, v, w, T, s, p, cp, a, M, // real part of the variables
                                             // imaginary part of the variables
    rho_im, u_im, v_im, w_im, T_im, s_im, p_im, cp_im, a_im, M_im;

    FlightCondition* flight_condition;
};




#endif /* defined(__MAST__frequency_domain_linearized_euler__) */
