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

#ifndef __MAST_linearized_fluid_system_h__
#define __MAST_linearized_fluid_system_h__

// libMesh includes
#include "libmesh/libmesh_config.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fem_system.h"
#include "libmesh/auto_ptr.h"
#include "libmesh/mesh_function.h"

// MAST includes
#include "FluidElems/fluid_elem_base.h"


class LinearizedFluidSystem : public libMesh::FEMSystem, public FluidElemBase
{
public:
    // Constructor
    LinearizedFluidSystem(libMesh::EquationSystems& es,
                const std::string& name_in,
                const unsigned int number_in)
    : libMesh::FEMSystem(es, name_in, number_in),
    FluidElemBase(*es.parameters.get<GetPot*>("input_file")),
    perturbed_surface_motion(NULL),
    _rho_norm_old(1.),
    _rho_norm_curr(1.),
    dc_recalculate_tolerance(1.0e-8),
    if_use_stored_dc_coeff(false)
    {}
    
    void init_data();
    
    virtual void init_context(libMesh::DiffContext &context);
    
    void localize_fluid_solution();

    virtual bool element_time_derivative (bool request_jacobian,
                                          libMesh::DiffContext &context);
    
    virtual bool side_time_derivative (bool request_jacobian,
                                       libMesh::DiffContext &context);
    
    virtual bool mass_residual (bool request_jacobian,
                                libMesh::DiffContext& context);
    
    virtual void postprocess();
    
    
    void evaluate_recalculate_dc_flag();
    
    /*!
     *    tolerance threshold for turning on/off the stored dc coeff flag
     */
    Real dc_recalculate_tolerance;
    
    /*!
     *     flag to tell the system to use the stored discontinuity capturing terms
     */
    bool if_use_stored_dc_coeff;
    
    std::vector<unsigned int> vars;
    
    /*!
     *   this defines the small disturbance surface motion on top of the
     *   steady surface motion that the body might have seen.
     */
    MAST::SurfaceMotionBase* perturbed_surface_motion;

protected:
    
    /*!
     *   localized reference solution for calculation of dc_val
     */
    std::auto_ptr<libMesh::NumericVector<Real> > _dc_ref_sol;
    
    /*!
     *    Current and old norms of density in the flow-field
     */
    Real _rho_norm_old, _rho_norm_curr;
    
    
    bool _if_localized_sol;
    
    libMesh::AutoPtr<libMesh::NumericVector<Real> > _local_fluid_solution;
};



#endif // __MAST_linearized_fluid_system_h__
