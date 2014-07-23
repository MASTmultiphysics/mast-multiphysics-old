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

#ifndef MAST_cfd_aerodynamic_model_h
#define MAST_cfd_aerodynamic_model_h

// C++ includes
#include <memory>

// MAST includes
#include "Aeroelasticity/aerodynamic_model.h"
#include "FluidElems/frequency_domain_linearized_fluid_system.h"
#include "BoundaryConditions/flexible_surface_motion.h"


class CFDAerodynamicModel: public AerodynamicModel
{
public:
    CFDAerodynamicModel(libMesh::System& nl_sys,
                        FrequencyDomainLinearizedFluidSystem& lin_sys):
    AerodynamicModel(),
    nonlinear_fluid_system(nl_sys),
    linearized_fluid_system(lin_sys)
    { }
    
    virtual ~CFDAerodynamicModel()
    { }
    
    /*!
     *    updates the aerodynamic damping matrix operator in \par a. This is
     *    applicable for only the time-domain methods. Returns \par false if
     *    the model does not have a damping term.
     */
    virtual bool get_damping_matrix(RealMatrixX& a)
    { libmesh_assert(false); }
    
    /*!
     *    updates the aerodynamic stiffness matrix operator in \par a. This is
     *    applicable for only the time-domain methods. Returns \par false if
     *    the model does not have a stiffness term.
     */
    virtual bool get_stiffness_matrix(RealMatrixX& a)
    { libmesh_assert(false); }
    

    /*!
     *   Nonlinear fluid system
     */
    libMesh::System& nonlinear_fluid_system;

    /*!
     *   The small-disturbance fluid system object that provides the basis
     *   of calculations for this model.
     */
    FrequencyDomainLinearizedFluidSystem& linearized_fluid_system;
    
};





#endif
