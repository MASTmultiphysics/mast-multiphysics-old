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

#ifndef MAST_surface_pressure_load_h
#define MAST_surface_pressure_load_h

// libMesh includes
#include "libmesh/mesh_function.h"
#include "libmesh/system.h"
#include "libmesh/dof_map.h"
#include "libmesh/equation_systems.h"
#include "libmesh/mesh_serializer.h"
#include "libmesh/numeric_vector.h"


// MAST includes
#include "BoundaryConditions/boundary_condition.h"
#include "Flight/flight_condition.h"
#include "FluidElems/frequency_domain_linearized_fluid_system.h"
#include "FluidElems/fluid_elem_base.h"


namespace MAST {
    
    class SurfacePressure: public MAST::BoundaryCondition {
    public:
        SurfacePressure():
        MAST::BoundaryCondition(MAST::SURFACE_PRESSURE)
        { }
        
        virtual ~SurfacePressure()
        { }
        
    };
    
    
    
    class SmallDisturbanceSurfacePressure: public MAST::BoundaryCondition
    {
    public:
        SmallDisturbanceSurfacePressure(libMesh::System& nl_sys,
                                        libMesh::System& lin_sys):
        MAST::BoundaryCondition(MAST::SMALL_DISTURBANCE_PRESSURE),
        nonlinear_sys(nl_sys),
        linearized_sys(lin_sys),
        _flt_cond(NULL),
        _dim(0)
        {
            
            libMesh::MeshBase& linear_sys_mesh = linearized_sys.get_mesh();
            _linear_mesh_serializer.reset(new libMesh::MeshSerializer(linear_sys_mesh, true));
            
            libMesh::MeshBase& nonlinear_sys_mesh = nonlinear_sys.get_mesh();
            _nonlinear_mesh_serializer.reset(new libMesh::MeshSerializer(nonlinear_sys_mesh, true));
            
            
            _dim = nonlinear_sys.n_vars()-2;
            
            // copy the pointer for flight condition data
            _flt_cond = dynamic_cast<FrequencyDomainLinearizedFluidSystem&>
            (lin_sys).flight_condition;
        }
        
        virtual ~SmallDisturbanceSurfacePressure()
        { }
        
        virtual void init(libMesh::NumericVector<Real>& nonlinear_sol,
                          libMesh::NumericVector<Real>& linearized_sol);
        
        // calculation in frequency domain
        template <typename ValType>
        void surface_pressure(const Real t,
                              const libMesh::Point& p,
                              Real& cp,
                              ValType& dcp);
        
    protected:
        
        /*!
         *   systems that this is attacehd to
         */
        libMesh::System& nonlinear_sys;
        libMesh::System& linearized_sys;
        
        /*!
         *   mesh function that interpolates the nonlinear solution
         */
        std::auto_ptr<libMesh::MeshFunction> _function_nonlinear;
        
        /*!
         *   mesh function that interpolates the linearized solution
         */
        std::auto_ptr<libMesh::MeshFunction> _function_linear;
        
        /*!
         *    numeric vector that stores the solution for nonlinear system
         */
        std::auto_ptr<libMesh::NumericVector<Real> > _sol_nonlinear;
        
        /*!
         *    numeric vector that stores the solution for linearized system
         */
        std::auto_ptr<libMesh::NumericVector<Real> > _sol_linear;
        
        /*!
         *    this provides the fluid values for calculation of cp
         */
        FlightCondition* _flt_cond;
        
        /*!
         *    dimension of the analysis mesh
         */
        unsigned int _dim;
        
        /*!
         *   this serializes the mesh for use in interpolation
         */
        std::auto_ptr<libMesh::MeshSerializer> _nonlinear_mesh_serializer;
        std::auto_ptr<libMesh::MeshSerializer> _linear_mesh_serializer;
        
    };
}



#endif
