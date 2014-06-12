//
//  surface_pressure_load.h
//  MAST
//
//  Created by Manav Bhatia on 8/1/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

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
//#ifdef LIBMESH_USE_COMPLEX_NUMBERS
            
            MeshBase& linear_sys_mesh = linearized_sys.get_mesh();
            _linear_mesh_serializer.reset(new MeshSerializer(linear_sys_mesh, true));
            
            MeshBase& nonlinear_sys_mesh = nonlinear_sys.get_mesh();
            _nonlinear_mesh_serializer.reset(new MeshSerializer(nonlinear_sys_mesh, true));
            
            
            _dim = linearized_sys.n_vars()-2;
            
            // copy the pointer for flight condition data
            _flt_cond = dynamic_cast<FrequencyDomainLinearizedFluidSystem&>
            (lin_sys).flight_condition;
//#endif
        }
        
        virtual ~SmallDisturbanceSurfacePressure()
        { }
        
        virtual void init(libMesh::NumericVector<libMesh::Real>& nonlinear_sol,
                          libMesh::NumericVector<libMesh::Real>& linearized_sol);
        
        // calculation in frequency domain
        template <typename ValType>
        void surface_pressure(const libMesh::Point& p,
                              libMesh::Real& cp,
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
        std::auto_ptr<MeshFunction> _function_nonlinear;
        
        /*!
         *   mesh function that interpolates the linearized solution
         */
        std::auto_ptr<MeshFunction> _function_linear;
        
        /*!
         *    numeric vector that stores the solution for nonlinear system
         */
        std::auto_ptr<libMesh::NumericVector<libMesh::Real> > _sol_nonlinear;
        
        /*!
         *    numeric vector that stores the solution for linearized system
         */
        std::auto_ptr<libMesh::NumericVector<libMesh::Real> > _sol_linear;
        
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
        std::auto_ptr<MeshSerializer> _nonlinear_mesh_serializer;
        std::auto_ptr<MeshSerializer> _linear_mesh_serializer;
        
    };
}



#endif
