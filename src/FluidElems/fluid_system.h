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

#ifndef __MAST__assembleEuler__
#define __MAST__assembleEuler__

// libmesh includes
#include "libmesh/libmesh_config.h"
#include "libmesh/fem_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/mesh_function.h"
#include "libmesh/mesh_serializer.h"

// MAST includes
#include "FluidElems/fluid_elem_base.h"
#include "Numerics/function_base.h"




Real euler_solution_value(const libMesh::Point& p,
                          const libMesh::Parameters& parameters,
                          const std::string& sys_name,
                          const std::string& var_name);



void init_euler_variables(libMesh::EquationSystems& es,
                          const std::string& system_name);



class FluidSystem : public libMesh::FEMSystem, public FluidElemBase
{
public:
    // Constructor
    FluidSystem(libMesh::EquationSystems& es,
                const std::string& name_in,
                const unsigned int number_in)
    : libMesh::FEMSystem(es, name_in, number_in),
    FluidElemBase(*es.parameters.get<GetPot*>("input_file")),
    _rho_norm_old(1.),
    _rho_norm_curr(1.),
    dc_recalculate_tolerance(1.0e-8),
    if_use_stored_dc_coeff(false)
    {}
    
    void init_data();
    
    virtual void init_context(libMesh::DiffContext &context);
    
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

protected:
    
    /*!
     *   localized reference solution for calculation of dc_val
     */
    std::auto_ptr<libMesh::NumericVector<Real> > _dc_ref_sol;
    
    /*!
     *    Current and old norms of density in the flow-field
     */
    Real _rho_norm_old, _rho_norm_curr;
};





class FluidPostProcessSystem : public libMesh::System
{
public:
    // Constructor
    FluidPostProcessSystem(libMesh::EquationSystems& es,
                           const std::string& name_in,
                           const unsigned int number_in)
    : libMesh::System(es, name_in, number_in)
    {}
    
    virtual void init_data();
    
    virtual void postprocess();
    
    unsigned int u, v, w, T, s, p, cp, a, M;
    
    FlightCondition* flight_condition;
};



/*!
 *   This provides a function to calculate the surface pressure
 *   from a fluid FEM solution.
 */
class CFDPressure: public MAST::FieldFunction<Real> {
public:
    CFDPressure(const std::string& nm,
                FluidSystem& sys):
    MAST::FieldFunction<Real>(nm),
    _system(sys),
    _primitive(new PrimitiveSolution) {
        
        libMesh::MeshBase& mesh = sys.get_mesh();
        _mesh_serializer.reset(new libMesh::MeshSerializer(mesh, true));
    }
    
    
    CFDPressure(const CFDPressure& press):
    MAST::FieldFunction<Real>(press),
    _system(press._system),
    _primitive(new PrimitiveSolution) {
        
        libMesh::MeshBase& mesh = _system.get_mesh();
        _mesh_serializer.reset(new libMesh::MeshSerializer(mesh, true));
    }
    
    
    virtual ~CFDPressure() {
        
    }
    
    
    /*!
     *   @returns a clone of the function
     */
    virtual std::auto_ptr<MAST::FieldFunction<Real> > clone() const {
        return std::auto_ptr<MAST::FieldFunction<Real> >
        (new CFDPressure(*this));
    }

    
    void init(libMesh::NumericVector<Real>& sol) {
        // first initialize the solution to the given vector
        if (!_sol.get())
        {
            _sol.reset(libMesh::NumericVector<Real>::build(_system.comm()).release());
            _sol->init(sol.size(), true, libMesh::SERIAL);
        }
        
        // now localize the give solution to this objects's vector
        sol.localize(*_sol);
        
        // if the mesh function has not been created so far, initialize it
        if (!_function.get()) {
            
            const unsigned int n_vars = _system.n_vars();
            std::vector<unsigned int> vars(n_vars);
            vars[0] = _system.variable_number("rho");
            vars[1] = _system.variable_number("rhoux");
            if (n_vars > 3) // > 1D
                vars[2] = _system.variable_number("rhouy");
            if (n_vars > 4) // > 2D
                vars[3] = _system.variable_number("rhouz");
            vars[n_vars-1] = _system.variable_number("rhoe");
            _function.reset(new libMesh::MeshFunction(_system.get_equation_systems(),
                                                      *_sol,
                                                      _system.get_dof_map(),
                                                      vars));
            _function->init();
        }
    }
    
    
    /*!
     *   returns the value of this function
     */
    virtual void operator() (const libMesh::Point& p,
                             const Real t,
                             Real& v) const {
        
        // zero the data structures
        _primitive->zero();
        
        // get the solution value from the fluid system
        DenseRealVector fluid_sol;
        (*_function)(p, t, fluid_sol);
        
        _primitive->init(_system.dim,
                         fluid_sol,
                         _system.flight_condition->gas_property.cp,
                         _system.flight_condition->gas_property.cv,
                         _system.if_viscous());
        
        v = _primitive->p - _system.flight_condition->gas_property.pressure;
    }
    
    
    /*!
     *   returns the sensitivity of this function
     */
    virtual void partial (const MAST::FieldFunctionBase& f,
                          const libMesh::Point& p,
                          const Real t,
                          Real& v) const {
        libmesh_error(); // to be implemented
    }
    
    
    /*!
     *   returns the sensitivity of this function. This is the same as partial
     *   sensitivity for a constant function.
     */
    virtual void total (const MAST::FieldFunctionBase& f,
                        const libMesh::Point& p,
                        const Real t,
                        Real& v) const {
        libmesh_error(); // to be implemented
    }
    
    
protected:
    
    /*!
     *   libmesh system that provides the fluid FEM solution space
     */
    FluidSystem& _system;
    
    /*!
     *   mesh function that interpolates the solution
     */
    std::auto_ptr<libMesh::MeshFunction> _function;
    
    /*!
     *    numeric vector that stores the solution
     */
    std::auto_ptr<libMesh::NumericVector<Real> > _sol;
    
    
    /*!
     *   this serializes the mesh for use in interpolation
     */
    std::auto_ptr<libMesh::MeshSerializer> _mesh_serializer;
    
    /*!
     *  object to calculate the primite variable values given conservative 
     *  variables
     */
    std::auto_ptr<PrimitiveSolution> _primitive;
    
};



#endif /* defined(__MAST__assembleEuler__) */
