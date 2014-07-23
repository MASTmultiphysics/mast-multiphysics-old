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

#ifndef __MAST_beam_postbuckling_sizing_optimization_h__
#define __MAST_beam_postbuckling_sizing_optimization_h__



// MAST includes
#include "Optimization/optimization_interface.h"
#include "StructuralElems/structural_system_assembly.h"
#include "PropertyCards/material_property_card_base.h"
#include "PropertyCards/solid_1d_section_element_property_card.h"
#include "PropertyCards/solid_2d_section_element_property_card.h"
#include "PropertyCards/multilayer_1d_section_element_property_card.h"
#include "PropertyCards/multilayer_2d_section_element_property_card.h"
#include "BoundaryConditions/temperature.h"
#include "BoundaryConditions/displacement_boundary_condition.h"
#include "Mesh/mesh_initializer.h"
#include "PropertyCards/isotropic_material_property_card.h"
#include "FluidElems/frequency_domain_linearized_fluid_system.h"
#include "Flight/flight_condition.h"
#include "FluidElems/fluid_system.h"
#include "FluidElems/frequency_domain_linearized_fluid_system.h"
#include "Aeroelasticity/ug_flutter_solver.h"
#include "Aeroelasticity/coupled_fluid_structure_system.h"
#include "BoundaryConditions/flexible_surface_motion.h"


// libmesh includes
#include "libmesh/getpot.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/nemesis_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/eigen_solver.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/parameter_vector.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/const_function.h"
#include "libmesh/zero_function.h"
#include "libmesh/nonlinear_solver.h"
#include "libmesh/nemesis_io.h"
#include "libmesh/sensitivity_data.h"
#include "libmesh/condensed_eigen_system.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/mesh_function.h"
#include "libmesh/steady_solver.h"
#include "libmesh/newton_solver.h"
#include "libmesh/euler2_solver.h"
#include "libmesh/system_norm.h"



namespace MAST {
    
    class MultilinearInterpolation: public MAST::FieldFunction<Real> {
    public:
        MultilinearInterpolation(const std::string& nm,
                                 std::map<Real, MAST::FieldFunction<Real>*>& values):
        MAST::FieldFunction<Real>(nm),
        _values(values) {
            
            // make sure that the size of the provided values is finite
            libmesh_assert(values.size() > 0);
            
            std::map<Real, MAST::FieldFunction<Real>*>::iterator
            it = values.begin(), end = values.end();

            // tell the function that it is dependent on the provided functions
            for ( ; it != end; it++)
                _functions.insert(it->second->master());
        }
        
        MultilinearInterpolation(const MAST::MultilinearInterpolation& o):
        MAST::FieldFunction<Real>(o),
        _values(o._values) {
            std::map<Real, MAST::FieldFunction<Real>*>::iterator
            it = _values.begin(), end = _values.end();
            
            // tell the function that it is dependent on the provided functions
            for ( ; it != end; it++)
                _functions.insert(it->second->master());
        }
        
        virtual std::auto_ptr<MAST::FieldFunction<Real> >
        clone() const {
            return std::auto_ptr<MAST::FieldFunction<Real> >
            (new MAST::MultilinearInterpolation(*this));
        }
        
        virtual ~MultilinearInterpolation() {
            
        }
        
    protected:
        
        std::map<Real, MAST::FieldFunction<Real>*> _values;
        
    public:
        
        virtual void operator() (const libMesh::Point& p, Real t, Real& v) const {
            
            //
            // the following is used for calculation of the return value
            //   f(x) is defined for x for each x0 < x < x1
            //   if   x <= x0,      f(x) = f(x0)
            //   if   x0 < x < x1,  f(x) is interpolated
            //   if   x >= x1,      f(x) = f(x1)
            //
            
            std::map<Real, MAST::FieldFunction<Real>*>::const_iterator
            it1, it2;
            std::map<Real, MAST::FieldFunction<Real>*>::const_reverse_iterator
            rit = _values.rbegin();
            it1  = _values.begin();
            
            // check the lower bound
            if (p(0) <=  it1->first) {
                (*it1->second)(p, t, v);
            }
            // check the upper bound
            else if (p(0) >=  rit->first) {
                (*rit->second)(p, t, v);
            }
            else {
                // if it gets here, the ordinate is in between the provided range
                it2 = _values.lower_bound(p(0));
                // this cannot be the first element of the map
                libmesh_assert(it2 != _values.begin());
                // it2 provides the upper bound. The lower bound is provided by the
                // preceding iterator
                it1 = it2--;
                Real f0 = 0., f1 = 0.;
                (*it1->second)(p, t, f0);
                (*it2->second)(p, t, f1);
                // now interpolate
                v =  (f0 +
                      (p(0) - it1->first)/(it2->first - it1->first) *
                      (f1-f0));
            }
        }
        
        virtual void partial(const MAST::FieldFunctionBase& f,
                             const libMesh::Point& p, Real t, Real& v) const {
            libmesh_error();
        }
        
        virtual void total(const MAST::FieldFunctionBase& f,
                           const libMesh::Point& p, Real t, Real& v) const {
            
            //
            // the following is used for calculation of the return value
            //   f(x) is defined for x for each x0 < x < x1
            //   if   x <= x0,      f(x) = f(x0)
            //   if   x0 < x < x1,  f(x) is interpolated
            //   if   x >= x1,      f(x) = f(x1)
            //
            
            std::map<Real, MAST::FieldFunction<Real>*>::const_iterator
            it1, it2;
            std::map<Real, MAST::FieldFunction<Real>*>::const_reverse_iterator
            rit = _values.rbegin();
            it1  = _values.begin();
            
            // check the lower bound
            if (p(0) <=  it1->first) {
                (*it1->second)(p, t, v);
            }
            // check the upper bound
            else if (p(0) >=  rit->first) {
                (*rit->second)(p, t, v);
            }
            else {
                // if it gets here, the ordinate is in between the provided range
                it2 = _values.lower_bound(p(0));
                // this cannot be the first element of the map
                libmesh_assert(it2 != _values.begin());
                // it2 provides the upper bound. The lower bound is provided by the
                // preceding iterator
                it1 = it2--;
                Real f0 = 0., f1 = 0.;
                it1->second->total(f, p, t, f0);
                it2->second->total(f, p, t, f1);
                // now interpolate
                v =  (f0 +
                      (p(0) - it1->first)/(it2->first - it1->first) *
                      (f1-f0));
            }
        }
    };
    
    
    class BeamOffset: public MAST::FieldFunction<Real> {
    public:
        BeamOffset(const std::string& nm,
                   MAST::FieldFunction<Real> *thickness):
        MAST::FieldFunction<Real>(nm),
        _dim(thickness) {
            _functions.insert(thickness->master());
        }
        
        BeamOffset(const MAST::BeamOffset& o):
        MAST::FieldFunction<Real>(o),
        _dim(o._dim->clone().release()) {
            _functions.insert(_dim->master());
        }
        
        virtual std::auto_ptr<MAST::FieldFunction<Real> >
        clone() const {
            return std::auto_ptr<MAST::FieldFunction<Real> >
            (new MAST::BeamOffset(*this));
        }
        
        virtual ~BeamOffset() {
            delete _dim;
        }
        
    protected:
        
        MAST::FieldFunction<Real> *_dim;
        
    public:
        
        virtual void operator() (const libMesh::Point& p, Real t, Real& v) const {
            (*_dim)(p, t, v);
            v *= 0.5;
        }
        
        virtual void partial(const MAST::FieldFunctionBase& f,
                             const libMesh::Point& p, Real t, Real& v) const {
            libmesh_error();
        }
        
        virtual void total(const MAST::FieldFunctionBase& f,
                           const libMesh::Point& p, Real t, Real& v) const {
            _dim->total(f, p, t, v);
            v *= 0.5;
        }
        
    };
    
    
    class Weight: public MAST::FieldFunction<Real> {
    public:
        Weight(libMesh::UnstructuredMesh& m,
               MAST::StructuralSystemAssembly& assembly):
        MAST::FieldFunction<Real>("Weight"),
        _mesh(m),
        _assembly(assembly) {
            
        }
        
        Weight(const MAST::Weight& w):
        MAST::FieldFunction<Real>(w),
        _mesh(w._mesh),
        _assembly(w._assembly) {
            
        }
        
        virtual std::auto_ptr<MAST::FieldFunction<Real> >
        clone() const {
            return std::auto_ptr<MAST::FieldFunction<Real> >
            (new MAST::Weight(*this));
        }
        
        virtual ~Weight() { }
        
    protected:
        
        libMesh::UnstructuredMesh &_mesh;
        
        MAST::StructuralSystemAssembly& _assembly;
        
    public:
        
        virtual void operator() (const libMesh::Point& p, Real t, Real& v) const {
            libMesh::MeshBase::const_element_iterator
            eit  = _mesh.active_local_elements_begin(),
            eend = _mesh.active_local_elements_end();
            
            Real h, rho, x0, x1, dx;
            v = 0.;
            
            libMesh::Point elem_p;
            const unsigned int n_sec = 3; // number of quadrature divs
            
            for ( ; eit != eend; eit++ ) {
                const libMesh::Elem* e = *eit;
                const MAST::ElementPropertyCardBase& prop =
                _assembly.get_property_card(*e);
                const MAST::Solid1DSectionElementPropertyCard& prop1d =
                dynamic_cast<const MAST::Solid1DSectionElementPropertyCard&>(prop);
                std::auto_ptr<MAST::FieldFunction<Real> >
                stiff_area (prop1d.section_property<MAST::FieldFunction<Real> >("A").release());
                const MAST::MaterialPropertyCardBase& mat =
                prop.get_material();
                const MAST::FieldFunction<Real> &rhof =
                mat.get<MAST::FieldFunction<Real> >("rho");

                // for each element iterate over the length and calculate the
                // weight from the section area and section density
                // use three point trapezoidal rule to calculate the integral
                x0 = e->point(0)(0);
                x1 = e->point(1)(0);
                dx = (x1-x0)/n_sec;
                for (unsigned int i=0; i<n_sec; i++) {
                    elem_p(0) = x0 + dx*(i+0.5);
                    (*stiff_area)(elem_p, 0., h);
                    rhof(elem_p, 0., rho);
                    v += h * rho * dx;
                }
                
            }
        }
        
        virtual void partial(const MAST::FieldFunctionBase& f,
                             const libMesh::Point& p, Real t, Real& v) const {
            libmesh_error();
        }
        
        virtual void total(const MAST::FieldFunctionBase& f,
                           const libMesh::Point& p, Real t, Real& v) const {
            libMesh::MeshBase::const_element_iterator
            eit  = _mesh.active_local_elements_begin(),
            eend = _mesh.active_local_elements_end();
            
            Real h, rho, dh, drho, x0, x1, dx;
            v = 0.;
            
            libMesh::Point elem_p;
            const unsigned int n_sec = 3; // number of quadrature divs
            
            for ( ; eit != eend; eit++ ) {
                const libMesh::Elem* e = *eit;
                const MAST::ElementPropertyCardBase& prop =
                _assembly.get_property_card(*e);
                const MAST::Solid1DSectionElementPropertyCard& prop1d =
                dynamic_cast<const MAST::Solid1DSectionElementPropertyCard&>(prop);
                std::auto_ptr<MAST::FieldFunction<Real> >
                stiff_area (prop1d.section_property<MAST::FieldFunction<Real> >("A").release());
                const MAST::MaterialPropertyCardBase& mat =
                prop.get_material();
                const MAST::FieldFunction<Real> &rhof =
                mat.get<MAST::FieldFunction<Real> >("rho");
                
                // for each element iterate over the length and calculate the
                // weight from the section area and section density
                // use three point trapezoidal rule to calculate the integral
                x0 = e->point(0)(0);
                x1 = e->point(1)(0);
                dx = (x1-x0)/n_sec;
                for (unsigned int i=0; i<n_sec; i++) {
                    elem_p(0) = x0 + dx*(i+0.5);
                    (*stiff_area)(elem_p, 0., h);
                    stiff_area->total(f, elem_p, 0., dh);
                    rhof(elem_p, 0., rho);
                    rhof.total(f, elem_p, 0., drho);
                    v += (dh * rho + h * drho) * dx;
                }
                
            }
        }
    };
    
    
    class SizingOptimization:
    public MAST::FunctionEvaluation {
        
    public:
        
        SizingOptimization(libMesh::LibMeshInit& init,
                           GetPot& input,
                           std::ostream& output):
        MAST::FunctionEvaluation(output),
        _libmesh_init(init),
        _infile(input),
        _fluid_input(GetPot("system_input.in")),
        _disp_0(0.),
        _vf_0(0.),
        _n_eig(0),
        _weight(NULL),
        _eq_systems(NULL),
        _static_system(NULL),
        _eigen_system(NULL),
        _static_structural_assembly(NULL),
        _eigen_structural_assembly(NULL),
        _mesh(NULL),
        _press(NULL),
        _zero_function(NULL)
        {
            _init();
        }
        
        virtual ~SizingOptimization() {
            
            delete _static_structural_assembly;
            
            delete _eigen_structural_assembly;
            
            delete _eq_systems;
            
            delete _fluid_eq_systems;
            
            delete _mesh;
            
            delete _fluid_mesh;
            
            for (unsigned int i=0; i<_bc.size(); i++)
                delete _bc[i];
            
            delete _weight;
            
            // delete the density station values
            std::map<Real, MAST::FieldFunction<Real>*>::iterator it, end;
            it  = _rho_station_vals.begin();
            end = _rho_station_vals.end();
            for (; it != end; it++)
                delete it->second;
            
            // delete the h_y station values
            it  = _h_y_station_vals.begin();
            end = _h_y_station_vals.end();
            for (; it != end; it++)
                delete it->second;
            
            for (unsigned int i=0; i<_n_vars; i++) {
                delete _disp_function_sens[i];
            }
        }
        
        
        
        virtual void init_dvar(std::vector<Real>& x,
                               std::vector<Real>& xmin,
                               std::vector<Real>& xmax);
        
        
        
        virtual void evaluate(const std::vector<Real>& dvars,
                              Real& obj,
                              bool eval_obj_grad,
                              std::vector<Real>& obj_grad,
                              std::vector<Real>& fvals,
                              std::vector<bool>& eval_grads,
                              std::vector<Real>& grads);
        
        virtual void output(unsigned int iter,
                            const std::vector<Real>& x,
                            Real obj,
                            const std::vector<Real>& fval) const;
        
    protected:
        
        
        void _init();
        
        
        libMesh::LibMeshInit& _libmesh_init;
        
        GetPot& _infile, _fluid_input;
        
        unsigned int _n_eig;
        
        Real _disp_0, _vf_0;
        
        MAST::Weight* _weight;
        
        libMesh::EquationSystems* _eq_systems, *_fluid_eq_systems;
        
        libMesh::NonlinearImplicitSystem *_static_system;
        
        libMesh::CondensedEigenSystem *_eigen_system;
        
        FluidSystem *_fluid_system_nonlin;
        
        FrequencyDomainLinearizedFluidSystem *_fluid_system_freq;
        
        MAST::StructuralSystemAssembly *_static_structural_assembly,
        *_eigen_structural_assembly;
        
        libMesh::UnstructuredMesh* _mesh, *_fluid_mesh;
        
        std::auto_ptr<ConstantFunction<Real> > _press;
        
        std::auto_ptr<libMesh::ZeroFunction<Real> > _zero_function;
        
        std::auto_ptr<MAST::ConstantFunction<Real> > _E, _nu, _alpha, _kappa;
        
        /*!
         *    map of density functions at discrete stations
         */
        std::map<Real, MAST::FieldFunction<Real>*> _rho_station_vals;
        
        /*!
         *    multilinear density function
         */
        std::auto_ptr<MAST::MultilinearInterpolation> _rho;
        
        /*!
         *    map of thickness functions at discrete stations
         */
        std::map<Real, MAST::FieldFunction<Real>*> _h_y_station_vals;
        
        /*!
         *    multilinear thickness function
         */
        std::auto_ptr<MAST::MultilinearInterpolation> _h_y;

        std::auto_ptr<MAST::ConstantFunction<Real> >  _h_z, _offset_h_z;
        
        std::auto_ptr<MAST::BeamOffset> _offset_h_y;
        
        std::auto_ptr<MAST::ConstantFunction<DenseRealMatrix> > _prestress;
        
        std::auto_ptr<MAST::ConstantFunction<Real> >_temperature, _ref_temperature;
        
        std::auto_ptr<MAST::Temperature> _temperature_bc;
        
        libMesh::ParameterVector _parameters;
        
        std::vector<MAST::ConstantFunction<Real>*> _parameter_functions;
        
        std::vector<MAST::DisplacementDirichletBoundaryCondition*> _bc;
        
        std::auto_ptr<MAST::MaterialPropertyCardBase> _materials;
        
        std::auto_ptr<MAST::ElementPropertyCardBase> _elem_properties;
        
        std::auto_ptr<libMesh::MeshFunction> _disp_function;

        std::vector<libMesh::MeshFunction*> _disp_function_sens;
        
        std::auto_ptr<FlightCondition> _flight_cond;
        
        std::auto_ptr<FEMStructuralModel> _structural_model;
        
        std::auto_ptr<CFDAerodynamicModel> _aero_model;
        
        std::auto_ptr<CoupledFluidStructureSystem> _coupled_system;
        
        std::auto_ptr<MAST::UGFlutterSolver> _flutter_solver;
        
        std::auto_ptr<MAST::FlexibleSurfaceMotion> _surface_motion;
        
        std::auto_ptr<MAST::ConstantFunction<Real> > _pressure;
        
        std::auto_ptr<CFDPressure> _cfd_pressure;

        std::auto_ptr<MAST::BoundaryCondition> _pressure_bc;
        
        std::vector<Real> _dv_scaling, _dv_low, _dv_init;
        
    };
}



inline
void
MAST::SizingOptimization::init_dvar(std::vector<Real>& x,
                                    std::vector<Real>& xmin,
                                    std::vector<Real>& xmax) {
    // one DV for each element
    x       = _dv_init;
    xmin    = _dv_low;
    xmax.resize(_n_vars);
    std::fill(xmax.begin(), xmax.end(), 1.);
}



inline
void
MAST::SizingOptimization::evaluate(const std::vector<Real>& dvars,
                                   Real& obj,
                                   bool eval_obj_grad,
                                   std::vector<Real>& obj_grad,
                                   std::vector<Real>& fvals,
                                   std::vector<bool>& eval_grads,
                                   std::vector<Real>& grads) {
    
    
    libmesh_assert_equal_to(dvars.size(), _n_vars);
    
    // set the parameter values equal to the DV value
    for (unsigned int i=0; i<_n_vars; i++)
        *_parameters[i] = dvars[i]*_dv_scaling[i];
    
    // DO NOT zero out the gradient vector, since GCMMA needs it for the
    // subproblem solution
    //std::fill(obj_grad.begin(), obj_grad.end(), 0.);
    //std::fill(grads.begin(), grads.end(), 0.);
    
    libMesh::Point pt; // dummy point object
    
    libMesh::out << "New Eval" << std::endl;

    // the optimization problem is defined as
    // max Vf subject to constraint on weight
    Real wt = 0., vf = 0., disp = 0.;
    // calculate weight
    (*_weight)(pt, 0., wt);
    
    // first zero the solution and init the Euler variables to undisturbed solution
    _static_system->solution->zero();
    _fluid_system_nonlin->time = 0.;
    init_euler_variables(*_fluid_eq_systems, "FluidSystem");
    
    // increase the load over several load steps
    const unsigned int n_load_steps = 000;
    Real temp_val, ref_temp;
    (*_temperature)(pt, 0., temp_val);
    (*_ref_temperature)(pt, 0., ref_temp);

    bool continue_fsi_iterations = true;

    for (unsigned int i=0; i<n_load_steps; i++) {
        
        libMesh::out
        << "+++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl
        << "Solving load step: " << i << std::endl
        << "+++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;

        // use this displacement to converge the fluid solver
        continue_fsi_iterations = true;
        
        while (continue_fsi_iterations) {
            // initialize the fluid solution
            _fluid_system_nonlin->time = 0.;

            _fluid_system_nonlin->print_residual_norms = true;
            _fluid_system_nonlin->print_residuals = false;
            _fluid_system_nonlin->print_jacobian_norms = false;
            _fluid_system_nonlin->print_jacobians = false;
            
            libMesh::NewtonSolver &solver = dynamic_cast<libMesh::NewtonSolver&>
            (*(_fluid_system_nonlin->time_solver->diff_solver().get()));
            solver.quiet = false;
            solver.verbose = !solver.quiet;
            solver.brent_line_search = false;
            solver.max_nonlinear_iterations = 1;
            solver.relative_step_tolerance = 1.0e-6;
            solver.relative_residual_tolerance = 1.0e-6;
            solver.absolute_residual_tolerance = 1.0e-6;
            solver.continue_after_backtrack_failure = true;
            solver.continue_after_max_iterations = true;
            solver.require_residual_reduction = false;
            
            // And the linear solver options
            solver.max_linear_iterations = 1000;
            solver.initial_linear_tolerance = 1.0e-8;


            // init with zero frequency, since we are working with steady solver
            _surface_motion->init(0., 0., *(_static_system->solution));

            // check if the DC operator needs to be reevaluated
            _fluid_system_nonlin->evaluate_recalculate_dc_flag();
            
            // now solve the fluid system
            _fluid_system_nonlin->solve();
            Real du = _fluid_system_nonlin->time_solver->du(libMesh::SystemNorm());
            libMesh::out << "du = " << du << std::endl;
            if (du <= 1.0e-0 ||  // this is a very relaxed set of requirements input
                (_fluid_system_nonlin->time > _fluid_system_nonlin->deltat * 25)) // max 100 iterations
                continue_fsi_iterations = false;
            _fluid_system_nonlin->time_solver->advance_timestep();
        }
        
        // increment the temperature to the next load step
        (*_temperature) = ref_temp + (temp_val-ref_temp)*i/(n_load_steps-1);
        
        // initialize the pressure boundary condition
        _cfd_pressure->init(*(_fluid_system_nonlin->solution));
        
        // now solve for this load step
        _static_system->solve();
        
        // write both the fluid and structural systems
        {
            std::set<std::string> names; names.insert("FluidSystem");
            libMesh::ExodusII_IO(*_fluid_mesh).write_equation_systems("fluid.exo", *_fluid_eq_systems, &names);
        }
        {
            std::set<std::string> names; names.insert("StaticStructuralSystem");
            libMesh::ExodusII_IO(*_mesh).write_equation_systems("str.exo", *_eq_systems, &names);
        }
    }
    
    // eigen analysis is performed with von-Karman strain
    _eigen_system->solve();
    
    // the number of converged eigenpairs could be different from the number asked
    unsigned int n_required = std::min(_n_eig, _eigen_system->get_n_converged());
    
    std::pair<Real, Real> val;
    Complex eigval;
    std::fill(fvals.begin(), fvals.end(), 0.);
    
    _structural_model->eigen_vals.resize(_n_eig);

    if (_eq_systems->parameters.get<bool>("if_exchange_AB_matrices"))
        // the total number of eigenvalues is _n_eig, but only n_required are usable
        for (unsigned int i=0; i<n_required; i++) {
            val = _eigen_system->get_eigenpair(i);
            eigval = std::complex<Real>(val.first, val.second);
            eigval = 1./eigval;
            //fvals[i] = 1.-eigval.real(); // g <= 0.
            
            // copy modal data into the structural model for flutter analysis
            _structural_model->eigen_vals(i) = eigval.real();
            
            libMesh::out << eigval.real() << std::endl;
            
            // get the mode
            std::ostringstream vec_nm;
            vec_nm << "mode_" << i;
            libMesh::NumericVector<Real>& vec = _eigen_system->get_vector(vec_nm.str());
            _eigen_system->solution->scale(sqrt(eigval.real()));
            vec = *(_eigen_system->solution);
            vec.close();
                        
            // rescale so that the inner product with the mass matrix is identity
            
            
            {
                std::ostringstream file_name;
                
                // We write the file in the ExodusII format.
                file_name << "out_"
                << std::setw(3)
                << std::setfill('0')
                << std::right
                << i
                << ".exo";
                
                // We write the file in the ExodusII format.
                libMesh::Nemesis_IO(*_mesh).write_equation_systems(file_name.str(),
                                                                   *_eq_systems);
            }
        }
    else
        libmesh_error(); // should not get here.
    
    
    // now get the displacement constraint
    pt(0) = 3.;
    DenseRealVector disp_vec;
    (*_disp_function)(pt, 0., disp_vec);
    // reference displacement value
    // w < w0 => w/w0 < 1. => w/w0 - 1. < 0
    disp = disp_vec(0);
    
    // flutter solution
    _flutter_solver->clear_solutions();
    _flutter_solver->scan_for_roots();
    if (!_libmesh_init.comm().rank())
        _flutter_solver->print_crossover_points();
    std::pair<bool, const MAST::FlutterRootBase*> root =
    _flutter_solver->find_critical_root(1.e-3, 5);
    if (root.first) {
        if (!_libmesh_init.comm().rank())
            _flutter_solver->print_sorted_roots();
        vf = root.second->V;
    }
    else
        // no root was identified. So, set a feasible constraint value
        vf = 1.0e6;
    
    std::vector<Real> grad_vals;
    
    // set the function and objective values
    // flutter objective
    // vf > v0 => vf/v0 > 1 => 1-vf/v0 < 0
    //
    obj = wt;
    fvals[0] = -wt;// disp/disp_0-1.;
    fvals[1] =   1.-vf/_vf_0;
    Real w_sens = 0.;
    
    // evaluate sensitivity if any of the function sensitivity is required
    bool if_sens = (false || eval_obj_grad);
    libMesh::out << eval_obj_grad << "  ";
    for (unsigned int i=0; i<eval_grads.size(); i++) {
        libMesh::out << eval_grads[i] << "  ";
        if_sens = (if_sens || eval_grads[i]);
    }
    libMesh::out << std::endl;
    
    if (if_sens) {
        // grad_k = dfi/dxj  ,  where k = j*NFunc + i
        /*// static analysis is performed without von-Karman strain
        _static_system->sensitivity_solve(_parameters);
        // eigen analysis is performed with von-Karman strain
        _eigen_system->sensitivity_solve(_parameters, grad_vals);
        // now correct this for the fact that the matrices were exchanged
        if (_eq_systems->parameters.get<bool>("if_exchange_AB_matrices"))
            for (unsigned int j=0; j<_n_vars; j++)
                for (unsigned int i=0; i<n_required; i++) {
                    val = _eigen_system->get_eigenpair(i);
                    grads[j*_n_ineq+i]  = grad_vals[j*_eigen_system->get_n_converged()+i];
                    grads[j*_n_ineq+i] /= -1. * -pow(val.first, 2); // sens = - d eig / dp
                }
        else
            libmesh_error(); // should not get here
        
        // get the displacement gradient
        pt(0) = 3.0;
        for (unsigned int j=0; j<_n_vars; j++) {
            disp_vec.zero();
            (*_disp_function_sens[j])(pt, 0., disp_vec);
            grads[(j+1)*_n_ineq-1] = disp_vec(0);
        }
         */
        // ask flutter solver for the sensitivity
        
        // set gradient of weight
        for (unsigned int i=0; i<_n_vars; i++) {
            _weight->total(*_parameter_functions[i], pt, 0., w_sens);
            obj_grad[i] = w_sens*_dv_scaling[i];
            grads[i*_n_ineq] = -w_sens*_dv_scaling[i];
        }
        
        // now get the flutter sensitivity
        // if the flutter was not found, then the sensitivity is zero.
        if (root.first) {
            // create a copy of the root and use it to calculate the
            // sensitivity of the flutter velocity
            MAST::FlutterRootBase root_sens(*root.second);
            for (unsigned int i=0; i<_n_vars; i++) {
                _flutter_solver->calculate_sensitivity(root_sens,
                                                       _parameters,
                                                       i);
                grads[i*_n_ineq+1] = -root_sens.V_sens/_vf_0*_dv_scaling[i];
            }
        }
    }
}


inline
void
MAST::SizingOptimization::_init() {
    
    _mesh = new libMesh::SerialMesh(_libmesh_init.comm());
    
    const unsigned int
    n_stations    = _infile("n_stations", 0),
    dim           = _infile("dimension",0),
    nx_divs       = _infile("nx_divs",0),
    ny_divs       = _infile("ny_divs",0),
    nz_divs       = _infile("nz_divs",0);
    libMesh::ElemType elem_type =
    libMesh::Utility::string_to_enum<libMesh::ElemType>(_infile("elem_type", "QUAD4"));
    
    std::vector<Real> x_div_loc(nx_divs+1), x_relative_dx(nx_divs+1),
    y_div_loc(ny_divs+1), y_relative_dx(ny_divs+1),
    z_div_loc(nz_divs+1), z_relative_dx(nz_divs+1);
    std::vector<unsigned int> x_divs(nx_divs), y_divs(ny_divs), z_divs(nz_divs);
    std::auto_ptr<MeshInitializer::CoordinateDivisions>
    x_coord_divs (new MeshInitializer::CoordinateDivisions),
    y_coord_divs (new MeshInitializer::CoordinateDivisions),
    z_coord_divs (new MeshInitializer::CoordinateDivisions);
    std::vector<MeshInitializer::CoordinateDivisions*> divs(dim);
    
    // now read in the values: x-coord
    if (nx_divs > 0)
    {
        for (unsigned int i_div=0; i_div<nx_divs+1; i_div++)
        {
            x_div_loc[i_div]     = _infile("x_div_loc", 0., i_div);
            x_relative_dx[i_div] = _infile( "x_rel_dx", 0., i_div);
            if (i_div < nx_divs) //  this is only till nx_divs
                x_divs[i_div] = _infile( "x_div_nelem", 0, i_div);
        }
        divs[0] = x_coord_divs.get();
        x_coord_divs->init(nx_divs, x_div_loc, x_relative_dx, x_divs);
    }
    
    // design data
    _n_eig = 6;
    _n_vars = n_stations*2;
    
    _n_eq = 0;
    _n_ineq = (1 +   // +1 for the displacement
               1);   // +1 for flutter speed
    _max_iters = 10000;
    
    // initialize the dv vector data
    Real
    th_l  = _infile("thickness_lower", 0.001),
    th_u  = _infile("thickness_upper", 0.2),
    th    = _infile("thickness", 0.01),
    rho_l = _infile("material_density_lower", 2500.),
    rho_u = _infile("material_density_upper", 10000.),
    rho   = _infile("material_density", 2700.);
    
    _disp_0 = _infile( "displacement_0", 0.);
    _vf_0   = _infile("flutter_speed_0", 0.);
    
    _dv_init.resize(_n_vars);
    _dv_scaling.resize(_n_vars);
    _dv_low.resize(_n_vars);
    
    for (unsigned int i=0; i<_n_vars/2; i++) {
        // first half is the thickness values
        _dv_init[i]    =   th/th_u;
        _dv_low[i]     = th_l/th_u;
        _dv_scaling[i] =      th_u;
        
        // next half is the density values
        _dv_init[i+_n_vars/2]    =  rho/rho_u;
        _dv_low[i+_n_vars/2]     =  rho_l/rho_u;
        _dv_scaling[i+_n_vars/2] =  rho_u;
    }

    
    
    // now initialize the mesh
    MeshInitializer().init(divs, *_mesh, elem_type);
    
    // Print information about the mesh to the screen.
    _mesh->print_info();
    
    // Create an equation systems object.
    _eq_systems = new libMesh::EquationSystems(*_mesh);
    _eq_systems->parameters.set<GetPot*>("input_file") = &_infile;
    
    // Declare the system
    _static_system = &_eq_systems->add_system<libMesh::NonlinearImplicitSystem> ("StaticStructuralSystem");
    _eigen_system = &_eq_systems->add_system<libMesh::CondensedEigenSystem> ("EigenStructuralSystem");
    
    unsigned int o = _infile("fe_order", 1);
    std::string fe_family = _infile("fe_family", std::string("LAGRANGE"));
    libMesh::FEFamily fefamily = libMesh::Utility::string_to_enum<libMesh::FEFamily>(fe_family);
    
    std::map<std::string, unsigned int> var_id;
    libMesh::Order order = static_cast<libMesh::Order>(o);
    var_id["ux"] = _static_system->add_variable ( "ux", order, fefamily);
    var_id["uy"] = _static_system->add_variable ( "uy", order, fefamily);
    var_id["uz"] = _static_system->add_variable ( "uz", order, fefamily);
    var_id["tx"] = _static_system->add_variable ( "tx", order, fefamily);
    var_id["ty"] = _static_system->add_variable ( "ty", order, fefamily);
    var_id["tz"] = _static_system->add_variable ( "tz", order, fefamily);
    
    _eigen_system->add_variable ( "ux", order, fefamily);
    _eigen_system->add_variable ( "uy", order, fefamily);
    _eigen_system->add_variable ( "uz", order, fefamily);
    _eigen_system->add_variable ( "tx", order, fefamily);
    _eigen_system->add_variable ( "ty", order, fefamily);
    _eigen_system->add_variable ( "tz", order, fefamily);
    
    _static_structural_assembly = new MAST::StructuralSystemAssembly(*_static_system,
                                                                     MAST::STATIC,
                                                                     _infile);
    
    _eigen_structural_assembly = new MAST::StructuralSystemAssembly(*_eigen_system,
                                                                    MAST::MODAL,
                                                                    _infile);
    
    _static_system->attach_assemble_object(*_static_structural_assembly);
    _static_system->attach_sensitivity_assemble_object(*_static_structural_assembly);
    
    _eigen_system->attach_assemble_object(*_eigen_structural_assembly);
    _eigen_system->attach_eigenproblem_sensitivity_assemble_object(*_eigen_structural_assembly);
    
    
    // temperature load
    _temperature.reset(new MAST::ConstantFunction<Real>
                       ("temp", _infile("panel_temperature", 353.15))); // K
    _ref_temperature.reset(new MAST::ConstantFunction<Real>
                           ("ref_temp", _infile("panel_ref_temperature", 303.15))); // K
    _temperature_bc.reset(new MAST::Temperature);
    _temperature_bc->set_function(*_temperature);
    _temperature_bc->set_reference_temperature_function(*_ref_temperature);
    //_static_structural_assembly->add_volume_load(0, *_temperature_bc);
    //_eigen_structural_assembly->add_volume_load(0, *_temperature_bc);
    
    // pressure boundary condition
    //_pressure.reset(new MAST::ConstantFunction<Real>("pressure", -0.));
    _pressure_bc.reset(new MAST::BoundaryCondition(MAST::SURFACE_PRESSURE));
    //_static_structural_assembly->add_volume_load(0, *_pressure_bc);
    //_eigen_structural_assembly->add_volume_load(0, *_pressure_bc);
    
    
    // apply the boundary conditions
    _zero_function.reset(new libMesh::ZeroFunction<Real>);
    // Pass the Dirichlet dof IDs to the libMesh::CondensedEigenSystem
    std::set<libMesh::boundary_id_type> dirichlet_boundary;
    // read and initialize the boundary conditions
    std::map<libMesh::boundary_id_type, std::vector<unsigned int> > boundary_constraint_map;
    unsigned int n_bc, b_id;
    // first read the boundaries for ux constraint
    
    for (std::map<std::string, unsigned int>::iterator it = var_id.begin();
         it != var_id.end(); it++) {
        
        std::string nm = "n_" + it->first + "_bc"; // name for # bcs
        n_bc = _infile(nm, 0);
        
        nm = it->first + "_bc";  // name for bc id vector
        
        for (unsigned int i=0; i<n_bc; i++) {
            b_id = _infile(nm, 0, i); // ith component of the bc id vector
            
            if (!boundary_constraint_map.count(b_id)) // add vector if it does not exist
                boundary_constraint_map[b_id] = std::vector<unsigned int>(0);
            
            boundary_constraint_map[b_id].push_back(var_id[it->first]);
        }
    }
    
    // now iterate over each boundary and create the boudnary condition object
    unsigned int counter=0;
    _bc.resize(boundary_constraint_map.size());
    for (std::map<libMesh::boundary_id_type, std::vector<unsigned int> >::iterator
         it = boundary_constraint_map.begin();
         it != boundary_constraint_map.end(); it++) {
        _bc[counter] = new MAST::DisplacementDirichletBoundaryCondition;
        _bc[counter]->init(it->first, it->second);
        _static_structural_assembly->add_side_load(it->first, *_bc[counter]);
        _eigen_structural_assembly->add_side_load(it->first, *_bc[counter]);
        counter++;
    }
    
    _eigen_system->set_eigenproblem_type(libMesh::GHEP);
    _eigen_system->eigen_solver->set_position_of_spectrum(libMesh::LARGEST_MAGNITUDE);
    //_eigen_structural_assembly->set_static_solution_system(_static_system);
    
    
    // initialize the fluid data structures
    
    // first the flight condition
    _flight_cond.reset(new FlightCondition);
    for (unsigned int i=0; i<3; i++)
    {
        _flight_cond->body_roll_axis(i)     = _fluid_input(    "body_roll_axis", 0., i);
        _flight_cond->body_pitch_axis(i)    = _fluid_input(   "body_pitch_axis", 0., i);
        _flight_cond->body_yaw_axis(i)      = _fluid_input(     "body_yaw_axis", 0., i);
        _flight_cond->body_euler_angles(i)  = _fluid_input( "body_euler_angles", 0., i);
        _flight_cond->body_angular_rates(i) = _fluid_input("body_angular_rates", 0., i);
    }
    _flight_cond->ref_chord       = _fluid_input("ref_c",   1.);
    _flight_cond->altitude        = _fluid_input( "alt",    0.);
    _flight_cond->mach            = _fluid_input("mach",    .5);
    _flight_cond->gas_property.cp = _fluid_input(  "cp", 1003.);
    _flight_cond->gas_property.cv = _fluid_input(  "cv",  716.);
    _flight_cond->gas_property.T  = _fluid_input("temp",  300.);
    _flight_cond->gas_property.rho= _fluid_input( "rho",  1.05);
    
    _flight_cond->init();

    // next the fluid mesh
    _fluid_mesh = new libMesh::SerialMesh(_libmesh_init.comm());

    const unsigned int fluid_dim     = _fluid_input("dimension",0),
    fluid_nx_divs = _fluid_input("nx_divs",0),
    fluid_ny_divs = _fluid_input("ny_divs",0),
    fluid_nz_divs = _fluid_input("nz_divs",0);
    divs.resize(fluid_dim);
    libMesh::ElemType fluid_elem_type =
    libMesh::Utility::string_to_enum<libMesh::ElemType>(_fluid_input("elem_type", "QUAD4"));
    
    x_div_loc.resize(fluid_nx_divs+1), x_relative_dx.resize(fluid_nx_divs+1),
    y_div_loc.resize(fluid_ny_divs+1), y_relative_dx.resize(fluid_ny_divs+1),
    z_div_loc.resize(fluid_nz_divs+1), z_relative_dx.resize(fluid_nz_divs+1);
    x_divs.resize(fluid_nx_divs), y_divs.resize(fluid_ny_divs), z_divs.resize(fluid_nz_divs);
    std::auto_ptr<MeshInitializer::CoordinateDivisions>
    fluid_x_coord_divs (new MeshInitializer::CoordinateDivisions),
    fluid_y_coord_divs (new MeshInitializer::CoordinateDivisions),
    fluid_z_coord_divs (new MeshInitializer::CoordinateDivisions);
    std::vector<MeshInitializer::CoordinateDivisions*> fluid_divs(fluid_dim);
    
    // now read in the values: x-coord
    if (fluid_nx_divs > 0)
    {
        for (unsigned int i_div=0; i_div<fluid_nx_divs+1; i_div++)
        {
            x_div_loc[i_div]     = _fluid_input("x_div_loc", 0., i_div);
            x_relative_dx[i_div] = _fluid_input( "x_rel_dx", 0., i_div);
            if (i_div < fluid_nx_divs) //  this is only till nx_divs
                x_divs[i_div] = _fluid_input( "x_div_nelem", 0, i_div);
        }
        divs[0] = x_coord_divs.get();
        x_coord_divs->init(fluid_nx_divs, x_div_loc, x_relative_dx, x_divs);
    }
    
    
    if (fluid_ny_divs > 0)
    {
        for (unsigned int i_div=0; i_div<fluid_ny_divs+1; i_div++)
        {
            y_div_loc[i_div]     = _fluid_input("y_div_loc", 0., i_div);
            y_relative_dx[i_div] = _fluid_input( "y_rel_dx", 0., i_div);
            if (i_div < fluid_ny_divs) //  this is only till ny_divs
                y_divs[i_div] = _fluid_input( "y_div_nelem", 0, i_div);
        }
        divs[1] = y_coord_divs.get();
        y_coord_divs->init(fluid_ny_divs, y_div_loc, y_relative_dx, y_divs);
    }


    const bool if_cos_bump = _fluid_input("if_cos_bump", false);
    const unsigned int n_max_bumps_x = _fluid_input("n_max_bumps_x", 1),
    n_max_bumps_y = _fluid_input("n_max_bumps_x", 1),
    panel_bc_id = _fluid_input("panel_bc_id", 10),
    symmetry_bc_id = _fluid_input("symmetry_bc_id", 11);
    const Real t_by_c =  _fluid_input("t_by_c", 0.0);
    
    PanelMesh2D().init(t_by_c, if_cos_bump, n_max_bumps_x,
                       panel_bc_id, symmetry_bc_id,
                       divs, *_fluid_mesh, fluid_elem_type);

    // Print information about the mesh to the screen.
    _fluid_mesh->print_info();
    
    // Create an equation systems object.
    _fluid_eq_systems = new libMesh::EquationSystems(*_fluid_mesh);
    _fluid_eq_systems->parameters.set<GetPot*>("input_file") = &_fluid_input;
    
    // Declare the system
    _fluid_system_nonlin =
    &(_fluid_eq_systems->add_system<FluidSystem>("FluidSystem"));
    _fluid_system_freq =
    &(_fluid_eq_systems->add_system<FrequencyDomainLinearizedFluidSystem>
    ("FrequencyDomainLinearizedFluidSystem"));
    
    _fluid_system_nonlin->flight_condition = _flight_cond.get();
    _fluid_system_freq->flight_condition = _flight_cond.get();
    
    _fluid_system_nonlin->attach_init_function(init_euler_variables);
    
    _fluid_system_nonlin->time_solver =
    libMesh::AutoPtr<libMesh::TimeSolver>(new libMesh::Euler2Solver(*_fluid_system_nonlin));
    _fluid_system_freq->time_solver =
    libMesh::AutoPtr<libMesh::TimeSolver>(new libMesh::SteadySolver(*_fluid_system_freq));
    _fluid_system_freq->time_solver->quiet = false;

    // now initilaize the nonlinear solution
    _fluid_system_freq->extra_quadrature_order =
    _fluid_input("extra_quadrature_order", 0);
    _fluid_eq_systems->parameters.set<bool>("if_reduced_freq") =
    _fluid_input("if_reduced_freq", false);
    _fluid_system_nonlin->deltat = _fluid_input("deltat", 1.0e-2);
    
    _eq_systems->init ();
    _fluid_eq_systems->init();
    
    _surface_motion.reset(new MAST::FlexibleSurfaceMotion(*_static_system));
    _fluid_system_nonlin->surface_motion = _surface_motion.get();

    libMesh::NewtonSolver &solver = dynamic_cast<libMesh::NewtonSolver&>
    (*(_fluid_system_freq->time_solver->diff_solver()));
    solver.quiet = _fluid_input("solver_quiet", true);
    solver.verbose = !solver.quiet;
    solver.brent_line_search = false;
    solver.max_nonlinear_iterations =
    _fluid_input("max_nonlinear_iterations", 15);
    solver.relative_step_tolerance =
    _fluid_input("relative_step_tolerance", 1.e-3);
    solver.relative_residual_tolerance =
    _fluid_input("relative_residual_tolerance", 0.0);
    solver.absolute_residual_tolerance =
    _fluid_input("absolute_residual_tolerance", 0.0);
    solver.continue_after_backtrack_failure =
    _fluid_input("continue_after_backtrack_failure", false);
    solver.continue_after_max_iterations =
    _fluid_input("continue_after_max_iterations", false);
    solver.require_residual_reduction =
    _fluid_input("require_residual_reduction", true);
    
    // And the linear solver options
    solver.max_linear_iterations =
    _fluid_input("max_linear_iterations", 50000);
    solver.initial_linear_tolerance =
    _fluid_input("initial_linear_tolerance", 1.e-3);

    
    
    
    // the frequency domain
    _fluid_system_freq->localize_fluid_solution();

    
    // Print information about the system to the screen.
    _eq_systems->print_info();
    _fluid_eq_systems->print_info();
    
    
    // now initialize the flutter solver data structures
    _structural_model.reset(new FEMStructuralModel(*_eigen_structural_assembly));
    _aero_model.reset(new CFDAerodynamicModel(*_fluid_system_nonlin,
                                              *_fluid_system_freq));
    _coupled_system.reset(new CoupledFluidStructureSystem
                          (*_aero_model, *_structural_model));

    
    // create the solvers
    _flutter_solver.reset(new MAST::UGFlutterSolver);
    std::string nm = "flutter_output.txt";
    if (!_libmesh_init.comm().rank())
        _flutter_solver->set_output_file(nm);
    _flutter_solver->aero_structural_model   = _coupled_system.get();
    _flutter_solver->flight_condition        = _flight_cond.get();
    _flutter_solver->ref_val_range.first     = _fluid_input("ug_lower_k", 0.0);
    _flutter_solver->ref_val_range.second    = _fluid_input("ug_upper_k", 0.35);
    _flutter_solver->n_ref_val_divs          = _fluid_input("ug_k_divs", 10);

    
    // Pass the Dirichlet dof IDs to the libMesh::CondensedEigenSystem
    std::set<unsigned int> dirichlet_dof_ids;
    _eq_systems->parameters.set<bool>("if_exchange_AB_matrices") = true;
    _eq_systems->parameters.set<unsigned int>("eigenpairs")    = _n_eig;
    _eq_systems->parameters.set<unsigned int>("basis vectors") = _n_eig*3;
    _eq_systems->parameters.set<unsigned int>("nonlinear solver maximum iterations") =
    _infile("max_nonlinear_iterations", 5);
    _eigen_structural_assembly->get_dirichlet_dofs(dirichlet_dof_ids);
    _eigen_system->initialize_condensed_dofs(dirichlet_dof_ids);
    
    
    
    // element and material properties
    _parameters.resize(_n_vars);
    _parameter_functions.resize(_n_vars);
    
    DenseRealMatrix prestress; prestress.resize(3,3);
    //prestress(0,0) = -1.31345e6;
    
    _E.reset(new MAST::ConstantFunction<Real>
             ("E", _infile("youngs_modulus", 72.e9))),
    _nu.reset(new MAST::ConstantFunction<Real>
              ("nu", _infile("poisson_ratio", 0.33))),
    _kappa.reset(new MAST::ConstantFunction<Real>
                 ("kappa", _infile("shear_corr_factor", 5./6.))),
    _alpha.reset(new MAST::ConstantFunction<Real>
                 ("alpha", _infile("expansion_coefficient", 2.31e-5))),
    _prestress.reset(new MAST::ConstantFunction<DenseRealMatrix >
                     ("prestress", prestress));
    _materials.reset(new MAST::IsotropicMaterialPropertyCard(0));
    
    const Real
    x0 = _infile("x_div_loc", 0., 0),  // panel LE
    x1 = _infile("x_div_loc", 0., 1);  // panel TE
    Real dx = (x1-x0)/(n_stations-1);
    
    // create the density variables
    for (unsigned int i=0; i<n_stations; i++) {
        std::ostringstream oss;
        oss << "rho_" << i;
        
        MAST::ConstantFunction<Real>* rho
        = new MAST::ConstantFunction<Real>(oss.str(), _infile("rho", 2700.));
        
        // add this to the density map
        _rho_station_vals.insert(std::pair<Real, MAST::FieldFunction<Real>*>
                                 (x0+i*dx, rho));
        
        // add the function to the parameter map
        _parameters[i+n_stations]          = rho->ptr();
        _parameter_functions[i+n_stations] = rho;
        
        // tell the assembly system about the sensitvity parameter
        _static_structural_assembly->add_parameter(*rho);
        _eigen_structural_assembly->add_parameter(*rho);
    }
    
    // now create the density and give it to the property card
    _rho.reset(new MAST::MultilinearInterpolation("rho", _rho_station_vals));

    MAST::MaterialPropertyCardBase& mat = *_materials;
    // add the properties to the cards
    mat.add(*_E);
    mat.add(*_nu);
    mat.add(*_alpha);
    mat.add(*_rho);
    mat.add(*_kappa);
    
    // create the thickness variables
    for (unsigned int i=0; i<n_stations; i++) {
        std::ostringstream oss;
        oss << "h_y_" << i;
        
        MAST::ConstantFunction<Real>* h_y
        = new MAST::ConstantFunction<Real>(oss.str(), _infile("thickness", 0.002));
        
        // add this to the density map
        _h_y_station_vals.insert(std::pair<Real, MAST::FieldFunction<Real>*>
                                 (x0+i*dx, h_y));
        
        // add the function to the parameter set
        _parameters[i]          = h_y->ptr();
        _parameter_functions[i] = h_y;
        
        // tell the assembly system about the sensitvity parameter
        _static_structural_assembly->add_parameter(*h_y);
        _eigen_structural_assembly->add_parameter(*h_y);
    }
    
    // now create the density and give it to the property card
    _h_y.reset(new MAST::MultilinearInterpolation("hy", _h_y_station_vals));
    
    // create values for each stiffener
    _h_z.reset(new MAST::ConstantFunction<Real>("hz", _infile("width", 0.002))),
    _offset_h_y.reset(new MAST::BeamOffset("hy_offset", _h_y->clone().release())),
    _offset_h_z.reset(new MAST::ConstantFunction<Real>("hz_offset", 0.));

    MAST::Solid1DSectionElementPropertyCard *p = new MAST::Solid1DSectionElementPropertyCard(2);
    p->add(*_h_y); // thickness
    p->add(*_h_z); // width
    p->add(*_offset_h_y); // thickness offset
    p->add(*_offset_h_z); // width offset
    p->y_vector()(1) = 1.; // x-vector along x, y along y
    
    p->set_material(mat);
    p->set_strain(MAST::VON_KARMAN_STRAIN);
    _elem_properties.reset(p);
    
    _static_structural_assembly->set_property_for_subdomain(0, *_elem_properties);
    _eigen_structural_assembly->set_property_for_subdomain(0, *_elem_properties);
    
    // create the function to calculate weight
    _weight = new MAST::Weight(*_mesh, *_static_structural_assembly);
    
    // create the mesh function to calculate the displacement
    std::vector<unsigned int> vars(1);
    vars[0] = _static_system->variable_number("uy");
    _disp_function.reset(new libMesh::MeshFunction(_static_system->get_equation_systems(),
                                                   *_static_system->solution,
                                                   _static_system->get_dof_map(),
                                                   vars));
    _disp_function->init();
    
    _disp_function_sens.resize(_n_vars);
    for (unsigned int i=0; i<_n_vars; i++) {
        _static_system->add_sensitivity_solution(i);
        _disp_function_sens[i] = new libMesh::MeshFunction(_static_system->get_equation_systems(),
                                                           _static_system->get_sensitivity_solution(i),
                                                           _static_system->get_dof_map(),
                                                           vars);
        _disp_function_sens[i]->init();
    }
    
    // now add the vectors for structural modes, and init the fem structural model
    for (unsigned int i=0; i<_n_eig; i++) {
        std::ostringstream vec_nm;
        vec_nm << "mode_" << i ;
        _eigen_system->add_vector(vec_nm.str());
    }
    _structural_model->eigen_vals.resize(_n_eig);
    _structural_model->init();
    
    // initialize the pressure boundary condition
    _cfd_pressure.reset(new CFDPressure("cfd_pressure",
                                        *_fluid_system_nonlin));
    
    _pressure_bc->set_function(*_cfd_pressure);
}


inline void
MAST::SizingOptimization::output(unsigned int iter,
                                 const std::vector<Real>& x,
                                 Real obj,
                                 const std::vector<Real>& fval) const {
    
    libmesh_assert_equal_to(x.size(), _n_vars);
    
    // set the parameter values equal to the DV value
    for (unsigned int i=0; i<_n_vars; i++)
        *_parameters[i] = x[i];
    
    // the number of converged eigenpairs could be different from the number asked
    unsigned int n_required = std::min(_n_eig, _eigen_system->get_n_converged());
    
    std::pair<Real, Real> val;
    Complex eigval;
    
    for (unsigned int i=0; i<n_required; i++) {
        std::ostringstream file_name;
        
        // We write the file in the ExodusII format.
        file_name << "out_iter"
        << std::setw(3)
        << std::setfill('0')
        << std::right
        << iter << "_mode_"
        << std::setw(3)
        << std::setfill('0')
        << std::right
        << i
        << ".exo";
        
        // now write the eigenvlaues
        _eigen_system->get_eigenpair(i);
        
        // We write the file in the ExodusII format.
        libMesh::Nemesis_IO(*_mesh).write_equation_systems(file_name.str(),
                                                  *_eq_systems);
    }
    
    MAST::FunctionEvaluation::output(iter, x, obj, fval);
}


#endif // __MAST_beam_postbuckling_sizing_optimization_h__
