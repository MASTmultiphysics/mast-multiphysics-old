//
//  beam_postbuckling_sizing.cpp
//  RealSolver
//
//  Created by Manav Bhatia on 6/16/14.
//  Copyright (c) 2014 Manav Bhatia. All rights reserved.
//

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
            
            Real h, rho;
            v = 0.;
            
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

                (*stiff_area)(p, 0., h);
                rhof(p, 0., rho);
                
                v += e->volume() * h * rho;
            }
        }
        
        virtual void partial(const MAST::FieldFunctionBase& f,
                             const libMesh::Point& p, Real t, Real& v) const {
            libmesh_error();
        }
        
        virtual void total(const MAST::FieldFunctionBase& f,
                           const libMesh::Point& p, Real t, Real& v) const {
            
            v = 0.;
            
            libMesh::MeshBase::const_element_iterator
            eit  = _mesh.active_local_elements_begin(),
            eend = _mesh.active_local_elements_end();
            
            Real h, rho;
            v = 0.;
            
            for ( ; eit != eend; eit++ ) {
                const libMesh::Elem* e = *eit;
                const MAST::ElementPropertyCardBase& prop =
                _assembly.get_property_card(*e);
                const MAST::MaterialPropertyCardBase& mat =
                prop.get_material();
                const MAST::FieldFunction<Real> &rhof =
                mat.get<MAST::FieldFunction<Real> >("rho");
                const MAST::Solid1DSectionElementPropertyCard& prop1d =
                dynamic_cast<const MAST::Solid1DSectionElementPropertyCard&>(prop);
                std::auto_ptr<MAST::FieldFunction<Real> >
                stiff_area (prop1d.section_property<MAST::FieldFunction<Real> >("A").release());
                
                stiff_area->total(f, p, 0., h);
                rhof(p, 0., rho);
                
                v += e->volume() * h * rho;
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
            
            delete _zero_function;
            
            delete _weight;
            
            for (unsigned int i=0; i<_elem_properties.size(); i++)
                delete _elem_properties[i];
            
            for (unsigned int i=0; i<_materials.size(); i++)
                delete _materials[i];
            
            delete _E;
            delete _nu;
            delete _alpha;
            delete _rho;
            delete _kappa;
            
            for (unsigned int i=0; i<_h_y.size(); i++) {
                delete _h_y[i];
                delete _h_z[i];
                delete _offset_h_y[i];
                delete _offset_h_z[i];
            }

            delete _prestress;
            delete _temperature;
            delete _ref_temperature;
            delete _temperature_bc;
            
            for (unsigned int i=0; i<_n_vars; i++) {
                delete _disp_function_sens[i];
            }
        }
        
        
        
        virtual void init_dvar(std::vector<Real>& x,
                               std::vector<Real>& xmin,
                               std::vector<Real>& xmax);
        
        
        
        virtual void evaluate_func(const std::vector<Real>& dvars,
                                   Real& obj,
                                   bool eval_obj_grad,
                                   std::vector<Real>& obj_grad,
                                   std::vector<Real>& fvals,
                                   std::vector<bool>& eval_grads,
                                   std::vector<Real>& grads);

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
        
        MAST::Weight* _weight;
        
        libMesh::EquationSystems* _eq_systems, *_fluid_eq_systems;
        
        libMesh::NonlinearImplicitSystem *_static_system;
        
        libMesh::CondensedEigenSystem *_eigen_system;
        
        FluidSystem *_fluid_system_nonlin;
        
        FrequencyDomainLinearizedFluidSystem *_fluid_system_freq;
        
        MAST::StructuralSystemAssembly *_static_structural_assembly,
        *_eigen_structural_assembly;
        
        libMesh::UnstructuredMesh* _mesh, *_fluid_mesh;
        
        ConstantFunction<Real>* _press;
        
        libMesh::ZeroFunction<Real>* _zero_function;
        
        MAST::ConstantFunction<Real> *_E, *_nu,  *_alpha, *_rho, *_kappa;
        
        std::vector<MAST::ConstantFunction<Real>*> _h_y, _h_z,
        _offset_h_z;
        
        std::vector<MAST::BeamOffset*> _offset_h_y;
        
        MAST::ConstantFunction<DenseRealMatrix > *_prestress;
        
        MAST::ConstantFunction<Real> *_temperature, *_ref_temperature;
        
        MAST::Temperature *_temperature_bc;
        
        std::auto_ptr<MAST::ConstantFunction<Real> > _pressure;

        std::auto_ptr<MAST::BoundaryCondition> _pressure_bc;
        
        libMesh::ParameterVector _parameters;
        
        std::vector<MAST::ConstantFunction<Real>*> _parameter_functions;
        
        std::vector<MAST::DisplacementDirichletBoundaryCondition*> _bc;
        
        std::vector<MAST::MaterialPropertyCardBase*> _materials;
        
        std::vector<MAST::ElementPropertyCardBase*> _elem_properties;
        
        std::auto_ptr<libMesh::MeshFunction> _disp_function;

        std::vector<libMesh::MeshFunction*> _disp_function_sens;
        
        std::auto_ptr<FlightCondition> _flight_cond;
        
        std::auto_ptr<FEMStructuralModel> _structural_model;
        
        std::auto_ptr<CFDAerodynamicModel> _aero_model;
        
        std::auto_ptr<CoupledFluidStructureSystem> _coupled_system;
        
        std::auto_ptr<MAST::UGFlutterSolver> _flutter_solver;
        
        std::auto_ptr<MAST::SurfaceMotionBase> _surface_motion;
    };
}



inline
void
MAST::SizingOptimization::init_dvar(std::vector<Real>& x,
                                    std::vector<Real>& xmin,
                                    std::vector<Real>& xmax) {
    // one DV for each element
    x.resize   (_n_vars);
    std::fill(x.begin(), x.end(), 0.08);
    for (unsigned int i=1; i<x.size(); i+=2)
        x[i] = .09;
    xmin.resize(_n_vars);
    std::fill(xmin.begin(), xmin.end(), 6.0e-3);
    xmax.resize(_n_vars);
    std::fill(xmax.begin(), xmax.end(), 0.2);
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
    std::vector<Real> dv0 = dvars;
    
    std::vector<bool> eval_grads_false = eval_grads;
    std::fill(eval_grads_false.begin(), eval_grads_false.end(), false);
    
    
    evaluate_func(dv0, obj, false, obj_grad, fvals, eval_grads_false, grads);
    
    Real delta = 1.0e-4, dobj = 0.;
    if (eval_obj_grad) {
        std::vector<Real> dfvals = fvals;
        std::fill(dfvals.begin(), dfvals.end(), 0.);
        dv0[0] += delta;
        
        evaluate_func(dv0, dobj, false, obj_grad, dfvals, eval_grads_false, grads);
        obj_grad[0] = (dobj-obj)/delta;
        
        for (unsigned int i=0; i<_n_ineq; i++)
            grads[i] = (dfvals[i]-fvals[i])/delta;
    }
    
}


inline
void
MAST::SizingOptimization::evaluate_func(const std::vector<Real>& dvars,
                                        Real& obj,
                                        bool eval_obj_grad,
                                        std::vector<Real>& obj_grad,
                                        std::vector<Real>& fvals,
                                        std::vector<bool>& eval_grads,
                                        std::vector<Real>& grads) {
    
    
    libmesh_assert_equal_to(dvars.size(), _n_vars);
    
    // set the parameter values equal to the DV value
    for (unsigned int i=0; i<_n_vars; i++)
        *_parameters[i] = dvars[i];
    
    libMesh::Point pt; // dummy point object
    
    // calculate weight
    (*_weight)(pt, 0., obj);
    
    // calculate sensitivity of weight if requested
    if (eval_obj_grad) {
        std::fill(obj_grad.begin(), obj_grad.end(), 0.);
        for (unsigned int i=0; i<_n_vars; i++)
            _weight->total(*_parameter_functions[i], pt, 0., obj_grad[i]);
    }
    
    // now solve the system
    libMesh::out << "New Eval" << std::endl;
    
    // first zero the solution and init the Euler variables to undisturbed solution
    _static_system->solution->zero();
    init_euler_variables(*_fluid_eq_systems, "FluidSystem");
    
    // static analysis is performed with von-Karman strain
    for (unsigned int p=0; p<_elem_properties.size(); p++)
        _elem_properties[p]->set_strain(MAST::VON_KARMAN_STRAIN);
    // increase the load over several load steps
    const unsigned int n_load_steps = 100;
    Real temp_val, ref_temp, p_val;
    (*_temperature)(pt, 0., temp_val);
    (*_ref_temperature)(pt, 0., ref_temp);
    (*_pressure)(pt, 0., p_val);
    bool continue_fsi_iterations = true;

    for (unsigned int i=0; i<n_load_steps; i++) {
        std::cout << "Solving load step: " << i << std::endl;

        // use this displacement to converge the fluid solver
        continue_fsi_iterations = true;
        while (continue_fsi_iterations) {
            _fluid_system_nonlin->time = 0.;
            _fluid_system_nonlin->solve();
            if (_fluid_system_nonlin->time_solver->du(libMesh::SystemNorm()) <= 1.0e-6)
                continue_fsi_iterations = false;
        }
        
        // increment the temperature to the next load step
        (*_temperature) = ref_temp + (temp_val-ref_temp)*i/(n_load_steps-1);
        
        // now solve for this load step
        _static_system->solve();
    }
    
    // eigen analysis is performed with von-Karman strain
    for (unsigned int p=0; p<_elem_properties.size(); p++)
        _elem_properties[p]->set_strain(MAST::VON_KARMAN_STRAIN);
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
            
            // get the mode
            std::ostringstream vec_nm;
            vec_nm << "mode_" << i;
            libMesh::NumericVector<Real>& vec = _eigen_system->get_vector(vec_nm.str());
            vec = *(_eigen_system->solution);
            vec.scale(sqrt(eigval.real()));
            vec.close();
                        
            // rescale so that the inner product with the mass matrix is identity
            
            
            //            {
            //                std::ostringstream file_name;
            //
            //                // We write the file in the ExodusII format.
            //                file_name << "out_"
            //                << std::setw(3)
            //                << std::setfill('0')
            //                << std::right
            //                << i
            //                << ".exo";
            //
            //                // now write the eigenvlaues
            //                _system->get_eigenpair(i);
            //
            //                // We write the file in the ExodusII format.
            //                libMesh::Nemesis_IO(*_mesh).write_equation_systems(file_name.str(),
            //                                                          *_eq_systems);
            //            }
        }
    else
        libmesh_error(); // should not get here.
    
    
    // now get the displacement constraint
    pt(0) = 3.0;
    DenseRealVector disp;
    (*_disp_function)(pt, 0., disp);
    // reference displacement value
    Real ref_disp = 0.01;
    // w < w0 => w/w0 < 1. => w/w0 - 1. < 0
    fvals[0] = disp(0)/ref_disp-1.;

    // flutter solution
    _flutter_solver->clear_solutions();
    _flutter_solver->scan_for_roots();
    if (!_libmesh_init.comm().rank())
        _flutter_solver->print_crossover_points();
    std::pair<bool, const MAST::FlutterRootBase*> root =
    _flutter_solver->find_critical_root();
    libmesh_assert(root.first);
    if (!_libmesh_init.comm().rank())
        _flutter_solver->print_sorted_roots();

    // flutter constraint
    Real ref_V = 1000.;
    // vf > v0 => vf/v0 > 1 => 1-vf/v0 < 0
    //
    fvals[1] = 1. - root.second->V/ref_V;;
    
    std::vector<Real> grad_vals;
    
    if (eval_grads[0]) {
        // grad_k = dfi/dxj  ,  where k = j*NFunc + i
        // static analysis is performed without von-Karman strain
        for (unsigned int p=0; p<_elem_properties.size(); p++)
            _elem_properties[p]->set_strain(MAST::VON_KARMAN_STRAIN);
        _static_system->sensitivity_solve(_parameters);
        // eigen analysis is performed with von-Karman strain
        for (unsigned int p=0; p<_elem_properties.size(); p++)
            _elem_properties[p]->set_strain(MAST::VON_KARMAN_STRAIN);
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
            disp.zero();
            (*_disp_function_sens[j])(pt, 0., disp);
            grads[(j+1)*_n_ineq-1] = disp(0);
        }
    }
}


inline
void
MAST::SizingOptimization::_init() {
    
    _mesh = new libMesh::SerialMesh(_libmesh_init.comm());
    
    const unsigned int
    dim     = _infile("dimension",0),
    nx_divs = _infile("nx_divs",0),
    ny_divs = _infile("ny_divs",0),
    nz_divs = _infile("nz_divs",0);
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
    _n_vars = 1;
    
    _n_eq = 0;
    _n_ineq = (1 +   // +1 for the displacement
               1);   // +1 for flutter constraint
    _max_iters = 10000;
    
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
    var_id["ux"] = _static_system->add_variable ( "sux", static_cast<libMesh::Order>(o), fefamily);
    var_id["uy"] = _static_system->add_variable ( "suy", static_cast<libMesh::Order>(o), fefamily);
    var_id["uz"] = _static_system->add_variable ( "suz", static_cast<libMesh::Order>(o), fefamily);
    var_id["tx"] = _static_system->add_variable ( "stx", static_cast<libMesh::Order>(o), fefamily);
    var_id["ty"] = _static_system->add_variable ( "sty", static_cast<libMesh::Order>(o), fefamily);
    var_id["tz"] = _static_system->add_variable ( "stz", static_cast<libMesh::Order>(o), fefamily);
    
    _eigen_system->add_variable ( "ux", static_cast<libMesh::Order>(o), fefamily);
    _eigen_system->add_variable ( "uy", static_cast<libMesh::Order>(o), fefamily);
    _eigen_system->add_variable ( "uz", static_cast<libMesh::Order>(o), fefamily);
    _eigen_system->add_variable ( "tx", static_cast<libMesh::Order>(o), fefamily);
    _eigen_system->add_variable ( "ty", static_cast<libMesh::Order>(o), fefamily);
    _eigen_system->add_variable ( "tz", static_cast<libMesh::Order>(o), fefamily);
    
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
    _temperature = new MAST::ConstantFunction<Real>("temp", _infile("panel_temperature", 303.15)); // K
    _ref_temperature = new MAST::ConstantFunction<Real>("ref_temp", _infile("panel_ref_temperature", 303.15)); // K
    _temperature_bc = new MAST::Temperature;
    _temperature_bc->set_function(*_temperature);
    _temperature_bc->set_reference_temperature_function(*_ref_temperature);
    _static_structural_assembly->add_volume_load(0, *_temperature_bc);
    _eigen_structural_assembly->add_volume_load(0, *_temperature_bc);
    
    // pressure boundary condition
    _pressure.reset(new MAST::ConstantFunction<Real>("pressure", -0.));
    _pressure_bc.reset(new MAST::BoundaryCondition(MAST::SURFACE_PRESSURE));
    _pressure_bc->set_function(*_pressure);
    _static_structural_assembly->add_volume_load(0, *_pressure_bc);
    _eigen_structural_assembly->add_volume_load(0, *_pressure_bc);
    
    
    // apply the boundary conditions
    _zero_function = new libMesh::ZeroFunction<Real>;
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
    _eigen_structural_assembly->set_static_solution_system(_static_system);
    
    
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
    
    _surface_motion.reset(new MAST::FlexibleSurfaceMotion(*_static_system));
    _fluid_system_nonlin->surface_motion = _surface_motion.get();
    
    _eq_systems->init ();
    _fluid_eq_systems->init();
    

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
    _materials.resize(1);
    _parameters.resize(_n_vars);
    _parameter_functions.resize(_n_vars);
    
    DenseRealMatrix prestress; prestress.resize(3,3);
    prestress(0,0) = -1.31345e6;
    
    _E = new MAST::ConstantFunction<Real>("E", _infile("youngs_modulus", 72.e9)),
    _nu = new MAST::ConstantFunction<Real>("nu", _infile("poisson_ratio", 0.33)),
    _rho = new MAST::ConstantFunction<Real>("rho", _infile("material_density", 2700.)),
    _kappa = new MAST::ConstantFunction<Real>("kappa", _infile("shear_corr_factor", 5./6.)),
    _alpha = new MAST::ConstantFunction<Real>("alpha", _infile("expansion_coefficient", 2.31e-5)),
    
    _prestress = new MAST::ConstantFunction<DenseRealMatrix >("prestress", prestress);
    _materials[0] = new MAST::IsotropicMaterialPropertyCard(0);
    
    MAST::MaterialPropertyCardBase& mat = *_materials[0];
    // add the properties to the cards
    mat.add(*_E);
    mat.add(*_nu);
    mat.add(*_alpha);
    mat.add(*_rho);
    mat.add(*_kappa);
    
    _elem_properties.resize(1);
    // stiffener modeled as beam
    _h_y.resize(_n_vars);
    _h_z.resize(_n_vars);
    _offset_h_y.resize(_n_vars);
    _offset_h_z.resize(_n_vars);
    
    _elem_properties[0] = new MAST::Solid2DSectionElementPropertyCard(1); // panel
    
    // create values for each stiffener
    for (unsigned int i=0; i<_n_vars; i++) {
        MAST::Solid1DSectionElementPropertyCard *p = new MAST::Solid1DSectionElementPropertyCard(2);
        _h_y[i] = new MAST::ConstantFunction<Real>("hy", _infile("thickness", 0.002));
        _h_z[i] = new MAST::ConstantFunction<Real>("hz", _infile("width", 0.002)),
        _offset_h_y[i] = new MAST::BeamOffset("hy_offset", _h_y[i]->clone().release()),
        _offset_h_z[i] = new MAST::ConstantFunction<Real>("hz_offset", 0.);
        p->add(*_h_y[i]); // thickness
        p->add(*_h_z[i]); // width
        p->add(*_offset_h_y[i]); // thickness offset
        p->add(*_offset_h_z[i]); // width offset
        p->y_vector()(1) = 1.; // x-vector along x, y along y
        
        p->set_material(mat);
        _elem_properties[i] = p;
        
        _parameters[i] = _h_y[i]->ptr(); // set thickness as a modifiable parameter
        
        _static_structural_assembly->add_parameter(*_h_y[i]);
        _eigen_structural_assembly->add_parameter(*_h_y[i]);
        _parameter_functions[i] = _h_y[i];
        
        _static_structural_assembly->set_property_for_subdomain(0, *_elem_properties[i]);
        _eigen_structural_assembly->set_property_for_subdomain(0, *_elem_properties[i]);
    }
    
    // create the function to calculate weight
    _weight = new MAST::Weight(*_mesh, *_static_structural_assembly);
    
    // create the mesh function to calculate the displacement
    std::vector<unsigned int> vars(1);
    vars[0] = _static_system->variable_number("suy");
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
