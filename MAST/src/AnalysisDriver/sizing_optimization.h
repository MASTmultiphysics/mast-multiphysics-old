//
//  sizing_optimization.h
//  MAST
//
//  Created by Manav Bhatia on 12/30/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_sizing_optimization_h__
#define __MAST_sizing_optimization_h__



// MAST includes
#include "Optimization/optimization_interface.h"
#include "StructuralElems/structural_system_assembly.h"
#include "PropertyCards/material_property_card_base.h"
#include "PropertyCards/solid_1d_section_element_property_card.h"
#include "PropertyCards/solid_2d_section_element_property_card.h"
#include "BoundaryConditions/temperature.h"

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
        
        virtual void operator() (const Point& p, Real t, Real& v) const {
            (*_dim)(p, t, v);
            v *= 0.5;
        }
        
        virtual void partial(const MAST::FieldFunctionBase& f,
                             const Point& p, Real t, Real& v) const {
            libmesh_error();
        }
        
        virtual void total(const MAST::FieldFunctionBase& f,
                           const Point& p, Real t, Real& v) const {
            _dim->total(f, p, t, v);
            v *= 0.5;
        }

    };
    
    
    class Weight: public MAST::FieldFunction<Real> {
    public:
        Weight(UnstructuredMesh& m,
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
        
        UnstructuredMesh &_mesh;
        
        MAST::StructuralSystemAssembly& _assembly;
        
    public:
        
        virtual void operator() (const Point& p, Real t, Real& v) const {
            MeshBase::const_element_iterator
            eit  = _mesh.active_local_elements_begin(),
            eend = _mesh.active_local_elements_end();
            
            Real h, rho;
            v = 0.;
            
            for ( ; eit != eend; eit++ ) {
                const Elem* e = *eit;
                const MAST::ElementPropertyCardBase& prop =
                _assembly.get_property_card(*e);
                const MAST::MaterialPropertyCardBase& mat =
                prop.get_material();
                const MAST::FieldFunction<Real> &rhof =
                mat.get<MAST::FieldFunction<Real> >("rho");
                
                if (e->dim() == 2) { // panel
                    const MAST::FieldFunction<Real> &hf =
                    prop.get<MAST::FieldFunction<Real> >("h");
                    
                    hf(p, 0., h);
                    rhof(p, 0., rho);
                    
                    v += e->volume() * h * rho;
                }
                else if (e->dim() == 1) { // stiffener
                    const MAST::Solid1DSectionElementPropertyCard& prop1d =
                    dynamic_cast<const MAST::Solid1DSectionElementPropertyCard&>(prop);
                    std::auto_ptr<MAST::FieldFunction<Real> >
                    stiff_area (prop1d.section_property<MAST::FieldFunction<Real> >("A").release());
                    
                    (*stiff_area)(p, 0., h);
                    rhof(p, 0., rho);
            
                    v += e->volume() * h * rho;
                    
                }
                
            }
        }
    
        virtual void partial(const MAST::FieldFunctionBase& f,
                             const Point& p, Real t, Real& v) const {
            libmesh_error();
        }
        
        virtual void total(const MAST::FieldFunctionBase& f,
                           const Point& p, Real t, Real& v) const {
            
            v = 0.;

            MeshBase::const_element_iterator
            eit  = _mesh.active_local_elements_begin(),
            eend = _mesh.active_local_elements_end();
            
            Real h, rho;
            v = 0.;
            
            for ( ; eit != eend; eit++ ) {
                const Elem* e = *eit;
                const MAST::ElementPropertyCardBase& prop =
                _assembly.get_property_card(*e);
                const MAST::MaterialPropertyCardBase& mat =
                prop.get_material();
                const MAST::FieldFunction<Real> &rhof =
                mat.get<MAST::FieldFunction<Real> >("rho");
                
                if (e->dim() == 2) { // panel
                    const MAST::FieldFunction<Real> &hf =
                    prop.get<MAST::FieldFunction<Real> >("h");
                    
                    hf.total(f, p, 0., h);
                    rhof(p, 0., rho);
                    
                    v += e->volume() * h * rho;
                }
                else if (e->dim() == 1) { // stiffener
                    const MAST::Solid1DSectionElementPropertyCard& prop1d =
                    dynamic_cast<const MAST::Solid1DSectionElementPropertyCard&>(prop);
                    std::auto_ptr<MAST::FieldFunction<Real> >
                    stiff_area (prop1d.section_property<MAST::FieldFunction<Real> >("A").release());
                    
                    stiff_area->total(f, p, 0., h);
                    rhof(p, 0., rho);
                    
                    v += e->volume() * h * rho;
                }
            }
        }
    };
    
    
    class SizingOptimization:
    public MAST::FunctionEvaluation {
        
    public:
        
        SizingOptimization(LibMeshInit& init,
                           GetPot& input,
                           std::ostream& output):
        MAST::FunctionEvaluation(output),
        _libmesh_init(init),
        _infile(input),
        _n_stiff(0),
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
            
            delete _mesh;
            
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
            delete _h;
            for (unsigned int i=0; i<_n_stiff; i++) {
                delete _h_stiff[i];
                delete _h_z[i];
                delete _offset_h_y[i];
                delete _offset_h_z[i];
            }
            
            delete _prestress;
            delete _temperature;
            delete _ref_temperature;
            delete _temperature_bc;
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
        
        
        LibMeshInit& _libmesh_init;
        
        GetPot& _infile;
        
        bool _beam_stiff;
        
        unsigned int _n_stiff, _n_eig;
        
        MAST::Weight* _weight;
        
        EquationSystems* _eq_systems;
        
        NonlinearImplicitSystem *_static_system;

        CondensedEigenSystem *_eigen_system;
        
        MAST::StructuralSystemAssembly* _static_structural_assembly,
        *_eigen_structural_assembly;
        
        UnstructuredMesh* _mesh;
        
        ConstantFunction<Real>* _press;
        
        ZeroFunction<Real>* _zero_function;
        
        MAST::ConstantFunction<Real> *_E, *_nu,  *_alpha, *_rho, *_kappa, *_h;
        
        std::vector<MAST::ConstantFunction<Real>*> _h_stiff, _h_z,
        _offset_h_y;

        std::vector<MAST::BeamOffset*> _offset_h_z;

        MAST::ConstantFunction<DenseMatrix<Real> > *_prestress;

        MAST::ConstantFunction<Real> *_temperature, *_ref_temperature;

        MAST::Temperature *_temperature_bc;
        
        ParameterVector _parameters;
        
        std::vector<MAST::ConstantFunction<Real>*> _parameter_functions;
        
        std::vector<MAST::DisplacementDirichletBoundaryCondition*> _bc;

        std::vector<MAST::MaterialPropertyCardBase*> _materials;
        
        std::vector<MAST::ElementPropertyCardBase*> _elem_properties;
    };
}



inline
void
MAST::SizingOptimization::init_dvar(std::vector<Real>& x,
                                      std::vector<Real>& xmin,
                                      std::vector<Real>& xmax) {
    // one DV for each element
    x.resize   (_n_vars);
    std::fill(x.begin(), x.end(), 0.008);
    for (unsigned int i=1; i<x.size(); i+=2)
        x[i] = .09;
    x[0] = 0.02;
    xmin.resize(_n_vars);
    std::fill(xmin.begin(), xmin.end(), 2.0e-3);
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
#ifndef LIBMESH_USE_COMPLEX_NUMBERS
    
    libmesh_assert_equal_to(dvars.size(), _n_vars);
    
    // set the parameter values equal to the DV value
    for (unsigned int i=0; i<_n_vars; i++)
        *_parameters[i] = dvars[i];
    
    Point p; // dummy point object

    // calculate weight
    (*_weight)(p, 0., obj);

    // calculate sensitivity of weight if requested
    if (eval_obj_grad) {
        std::fill(obj_grad.begin(), obj_grad.end(), 0.);
        for (unsigned int i=0; i<_n_vars; i++)
            _weight->total(*_parameter_functions[i], p, 0., obj_grad[i]);
    }

    // now solve the system
    libMesh::out << "New Eval" << std::endl;
    _static_system->solve();
    _eigen_system->solve();

    // the number of converged eigenpairs could be different from the number asked
    unsigned int n_required = std::min(_n_eig, _eigen_system->get_n_converged());
    
    std::pair<Real, Real> val;
    Complex eigval;
    std::fill(fvals.begin(), fvals.end(), 0.);
    if (_eq_systems->parameters.get<bool>("if_exchange_AB_matrices"))
        // the total number of constraints is _n_eig, but only n_required are usable
        for (unsigned int i=0; i<n_required; i++) {
            val = _eigen_system->get_eigenpair(i);
            eigval = std::complex<Real>(val.first, val.second);
            eigval = 1./eigval;
            fvals[i] = 1.-eigval.real(); // g <= 0.
            
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
//                Nemesis_IO(*_mesh).write_equation_systems(file_name.str(),
//                                                          *_eq_systems);
//            }
        }
    else
        libmesh_error(); // should not get here.
    
    std::vector<Real> grad_vals;
    
    if (eval_grads[0]) {
        // grad_k = dfi/dxj  ,  where k = j*NFunc + i
        _static_system->sensitivity_solve(_parameters);
        _eigen_system->sensitivity_solve(_parameters, grad_vals);
        // now correct this for the fact that the matrices were exchanged
        if (_eq_systems->parameters.get<bool>("if_exchange_AB_matrices"))
            for (unsigned int j=0; j<_n_vars; j++)
                for (unsigned int i=0; i<n_required; i++) {
                    val = _eigen_system->get_eigenpair(i);
                    grads[j*_n_eig+i]  = grad_vals[j*_eigen_system->get_n_converged()+i];
                    grads[j*_n_eig+i] /= -1. * -pow(val.first, 2); // sens = - d eig / dp
                }
        else
            libmesh_error(); // should not get here
    }
    
#endif // LIBMESH_USE_COMPLEX_NUMBERS
}


inline
void
MAST::SizingOptimization::_init() {
    
#ifndef LIBMESH_USE_COMPLEX_NUMBERS
    _mesh = new SerialMesh(_libmesh_init.comm());
    
    const unsigned int
    dim     = _infile("dimension",0),
    nx_divs = _infile("nx_divs",0),
    ny_divs = _infile("ny_divs",0),
    nz_divs = _infile("nz_divs",0);
    ElemType elem_type =
    Utility::string_to_enum<ElemType>(_infile("elem_type", "QUAD4"));
    
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
    
    // now read in the values: y-coord
    if ((dim > 1) && (ny_divs > 0))
    {
        for (unsigned int i_div=0; i_div<ny_divs+1; i_div++)
        {
            y_div_loc[i_div]     = _infile("y_div_loc", 0., i_div);
            y_relative_dx[i_div] = _infile( "y_rel_dx", 0., i_div);
            if (i_div < ny_divs) //  this is only till ny_divs
                y_divs[i_div] = _infile( "y_div_nelem", 0, i_div);
        }
        divs[1] = y_coord_divs.get();
        y_coord_divs->init(ny_divs, y_div_loc, y_relative_dx, y_divs);
    }
    
    // now read in the values: z-coord
    if ((dim == 3) && (nz_divs > 0))
    {
        for (unsigned int i_div=0; i_div<nz_divs+1; i_div++)
        {
            z_div_loc[i_div]     = _infile("z_div_loc", 0., i_div);
            z_relative_dx[i_div] = _infile( "z_rel_dx", 0., i_div);
            if (i_div < nz_divs) //  this is only till nz_divs
                z_divs[i_div] = _infile( "z_div_nelem", 0, i_div);
        }
        divs[2] = z_coord_divs.get();
        z_coord_divs->init(nz_divs, z_div_loc, z_relative_dx, z_divs);
    }
    
    // design data
    _n_eig = 10;
    _n_stiff = divs[0]->n_divs() + divs[1]->n_divs() - 2;
    _n_vars = _n_stiff*2+1;
    _n_eq = 0;
    _n_ineq = _n_eig;
    _max_iters = 10000;

    // now initialize the mesh
    _beam_stiff = _infile("beam_stiffeners", false);
    MAST::StiffenedPanel().init(divs, *_mesh, elem_type, _beam_stiff);

    // Print information about the mesh to the screen.
    _mesh->print_info();
    
    // Create an equation systems object.
    _eq_systems = new EquationSystems(*_mesh);
    _eq_systems->parameters.set<GetPot*>("input_file") = &_infile;
    
    // Declare the system
    _static_system = &_eq_systems->add_system<NonlinearImplicitSystem> ("StaticStructuralSystem");
    _eigen_system = &_eq_systems->add_system<CondensedEigenSystem> ("EigenStructuralSystem");
    
    unsigned int o = _infile("fe_order", 1);
    std::string fe_family = _infile("fe_family", std::string("LAGRANGE"));
    FEFamily fefamily = Utility::string_to_enum<FEFamily>(fe_family);
    
    std::map<std::string, unsigned int> var_id;
    var_id["ux"] = _static_system->add_variable ( "ux", static_cast<Order>(o), fefamily);
    var_id["uy"] = _static_system->add_variable ( "uy", static_cast<Order>(o), fefamily);
    var_id["uz"] = _static_system->add_variable ( "uz", static_cast<Order>(o), fefamily);
    var_id["tx"] = _static_system->add_variable ( "tx", static_cast<Order>(o), fefamily);
    var_id["ty"] = _static_system->add_variable ( "ty", static_cast<Order>(o), fefamily);
    var_id["tz"] = _static_system->add_variable ( "tz", static_cast<Order>(o), fefamily);
    
    _eigen_system->add_variable ( "ux", static_cast<Order>(o), fefamily);
    _eigen_system->add_variable ( "uy", static_cast<Order>(o), fefamily);
    _eigen_system->add_variable ( "uz", static_cast<Order>(o), fefamily);
    _eigen_system->add_variable ( "tx", static_cast<Order>(o), fefamily);
    _eigen_system->add_variable ( "ty", static_cast<Order>(o), fefamily);
    _eigen_system->add_variable ( "tz", static_cast<Order>(o), fefamily);
    
    _static_structural_assembly = new MAST::StructuralSystemAssembly(*_static_system,
                                                                     MAST::STATIC,
                                                                     _infile);

    _eigen_structural_assembly = new MAST::StructuralSystemAssembly(*_eigen_system,
                                                                    MAST::BUCKLING,
                                                                    _infile);

    _static_system->attach_assemble_object(*_static_structural_assembly);
    _static_system->attach_sensitivity_assemble_object(*_static_structural_assembly);

    _eigen_system->attach_assemble_object(*_eigen_structural_assembly);
    _eigen_system->attach_eigenproblem_sensitivity_assemble_object(*_eigen_structural_assembly);
    
    
    // temperature load
    _temperature = new MAST::ConstantFunction<Real>("temp", 505.95); // K
    _ref_temperature = new MAST::ConstantFunction<Real>("ref_temp", 232.15); // K
    _temperature_bc = new MAST::Temperature;
    _temperature_bc->set_function(*_temperature);
    _temperature_bc->set_reference_temperature_function(*_ref_temperature);
    _static_structural_assembly->add_volume_load(0, *_temperature_bc);
    _eigen_structural_assembly->add_volume_load(0, *_temperature_bc);
    
    
    // apply the boundary conditions
    _zero_function = new ZeroFunction<Real>;
    // Pass the Dirichlet dof IDs to the CondensedEigenSystem
    std::set<boundary_id_type> dirichlet_boundary;
    // read and initialize the boundary conditions
    std::map<boundary_id_type, std::vector<unsigned int> > boundary_constraint_map;
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
    for (std::map<boundary_id_type, std::vector<unsigned int> >::iterator
         it = boundary_constraint_map.begin();
         it != boundary_constraint_map.end(); it++) {
        _bc[counter] = new MAST::DisplacementDirichletBoundaryCondition;
        _bc[counter]->init(it->first, it->second);
        _static_structural_assembly->add_side_load(it->first, *_bc[counter]);
        _eigen_structural_assembly->add_side_load(it->first, *_bc[counter]);
        counter++;
    }

    _eigen_system->set_eigenproblem_type(GHEP);
    _eigen_system->eigen_solver->set_position_of_spectrum(LARGEST_MAGNITUDE);
    _eigen_structural_assembly->set_static_solution_system(_static_system);
    
    _eq_systems->init ();
    
    
    // Print information about the system to the screen.
    _eq_systems->print_info();

    // Pass the Dirichlet dof IDs to the CondensedEigenSystem
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
    _elem_properties.resize(_n_stiff+1);
    _parameters.resize(_n_vars);
    _parameter_functions.resize(_n_vars);
    
    DenseMatrix<Real> prestress; prestress.resize(3,3);
    prestress(0,0) = -1.31345e6;

    _E = new MAST::ConstantFunction<Real>("E", _infile("youngs_modulus", 72.e9)),
    _nu = new MAST::ConstantFunction<Real>("nu", _infile("poisson_ratio", 0.33)),
    _rho = new MAST::ConstantFunction<Real>("rho", _infile("material_density", 2700.)),
    _kappa = new MAST::ConstantFunction<Real>("kappa", _infile("shear_corr_factor", 5./6.)),
    _h = new MAST::ConstantFunction<Real>("h", _infile("thickness", 0.002));
    _alpha = new MAST::ConstantFunction<Real>("alpha", _infile("expansion_coefficient", 2.31e-5)),

    _prestress = new MAST::ConstantFunction<DenseMatrix<Real> >("prestress", prestress);
    
    _materials[0] = new MAST::IsotropicMaterialPropertyCard(0);
    _elem_properties[0] = new MAST::Solid2DSectionElementPropertyCard(1); // panel
    
    MAST::MaterialPropertyCardBase& mat = *_materials[0];
    // add the properties to the cards
    mat.add(*_E);
    mat.add(*_nu);
    mat.add(*_alpha);
    mat.add(*_rho);
    mat.add(*_kappa);

    // panel thickness as a design variable
    _parameters[0] = _h->ptr();
    _static_structural_assembly->add_parameter(*_h);
    _eigen_structural_assembly->add_parameter(*_h);
    _parameter_functions[0] = _h;

    _h_stiff.resize(_n_stiff);
    _h_z.resize(_n_stiff);
    _offset_h_y.resize(_n_stiff);
    _offset_h_z.resize(_n_stiff);

    if (_beam_stiff) {
        // create values for each stiffener
        for (unsigned int i=0; i<_n_stiff; i++) {
            MAST::Solid1DSectionElementPropertyCard *p = new MAST::Solid1DSectionElementPropertyCard(2);
            _h_stiff[i] = new MAST::ConstantFunction<Real>("hy", _infile("thickness", 0.002));
            _h_z[i] = new MAST::ConstantFunction<Real>("hz", z_div_loc[1]),
            _offset_h_y[i] = new MAST::ConstantFunction<Real>("hy_offset", 0.),
            _offset_h_z[i] = new MAST::BeamOffset("hz_offset", _h_z[i]->clone().release());
            p->add(*_h_stiff[i]); // width
            p->add(*_h_z[i]); // height
            p->add(*_offset_h_y[i]); // width offset
            p->add(*_offset_h_z[i]); // height offset
            if (i < ny_divs-1) // since stiffened_panel mesh first creates the x_stiffener
                p->y_vector()(1) = 1.; // x-vector along x, y along y
            else
                p->y_vector()(0) = -1.; // x-vector along y, y along -x
                
            p->set_material(mat);
            p->set_diagonal_mass_matrix(false);
            p->set_strain(MAST::VON_KARMAN_STRAIN);

            _elem_properties[i+1] = p;
            
            _parameters[i*2+1] = _h_z[i]->ptr(); // set thickness as a modifiable parameter
            _parameters[i*2+2] = _h_stiff[i]->ptr(); // set thickness as a modifiable parameter

            _static_structural_assembly->add_parameter(*_h_z[i]);
            _static_structural_assembly->add_parameter(*_h_stiff[i]);
            _eigen_structural_assembly->add_parameter(*_h_z[i]);
            _eigen_structural_assembly->add_parameter(*_h_stiff[i]);
            _parameter_functions[i*2+1] = _h_z[i];
            _parameter_functions[i*2+2] = _h_stiff[i];
            
            // temperature load for each stiffener domain
            _static_structural_assembly->add_volume_load(i+1, *_temperature_bc);
            _eigen_structural_assembly->add_volume_load(i+1, *_temperature_bc);
        }
    }
    else {
        libmesh_error(); // currently not implemented for non-beam stiffeners
//        _h_stiff = new MAST::ConstantFunction<Real>("h", _infile("thickness", 0.002));
//        MAST::Solid2DSectionElementPropertyCard *p = new MAST::Solid2DSectionElementPropertyCard(2);
//        p->add(*_h_stiff);
//        p->set_material(mat);
//        p->set_diagonal_mass_matrix(false);
//        p->set_strain(MAST::VON_KARMAN_STRAIN);
//        
//        _elem_properties[1] = p;
//        
//        _parameters[0] = _h->ptr(); // set thickness as a modifiable parameter
//        _parameters[1] = _h_stiff->ptr(); // set thickness as a modifiable parameter
//        
//        _structural_assembly->add_parameter(*_h);
//        _structural_assembly->add_parameter(*_h_stiff);
    }
    

    MAST::Solid2DSectionElementPropertyCard &prop2d =
    dynamic_cast<MAST::Solid2DSectionElementPropertyCard&> (*_elem_properties[0]);
    prop2d.set_material(mat);
    prop2d.set_diagonal_mass_matrix(false);
    prop2d.add(*_h);
    //prop2d.add(*_prestress); // no prestress for stiffener
    
    prop2d.set_strain(MAST::VON_KARMAN_STRAIN);

    _static_structural_assembly->set_property_for_subdomain(0, prop2d);
    _eigen_structural_assembly->set_property_for_subdomain(0, prop2d);
    for (unsigned int i=1; i<=_n_stiff; i++) {
        _static_structural_assembly->set_property_for_subdomain(i, *_elem_properties[i]);
        _eigen_structural_assembly->set_property_for_subdomain(i, *_elem_properties[i]);
    }
    
    // create the function to calculate weight
    _weight = new MAST::Weight(*_mesh, *_static_structural_assembly);
    
#endif // LIBMESH_USE_COMPLEX_NUMBERS
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
        Nemesis_IO(*_mesh).write_equation_systems(file_name.str(),
                                                  *_eq_systems);
    }
    
    MAST::FunctionEvaluation::output(iter, x, obj, fval);
}


#endif // __MAST_sizing_optimization_h__
