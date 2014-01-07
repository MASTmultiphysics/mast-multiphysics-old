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
#include "PropertyCards/element_property_card_2D.h"

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


namespace MAST {
    
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
        _eq_systems(NULL),
        _system(NULL),
        _structural_assembly(NULL),
        _mesh(NULL),
        _press(NULL),
        _zero_function(NULL)
        {
            _init();
        }
        
        virtual ~SizingOptimization() {
            
            delete _structural_assembly;
            
            delete _eq_systems;
            
            delete _mesh;
            
            for (unsigned int i=0; i<_bc.size(); i++)
                delete _bc[i];
            
            delete _zero_function;
            
            
            for (unsigned int i=0; i<_elem_properties.size(); i++)
                delete _elem_properties[i];
            
            for (unsigned int i=0; i<_materials.size(); i++)
                delete _materials[i];
            
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
        
    protected:
        
        
        void _init();
        
        
        LibMeshInit& _libmesh_init;
        
        GetPot& _infile;
        
        unsigned int _n_stiff, _n_eig;
        
        EquationSystems* _eq_systems;
        
        CondensedEigenSystem* _system;
        
        MAST::StructuralSystemAssembly* _structural_assembly;
        
        UnstructuredMesh* _mesh;
        
        ConstFunction<Real>* _press;
        
        ZeroFunction<Real>* _zero_function;
        
        ParameterVector _parameters;
        
        std::vector<MAST::DisplacementBoundaryCondition*> _bc;

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
    std::fill(x.begin(), x.end(), 0.003);           // start with a solid material
    xmin.resize(_n_vars);
    std::fill(xmin.begin(), xmin.end(), 2.0e-3); // lower limit is a small value.
    xmax.resize(_n_vars);
    std::fill(xmax.begin(), xmax.end(), 0.02);     // upper limit is 1.
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
    *_parameters[0] = dvars[0];
    *_parameters[1] = dvars[1];

    MeshBase::const_element_iterator
    eit  = _mesh->active_local_elements_begin(),
    eend = _mesh->active_local_elements_end();
    
    obj = 0.;
    
    if (eval_obj_grad)
        std::fill(obj_grad.begin(), obj_grad.end(), 0.);
    
    for ( ; eit != eend; eit++ ) {
        const Elem* e = *eit;
        const MAST::ElementPropertyCardBase& prop = _structural_assembly->get_property_card(*e);
        const MAST::MaterialPropertyCardBase& mat = prop.get_material();

        _structural_assembly->get_property_card(0);
        _structural_assembly->get_property_card(1);
        _structural_assembly->get_property_card(2);

        obj += e->volume() * prop.get<Real>("h")() * mat.get<Real>("rho")();
        
        if (eval_obj_grad)
            for (unsigned int i=0; i<_n_vars; i++)
                if (prop.depends_on(_parameters[i]))
                    obj_grad[i] += e->volume() * mat.get<Real>("rho")();
    }
    
    // now solve the system
    libMesh::out << "New Eval" << std::endl;
    _system->solve();

    // the number of converged eigenpairs could be different from the number asked
    unsigned int n_required = std::min(_n_eig, _system->get_n_converged());
    
    std::pair<Real, Real> val;
    Complex eigval;
    std::fill(fvals.begin(), fvals.end(), 0.);
    if (_eq_systems->parameters.get<bool>("if_exchange_AB_matrices"))
        // the total number of constraints is _n_eig, but only n_required are usable
        for (unsigned int i=0; i<n_required; i++) {
            val = _system->get_eigenpair(i);
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
        // grad k = dfi/dxj  ,  where k = j*NDV + i
        _system->sensitivity_solve(_parameters, grad_vals);
        // now correct this for the fact that the matrices were exchanged
        if (_eq_systems->parameters.get<bool>("if_exchange_AB_matrices"))
            for (unsigned int j=0; j<_n_vars; j++)
                for (unsigned int i=0; i<n_required; i++) {
                    val = _system->get_eigenpair(i);
                    grads[j*_n_eig+i]  = grad_vals[j*_system->get_n_converged()+i];
                    grads[j*_n_eig+i] /= -pow(val.first, 2);
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
    _n_vars = 2;
    _n_eq = 0;
    _n_ineq = _n_eig;
    _max_iters = 10000;

    // now initialize the mesh
    bool beam_stiff = _infile("beam_stiffeners", false);
    MAST::StiffenedPanel().init(divs, *_mesh, elem_type, beam_stiff);

    // Print information about the mesh to the screen.
    _mesh->print_info();
    
    // Create an equation systems object.
    _eq_systems = new EquationSystems(*_mesh);
    _eq_systems->parameters.set<GetPot*>("input_file") = &_infile;
    
    // Declare the system
    _system = &_eq_systems->add_system<CondensedEigenSystem> ("StructuralSystem");
    
    unsigned int o = _infile("fe_order", 1);
    std::string fe_family = _infile("fe_family", std::string("LAGRANGE"));
    FEFamily fefamily = Utility::string_to_enum<FEFamily>(fe_family);
    
    std::map<std::string, unsigned int> var_id;
    var_id["ux"] = _system->add_variable ( "ux", static_cast<Order>(o), fefamily);
    var_id["uy"] = _system->add_variable ( "uy", static_cast<Order>(o), fefamily);
    var_id["uz"] = _system->add_variable ( "uz", static_cast<Order>(o), fefamily);
    var_id["tx"] = _system->add_variable ( "tx", static_cast<Order>(o), fefamily);
    var_id["ty"] = _system->add_variable ( "ty", static_cast<Order>(o), fefamily);
    var_id["tz"] = _system->add_variable ( "tz", static_cast<Order>(o), fefamily);
    
    _structural_assembly = new MAST::StructuralSystemAssembly(*_system,
                                                              MAST::BUCKLING,
                                                              _infile);
    
    _system->attach_assemble_object(*_structural_assembly);
    _system->attach_eigenproblem_sensitivity_assemble_object(*_structural_assembly);
    
    
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
        _bc[counter] = new MAST::DisplacementBoundaryCondition;
        _bc[counter]->init(it->first, it->second);
        _structural_assembly->add_side_load(it->first, *_bc[counter]);
        counter++;
    }

    _system->set_eigenproblem_type(GHEP);
    _system->eigen_solver->set_position_of_spectrum(LARGEST_MAGNITUDE);

    _eq_systems->init ();
    
    
    // Print information about the system to the screen.
    _eq_systems->print_info();

    // Pass the Dirichlet dof IDs to the CondensedEigenSystem
    std::set<unsigned int> dirichlet_dof_ids;
    _eq_systems->parameters.set<bool>("if_exchange_AB_matrices") = true;
    _eq_systems->parameters.set<unsigned int>("eigenpairs")    = _n_eig;
    _eq_systems->parameters.set<unsigned int>("basis vectors") = _n_eig*3;
    _structural_assembly->get_dirichlet_dofs(dirichlet_dof_ids);
    _system->initialize_condensed_dofs(dirichlet_dof_ids);

    // element and material properties
    _materials.resize(1);
    _elem_properties.resize(_n_vars);
    _parameters.resize(_n_vars);
    
    _materials[0] = new MAST::IsotropicMaterialPropertyCard(0);
    _elem_properties[0] = new MAST::Solid2DSectionElementPropertyCard(1);
    _elem_properties[1] = new MAST::Solid2DSectionElementPropertyCard(2);
    
    MAST::MaterialPropertyCardBase& mat = *_materials[0];
    MAST::Solid2DSectionElementPropertyCard
    &prop2d = dynamic_cast<MAST::Solid2DSectionElementPropertyCard&> (*_elem_properties[0]),
    &prop2d_stiff = dynamic_cast<MAST::Solid2DSectionElementPropertyCard&> (*_elem_properties[1]);
    
    MAST::FunctionValue<Real>& E = mat.add<Real>("E", MAST::CONSTANT_FUNCTION),
    &nu = mat.add<Real>("nu", MAST::CONSTANT_FUNCTION),
    &rho = mat.add<Real>("rho", MAST::CONSTANT_FUNCTION),
    &kappa = mat.add<Real>("kappa", MAST::CONSTANT_FUNCTION),
    &h =  prop2d.add<Real>("h", MAST::CONSTANT_FUNCTION),
    &h_stiff =  prop2d_stiff.add<Real>("h", MAST::CONSTANT_FUNCTION);
    
    
    E     = _infile("youngs_modulus", 72.e9);
    nu    = _infile("poisson_ratio", 0.33);
    rho   = _infile("material_density", 2700.);
    kappa = _infile("shear_corr_factor", 5./6.);
    h     = _infile("thickness", 0.002);
    h_stiff = _infile("thickness", 0.002);
    
    DenseVector<Real> prestress; prestress.resize(6);
    prestress(0) = -1.31345e6;
    
    prop2d.set_material(mat); prop2d_stiff.set_material(mat);
    prop2d.set_diagonal_mass_matrix(false); prop2d_stiff.set_diagonal_mass_matrix(false);
    prop2d.prestress(prestress); // no prestress for stiffener
    
    prop2d.set_strain(MAST::VON_KARMAN_STRAIN); prop2d_stiff.set_strain(MAST::VON_KARMAN_STRAIN);
    _parameters[0] = h.ptr(); // set thickness as a modifiable parameter
    _parameters[1] = h_stiff.ptr(); // set thickness as a modifiable parameter
    
    
    _structural_assembly->set_property_for_subdomain(0, prop2d);
    for (unsigned int i=1; i<=_n_stiff; i++)
        _structural_assembly->set_property_for_subdomain(i, prop2d_stiff);
    _structural_assembly->add_parameter(h.ptr(), &h);
    _structural_assembly->add_parameter(h_stiff.ptr(), &h_stiff);
    
#endif // LIBMESH_USE_COMPLEX_NUMBERS
}



#endif // __MAST_sizing_optimization_h__
