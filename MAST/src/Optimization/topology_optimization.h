//
//  topology_optimization.h
//  MAST
//
//  Created by Manav Bhatia on 11/14/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_topology_optimization_h__
#define __MAST_topology_optimization_h__

// MAST includes
#include "Optimization/optimization_interface.h"
#include "StructuralElems/structural_system_assembly.h"

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


namespace MAST {
    
    class TopologyOptimization:
    public MAST::FunctionEvaluation {

    public:
        
        TopologyOptimization(LibMeshInit& init,
                             GetPot& input):
        MAST::FunctionEvaluation(),
        _libmesh_init(init),
        _infile(input),
        _eq_systems(NULL),
        _system(NULL),
        _structural_assembly(NULL),
        _mesh(NULL),
        _press(NULL),
        _zero_function(NULL),
        _bc(NULL),
        _penalty(2.),
        E0(72.0e9),
        V0(0.3)
        {
            _init();
        }

        virtual ~TopologyOptimization() {
            
            delete _structural_assembly;
            
            delete _eq_systems;
            
            delete _mesh;
            
            delete _bc;
            
            delete _press;
            
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
        
        EquationSystems* _eq_systems;
        
        NonlinearImplicitSystem* _system;
        
        MAST::StructuralSystemAssembly* _structural_assembly;
        
        UnstructuredMesh* _mesh;
        
        ConstFunction<Real>* _press;
        
        ZeroFunction<Real>* _zero_function;
        
        MAST::BoundaryCondition* _bc;
        
        Real _penalty, E0, V0;

        ParameterVector _parameters;
        
        std::vector<MAST::MaterialPropertyCardBase*> _materials;
        
        std::vector<MAST::ElementPropertyCardBase*> _elem_properties;
        
        // relative volume of each element
        std::vector<Real> _elem_vol;
    };
}



inline
void
MAST::TopologyOptimization::init_dvar(std::vector<Real>& x,
                                      std::vector<Real>& xmin,
                                      std::vector<Real>& xmax) {
    // one DV for each element
    x.resize   (_n_vars);
    std::fill(x.begin(), x.end(), 1.);           // start with a solid material
    xmin.resize(_n_vars);
    std::fill(xmin.begin(), xmin.end(), 1.0e-3); // lower limit is a small value.
    xmax.resize(_n_vars);
    std::fill(xmax.begin(), xmax.end(), 1.);     // upper limit is 1.
}


inline
void
MAST::TopologyOptimization::evaluate(const std::vector<Real>& dvars,
                                     Real& obj,
                                     bool eval_obj_grad,
                                     std::vector<Real>& obj_grad,
                                     std::vector<Real>& fvals,
                                     std::vector<bool>& eval_grads,
                                     std::vector<Real>& grads) {
 
    libmesh_assert_equal_to(dvars.size(), _mesh->n_elem());
        
    // set the parameter values equal to the DV value
    fvals[0] = -V0;
    for (unsigned int i=0; i<dvars.size(); i++) {
        *_parameters[i] = E0 * pow(dvars[i], _penalty); // Ei = E0 * x^p
        fvals[0] += dvars[i] * _elem_vol[i]; // constraint:  xi vi - V <= 0
    }
    
    
    // now solve the system
    _system->solve();
    _system->sensitivity_solve(_parameters);
    
    // now calculate the objective function
    std::auto_ptr<NumericVector<Number> >
    vec(NumericVector<Number>::build(_system->comm()).release());
    vec->init(*_system->solution);
    _system->matrix->vector_mult(*vec, *_system->solution);
    obj = -0.5 * vec->dot(*_system->solution); // negative, since J = -K
    
    if (eval_obj_grad) { // objective gradients
        // iterate over each dv and calculate the sensitivity
        for (unsigned int i=0; i<dvars.size(); i++) {
            NumericVector<Number>& sens = _system->get_sensitivity_solution(i);
            // df0/dxi = df0/dEi * dEi/dxi
            // dEi/dxi = E0 * p * x^(p-1)
            obj_grad[i] = vec->dot(sens) * E0 * _penalty * pow(dvars[i], _penalty-1.);
        }
    }
    
    // now calculate the constraint gradient
    // single volume constraint
    if (eval_grads[0])
        for (unsigned int i=0; i<dvars.size(); i++)
            grads[i] = _elem_vol[i];
}


inline
void
MAST::TopologyOptimization::_init() {

    _mesh = new SerialMesh(_libmesh_init.comm());
    _mesh->set_mesh_dimension(2);
    
    MeshTools::Generation::build_square (*_mesh, 10, 10, 0., 1., 0., .1, TRI3);
    
    _mesh->prepare_for_use();
    
    // Print information about the mesh to the screen.
    _mesh->print_info();
    
    // Create an equation systems object.
    _eq_systems = new EquationSystems(*_mesh);
    _eq_systems->parameters.set<GetPot*>("input_file") = &_infile;
    
    // Declare the system
    _system = &_eq_systems->add_system<NonlinearImplicitSystem> ("StructuralSystem");
    
    unsigned int o = _infile("fe_order", 1);
    std::string fe_family = _infile("fe_family", std::string("LAGRANGE"));
    FEFamily fefamily = Utility::string_to_enum<FEFamily>(fe_family);
    
    _system->add_variable ( "ux", static_cast<Order>(o), fefamily);
    _system->add_variable ( "uy", static_cast<Order>(o), fefamily);
    _system->add_variable ( "uz", static_cast<Order>(o), fefamily);
    _system->add_variable ( "tx", static_cast<Order>(o), fefamily);
    _system->add_variable ( "ty", static_cast<Order>(o), fefamily);
    _system->add_variable ( "tz", static_cast<Order>(o), fefamily);
    
    _structural_assembly = new MAST::StructuralSystemAssembly(*_system,
                                                              MAST::STATIC,
                                                              _infile);
    
    _press = new ConstFunction<Real>(1.e2);
    _bc = new MAST::BoundaryCondition(MAST::SURFACE_PRESSURE);
    _bc->set_function(*_press);
    _structural_assembly->add_side_load(2, *_bc);
    
    
    _system->attach_assemble_object(*_structural_assembly);
    _system->attach_sensitivity_assemble_object(*_structural_assembly);
    
    // Pass the Dirichlet dof IDs to the CondensedEigenSystem
    
    // apply the boundary conditions
    _zero_function = new ZeroFunction<Real>;
    std::vector<unsigned int> vars(6);
    for (unsigned int i=0; i<6; i++)
        vars[i] = i;
    std::set<boundary_id_type> dirichlet_boundary;
    //dirichlet_boundary.insert(0); // bottom
    dirichlet_boundary.insert(1); // right
    dirichlet_boundary.insert(3); // left
    _system->get_dof_map().add_dirichlet_boundary(DirichletBoundary(dirichlet_boundary, vars,
                                                                    _zero_function));
    _eq_systems->init ();
    
    
    // Print information about the system to the screen.
    _eq_systems->print_info();
    
    
    unsigned int n_elems = _mesh->n_active_local_elem();
    _n_vars = n_elems;
    _n_eq = 0;
    _n_ineq = 1;
    _max_iters = 100;
    
    _parameters.resize(n_elems);
    _elem_properties.resize(n_elems);
    _materials.resize(n_elems);
    _elem_vol.resize(n_elems);
    
    // now iterate over each element and create property cards for each element
    MeshBase::const_element_iterator
    eit  = _mesh->active_local_elements_begin(),
    eend = _mesh->active_local_elements_end();

    Real total_vol = 0.;
    
    unsigned int counter = 0;
    for ( ; eit != eend; eit++ ) {
        const Elem* e = *eit;
        _elem_vol[counter] = e->volume();
        total_vol += _elem_vol[counter];
        
        // create material properties for the element
        MAST::IsotropicMaterialPropertyCard* mat = new MAST::IsotropicMaterialPropertyCard;
        MAST::Solid2DSectionElementPropertyCard* prop2d = new MAST::Solid2DSectionElementPropertyCard;
        _materials[counter] = mat;
        _elem_properties[counter] = prop2d;
        
        MAST::FunctionValue<Real>& E = mat->add<Real>("E", MAST::CONSTANT_FUNCTION),
        &nu = mat->add<Real>("nu", MAST::CONSTANT_FUNCTION),
        &rho = mat->add<Real>("rho", MAST::CONSTANT_FUNCTION),
        &kappa = mat->add<Real>("kappa", MAST::CONSTANT_FUNCTION),
        &h =  prop2d->add<Real>("h", MAST::CONSTANT_FUNCTION);
        E  = E0;
        nu = 0.33;
        rho = 2700.;
        kappa = 5./6.;
        h  = 0.002;

        prop2d->set_material(*mat);
        _structural_assembly->set_property_for_elem(*e, *prop2d);
        _structural_assembly->add_parameter(E.ptr(), &E);
        _parameters[counter] = E.ptr(); // set Young's modulus as a modifiable parameter
        
        // increment counter for storing objects for next element
        counter++;
    }
    
    // now calculate relative element volumes
    for (unsigned int i=0; i<_elem_vol.size(); i++)
        _elem_vol[i] /= total_vol;
}


#endif // __MAST_topology_optimization_h__
