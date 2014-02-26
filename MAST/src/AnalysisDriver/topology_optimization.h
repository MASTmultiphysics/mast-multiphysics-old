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
#include "PropertyCards/isotropic_material_property_card.h"
#include "PropertyCards/solid_2d_section_element_property_card.h"
#include "Numerics/constant_function.h"

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


namespace MAST {
    
    class TopologyOptimization:
    public MAST::FunctionEvaluation {

    public:
        
        TopologyOptimization(LibMeshInit& init,
                             GetPot& input,
                             std::ostream& output):
        MAST::FunctionEvaluation(output),
        _libmesh_init(init),
        _infile(input),
        _eq_systems(NULL),
        _system(NULL),
        _rho_system(NULL),
        _structural_assembly(NULL),
        _compliance(NULL),
        _mesh(NULL),
        _press(NULL),
        _zero_function(NULL),
        _bc(NULL),
        _penalty(3.),
        E0(72.0e9),
        V0(0.3)
        {
            _init();
        }

        virtual ~TopologyOptimization() {
            
            delete _structural_assembly;
            
            delete _compliance;
            
            delete _eq_systems;
            
            delete _mesh;
            
            delete _bc;
            
            delete _press;
            
            delete _zero_function;
            
            
            for (unsigned int i=0; i<_elem_properties.size(); i++)
                delete _elem_properties[i];
            
            for (unsigned int i=0; i<_materials.size(); i++)
                delete _materials[i];

            for (unsigned int i=0; i<_E.size(); i++)
                delete _E[i];
            
            delete _nu;
            delete _rho;
            delete _kappa;
            delete _h;
            delete _h_stiff;
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

        class Compliance:
        public System::QOI,
        public System::QOIDerivative {
        public:
            Compliance(NonlinearImplicitSystem& sys):
            System::QOI(),
            System::QOIDerivative(),
            _system(sys)
            { }
            
            virtual void qoi(const QoISet& qoi_indices);
            
            virtual void qoi_derivative (const QoISet& qoi_indices);
            
        protected:
            NonlinearImplicitSystem& _system;
        };
        

        void _init();
        

        LibMeshInit& _libmesh_init;
        
        GetPot& _infile;
        
        EquationSystems* _eq_systems;
        
        NonlinearImplicitSystem* _system;

        System* _rho_system;

        MAST::StructuralSystemAssembly* _structural_assembly;
        
        MAST::TopologyOptimization::Compliance* _compliance;
        
        UnstructuredMesh* _mesh;
        
        ConstantFunction<Real>* _press;
        
        ZeroFunction<Real>* _zero_function;
        
        MAST::BoundaryCondition* _bc;
        
        Real _penalty, E0, V0;

        std::vector<MAST::ConstantFunction<Real>*> _E;
        MAST::ConstantFunction<Real> *_nu, *_rho, *_kappa, *_h, *_h_stiff;

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
    std::fill(x.begin(), x.end(), 0.3);           // start with a solid material
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
 #ifndef LIBMESH_USE_COMPLEX_NUMBERS
    
    libmesh_assert_equal_to(dvars.size(), _mesh->n_elem());
        
    // set the parameter values equal to the DV value
    fvals[0] = -V0;
    for (unsigned int i=0; i<dvars.size(); i++) {
        *_parameters[i] = E0 * pow(dvars[i], _penalty); // Ei = E0 * x^p
        fvals[0] += dvars[i] * _elem_vol[i]; // constraint:  xi vi - V <= 0
    }
    
    
    // now solve the system
    std::cout << "New Eval" << std::endl;
    QoISet qoi_set(*_system);
    qoi_set.add_index(0);
    _system->solution->zero();
    _system->solve();
    _system->assemble_qoi(qoi_set);
    
    obj = _system->qoi[0]; // compliance
    
    if (eval_obj_grad) { // objective gradients
        SensitivityData sens;
        _system->set_adjoint_already_solved(false);
        _system->adjoint_qoi_parameter_sensitivity(qoi_set, _parameters, sens);
        // df0/dxi = df0/dEi * dEi/dxi
        // dEi/dxi = E0 * p * x^(p-1)
        for (unsigned int i=0; i<dvars.size(); i++)
            obj_grad[i] = sens[0][i] * E0 * _penalty * pow(dvars[i], _penalty-1.);
        _system->set_adjoint_already_solved(true);
    }
    
    // now calculate the constraint gradient
    // single volume constraint
    if (eval_grads[0])
        for (unsigned int i=0; i<dvars.size(); i++)
            grads[i] = _elem_vol[i];

#endif // LIBMESH_USE_COMPLEX_NUMBERS
}


inline
void
MAST::TopologyOptimization::_init() {

#ifndef LIBMESH_USE_COMPLEX_NUMBERS
    _mesh = new SerialMesh(_libmesh_init.comm());
    _mesh->set_mesh_dimension(2);
    
    MeshTools::Generation::build_square (*_mesh, 50, 20, 0., 1., 0., .1, QUAD4);
    
    _mesh->prepare_for_use();
    
    // Print information about the mesh to the screen.
    _mesh->print_info();
    
    // Create an equation systems object.
    _eq_systems = new EquationSystems(*_mesh);
    _eq_systems->parameters.set<GetPot*>("input_file") = &_infile;
    
    // Declare the system
    _system = &_eq_systems->add_system<NonlinearImplicitSystem> ("StructuralSystem");
    _rho_system = &_eq_systems->add_system<System> ("RhoSystem");
    _rho_system->add_variable("rho", FEType(CONSTANT, MONOMIAL));

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
    _compliance = new MAST::TopologyOptimization::Compliance(*_system);
    
    _press = new MAST::ConstantFunction<Real>("pressure",1.e5);
    _bc = new MAST::BoundaryCondition(MAST::SURFACE_PRESSURE);
    _bc->set_function(*_press);
    _structural_assembly->add_side_load(2, *_bc);
    
    
    _system->attach_assemble_object(*_structural_assembly);
    _system->attach_sensitivity_assemble_object(*_structural_assembly);
    _system->attach_QOI_object(*_compliance);
    _system->attach_QOI_derivative_object(*_compliance);
    _system->qoi.resize(1);

    
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
    _max_iters = 10000;
    
    _parameters.resize(n_elems);
    _E.resize(n_elems);
    _elem_properties.resize(n_elems);
    _materials.resize(n_elems);
    _elem_vol.resize(n_elems);
    
    // now iterate over each element and create property cards for each element
    MeshBase::const_element_iterator
    eit  = _mesh->active_local_elements_begin(),
    eend = _mesh->active_local_elements_end();

    Real total_vol = 0.;
    
    _nu = new MAST::ConstantFunction<Real>("nu", _infile("poisson_ratio", 0.33)),
    _rho = new MAST::ConstantFunction<Real>("rho", _infile("material_density", 2700.)),
    _kappa = new MAST::ConstantFunction<Real>("kappa", _infile("shear_corr_factor", 5./6.)),
    _h = new MAST::ConstantFunction<Real>("h", _infile("thickness", 0.002)),
    *_nu = 0.33;
    *_rho = 2700.;
    *_kappa = 5./6.;
    *_h  = 0.002;

    unsigned int counter = 0;
    for ( ; eit != eend; eit++ ) {
        const Elem* e = *eit;
        _elem_vol[counter] = e->volume();
        total_vol += _elem_vol[counter];
        
        // create material properties for the element
        MAST::IsotropicMaterialPropertyCard* mat = new MAST::IsotropicMaterialPropertyCard(0);
        MAST::Solid2DSectionElementPropertyCard* prop2d = new MAST::Solid2DSectionElementPropertyCard(0);
        _materials[counter] = mat;
        _elem_properties[counter] = prop2d;
        
        _E[counter] = new MAST::ConstantFunction<Real>("E", E0);
        
        prop2d->set_material(*mat);
        _structural_assembly->set_property_for_subdomain(0, *prop2d);
        _structural_assembly->add_parameter(*_E[counter]);
        _parameters[counter] = _E[counter]->ptr(); // set Young's modulus as a modifiable parameter
        
        // increment counter for storing objects for next element
        counter++;
    }
    
    // now calculate relative element volumes
    for (unsigned int i=0; i<_elem_vol.size(); i++)
        _elem_vol[i] /= total_vol;
    
#endif // LIBMESH_USE_COMPLEX_NUMBERS
}


inline void
MAST::TopologyOptimization::output(unsigned int iter,
                                   const std::vector<Real>& x,
                                   Real obj,
                                   const std::vector<Real>& fval) const  {
    for (unsigned int i=0; i<x.size(); i++)
        _rho_system->solution->set(_mesh->elem(i)->dof_number(_rho_system->number(), 0, 0),
                                   x[i]);
    
    
    _rho_system->solution->close();
    _rho_system->update();
    
    Nemesis_IO(*_mesh).write_equation_systems("out.exo", *_eq_systems);
    
    MAST::FunctionEvaluation::output(iter, x, obj, fval);
}


inline void
MAST::TopologyOptimization::Compliance::qoi(const QoISet& qoi_indices){
#ifndef LIBMESH_USE_COMPLEX_NUMBERS
    if (qoi_indices.has_index(0)) {
        
        std::auto_ptr<NumericVector<Number> >
        vec(NumericVector<Number>::build(_system.comm()).release());
        vec->init(*_system.solution);
        
        _system.matrix->vector_mult(*vec, *_system.solution);
        _system.qoi[0] = -1 * vec->dot(*_system.solution); // negative, since J = -K
    }
#endif // LIBMESH_USE_COMPLEX_NUMBERS
}


inline void
MAST::TopologyOptimization::Compliance::qoi_derivative (const QoISet& qoi_indices) {
    if (qoi_indices.has_index(0))
        _system.matrix->vector_mult(_system.get_adjoint_rhs(), *_system.solution);
    _system.get_adjoint_rhs().scale(-1.);
}


#endif // __MAST_topology_optimization_h__
