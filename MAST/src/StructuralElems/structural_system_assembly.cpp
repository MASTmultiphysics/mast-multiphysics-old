//
//  structural_system_assembly.h
//  MAST
//
//  Created by Manav Bhatia on 12/16/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

// libMesh include files
#include "libmesh/nonlinear_solver.h"
#include "libmesh/condensed_eigen_system.h"
#include "libmesh/fe_interface.h"
#include "libmesh/parameter_vector.h"
#include "libmesh/dof_map.h"

// MAST includes
#include "StructuralElems/structural_system_assembly.h"
#include "StructuralElems/structural_elem_base.h"
#include "PropertyCards/element_property_card_base.h"
#include "BoundaryConditions/boundary_condition.h"
#include "BoundaryConditions/displacement_boundary_condition.h"


MAST::StructuralSystemAssembly::StructuralSystemAssembly(System& sys,
                                                         MAST::StructuralAnalysisType t,
                                                         GetPot& infile):

_system(sys),
_analysis_type(t),
_infile(infile),
_static_sol_system(NULL)
{
    // depending on the analysis type, forward to the appropriate function
    switch (_analysis_type) {
        case MAST::STATIC:
            dynamic_cast<NonlinearImplicitSystem&>(sys).
            nonlinear_solver->residual_and_jacobian_object =
            dynamic_cast<NonlinearImplicitSystem::ComputeResidualandJacobian*>(this);
            break;
            
        case MAST::MODAL:
        case MAST::BUCKLING:
            // nothing to be done here for now
            break;
            
        case MAST::DYNAMIC:
        default:
            libmesh_error();
            break;
    }
}



void
MAST::StructuralSystemAssembly::clear_loads() {
    _side_bc_map.clear();
    _vol_bc_map.clear();
    _static_sol_system = NULL;
}


void
MAST::StructuralSystemAssembly::add_side_load(boundary_id_type bid,
                                              MAST::BoundaryCondition& load) {
    // make sure that this boundary and load haven't already been applied
    std::pair<std::multimap<boundary_id_type, MAST::BoundaryCondition*>::const_iterator,
    std::multimap<boundary_id_type, MAST::BoundaryCondition*>::const_iterator> it =
    _side_bc_map.equal_range(bid);
    
    for ( ; it.first != it.second; it.first++) {
        libmesh_assert(it.first->second != &load);
        // only one displacement boundary condition is allowed per boundary
        if (load.type() == MAST::DISPLACEMENT_DIRICHLET)
            libmesh_assert(it.first->second->type() != MAST::DISPLACEMENT_DIRICHLET);
    }
    
    // displacement boundary condition needs to be hadled separately
    if (load.type() == MAST::DISPLACEMENT_DIRICHLET) {
        
        // get the Dirichlet boundary condition object
        DirichletBoundary& dirichlet_b =
        (dynamic_cast<MAST::DisplacementDirichletBoundaryCondition&>(load)).dirichlet_boundary();
        
        // add an entry for each boundary of this dirichlet object.
        for (std::set<boundary_id_type>::const_iterator it = dirichlet_b.b.begin();
             it != dirichlet_b.b.end(); it++)
            _side_bc_map.insert(std::multimap<boundary_id_type, MAST::BoundaryCondition*>::value_type
                                (*it, &load));
        
        // now add this to the dof_map for this system
        _system.get_dof_map().add_dirichlet_boundary(dirichlet_b);
    }
    else
        _side_bc_map.insert(std::multimap<boundary_id_type, MAST::BoundaryCondition*>::value_type
                            (bid, &load));
}



void
MAST::StructuralSystemAssembly::add_volume_load(subdomain_id_type bid,
                                                MAST::BoundaryCondition& load) {
    std::pair<std::multimap<subdomain_id_type, MAST::BoundaryCondition*>::const_iterator,
    std::multimap<subdomain_id_type, MAST::BoundaryCondition*>::const_iterator> it =
    _vol_bc_map.equal_range(bid);
    
    for ( ; it.first != it.second; it.first++)
        libmesh_assert(it.first->second != &load);
    
    _vol_bc_map.insert(std::multimap<boundary_id_type, MAST::BoundaryCondition*>::value_type
                       (bid, &load));
}



void
MAST::StructuralSystemAssembly::set_property_for_subdomain(const subdomain_id_type sid,
                                                           const MAST::ElementPropertyCardBase& prop) {
    std::map<subdomain_id_type, const MAST::ElementPropertyCardBase*>::const_iterator
    elem_p_it = _element_property.find(sid);
    libmesh_assert(elem_p_it == _element_property.end());

    _element_property[sid] = &prop;
}



const MAST::ElementPropertyCardBase&
MAST::StructuralSystemAssembly::get_property_card(const unsigned int i) const {
    
    std::map<subdomain_id_type, const MAST::ElementPropertyCardBase*>::const_iterator
    elem_p_it = _element_property.find(i);
    libmesh_assert(elem_p_it != _element_property.end());
    
    return *elem_p_it->second;
}



const MAST::ElementPropertyCardBase&
MAST::StructuralSystemAssembly::get_property_card(const Elem& elem) const {
    
    std::map<subdomain_id_type, const MAST::ElementPropertyCardBase*>::const_iterator
    elem_p_it = _element_property.find(elem.subdomain_id());
    libmesh_assert(elem_p_it != _element_property.end());
    
    return *elem_p_it->second;
}


void
MAST::StructuralSystemAssembly::add_parameter(MAST::ConstantFunction<Real>& f) {
    Real* par = f.ptr();
    // make sure it does not already exist in the map
    libmesh_assert(!_parameter_map.count(par));
    
    // now add this to the map
    bool insert_success = _parameter_map.insert
    (std::map<Real*, MAST::FieldFunctionBase*>::value_type(par, &f)).second;
    
    libmesh_assert(insert_success);
}



const MAST::FieldFunctionBase*
MAST::StructuralSystemAssembly::get_parameter(Real* par) const {
    // make sure valid values are given
    libmesh_assert(par);
    
    std::map<Real*, const MAST::FieldFunctionBase*>::const_iterator
    it = _parameter_map.find(par);
    
    // make sure it does not already exist in the map
    libmesh_assert(it != _parameter_map.end());
    
    return it->second;
}



void
MAST::StructuralSystemAssembly::residual_and_jacobian (const NumericVector<Number>& X,
                                                       NumericVector<Number>* R,
                                                       SparseMatrix<Number>*  J,
                                                       NonlinearImplicitSystem& S) {
    
    if (R) R->zero();
    if (J) J->zero();
    
    switch (_analysis_type) {
        case MAST::STATIC:
        case MAST::DYNAMIC:
            libmesh_assert(_system.system_type() == "NonlinearImplicit");
            _assemble_residual_and_jacobian(X, R, J, S, NULL);
            break;
            
        default:
            libmesh_error();
            break;
    }
    
    if (R) R->close();
    if (J) J->close();
}




bool
MAST::StructuralSystemAssembly::sensitivity_assemble (const ParameterVector& params,
                                                      const unsigned int i,
                                                      NumericVector<Number>& sensitivity_rhs) {
#ifndef LIBMESH_USE_COMPLEX_NUMBERS
    
    const MAST::FieldFunctionBase* f = this->get_parameter(params[i]);
    
    sensitivity_rhs.zero();
    
    switch (_analysis_type) {
        case MAST::STATIC:
        case MAST::DYNAMIC:
            libmesh_assert(_system.system_type() == "NonlinearImplicit");
            _assemble_residual_and_jacobian(*_system.solution, &sensitivity_rhs,
                                            NULL,
                                            dynamic_cast<NonlinearImplicitSystem&>(_system),
                                            f);
            break;
            
        default:
            libmesh_error();
            break;
    }
    
    sensitivity_rhs.close();
    // currently, all relevant parameter sensitivities are calculated
    return true;
    
#endif // LIBMESH_USE_COMPLEX_NUMBERS
}



void
MAST::StructuralSystemAssembly::assemble() {
    
    SparseMatrix<Number>&  matrix_A = *(dynamic_cast<EigenSystem&>(_system).matrix_A);
    SparseMatrix<Number>&  matrix_B = *(dynamic_cast<EigenSystem&>(_system).matrix_B);
    
    matrix_A.zero();
    matrix_B.zero();

    NumericVector<Number> *sol = NULL;
    
    // if a system is provided for getting the static solution, then get the vectors
    if (_static_sol_system)
        sol = _static_sol_system->solution.get();

    switch (_analysis_type) {
        case MAST::MODAL:
            _system.solution->zero();
            _assemble_matrices_for_modal_analysis(matrix_A,
                                                  matrix_B,
                                                  NULL,
                                                  sol, NULL);
            break;
            
        case MAST::BUCKLING:
            _system.solution->zero();
            _assemble_matrices_for_buckling_analysis(matrix_A,
                                                     matrix_B,
                                                     NULL,
                                                     sol, NULL);
            break;
            
        default:
            libmesh_error();
            break;
    }
    
    matrix_A.close();
    matrix_B.close();
}




bool
MAST::StructuralSystemAssembly::sensitivity_assemble (const ParameterVector& params,
                                                      const unsigned int i,
                                                      SparseMatrix<Number>* sensitivity_A,
                                                      SparseMatrix<Number>* sensitivity_B) {
#ifndef LIBMESH_USE_COMPLEX_NUMBERS
    
    const MAST::FieldFunctionBase* f = this->get_parameter(params[i]);
    
    sensitivity_A->zero();
    sensitivity_B->zero();
    
    NumericVector<Number> *sol = NULL, *sol_sens = NULL;
    
    // if a system is provided for getting the static solution, then get the vectors
    if (_static_sol_system) {
        sol = _static_sol_system->solution.get();
        sol_sens = &(_static_sol_system->get_sensitivity_solution(i));
    }
    
    switch (_analysis_type) {
        case MAST::MODAL:
            _assemble_matrices_for_modal_analysis(*sensitivity_A,
                                                  *sensitivity_B,
                                                  f,
                                                  sol, sol_sens);
            break;
            
        case MAST::BUCKLING:
            _assemble_matrices_for_buckling_analysis(*sensitivity_A,
                                                     *sensitivity_B,
                                                     f,
                                                     sol, sol_sens);
            break;
            
        default:
            libmesh_error();
            break;
    }
    
    sensitivity_A->close();
    sensitivity_B->close();

    // currently, all relevant parameter sensitivities are calculated
#endif // LIBMESH_USE_COMPLEX_NUMBERS
    return true;
}



void
MAST::StructuralSystemAssembly::_assemble_residual_and_jacobian (const NumericVector<Number>& X,
                                                                 NumericVector<Number>* R,
                                                                 SparseMatrix<Number>*  J,
                                                                 NonlinearImplicitSystem& S,
                                                                 const MAST::FieldFunctionBase* param) {
#ifndef LIBMESH_USE_COMPLEX_NUMBERS
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    DenseVector<Real> vec, sol;
    DenseMatrix<Real> mat;
    std::vector<dof_id_type> dof_indices;
    const DofMap& dof_map = _system.get_dof_map();
    std::auto_ptr<MAST::StructuralElementBase> structural_elem;
    
    AutoPtr<NumericVector<Number> > localized_solution =
    NumericVector<Number>::build(_system.comm());
    localized_solution->init(_system.n_dofs(), _system.n_local_dofs(),
                             _system.get_dof_map().get_send_list(),
                             false, GHOSTED);
    X.localize(*localized_solution, _system.get_dof_map().get_send_list());
    
    MeshBase::const_element_iterator       el     = _system.get_mesh().active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = _system.get_mesh().active_local_elements_end();
    
    for ( ; el != end_el; ++el) {
        
        const Elem* elem = *el;
        
        dof_map.dof_indices (elem, dof_indices);
        
        const MAST::ElementPropertyCardBase& p_card = this->get_property_card(*elem);
        
        // create the structural element for analysis
        structural_elem.reset(MAST::build_structural_element
                              (_system, *elem, p_card).release());
        
        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        sol.resize(ndofs);
        vec.resize(ndofs);
        mat.resize(ndofs, ndofs);
        
        for (unsigned int i=0; i<dof_indices.size(); i++)
            sol(i) = (*localized_solution)(dof_indices[i]);
        structural_elem->local_solution.resize(sol.size());
        structural_elem->transform_to_local_system(sol, structural_elem->local_solution);
        
        // now get the vector values
        if (!param) {
            structural_elem->internal_force(J!=NULL?true:false,
                                            vec, mat, false);
            structural_elem->prestress_force(J!=NULL?true:false,
                                             vec, mat);
            if (_analysis_type == MAST::DYNAMIC)
                structural_elem->inertial_force(J!=NULL?true:false,
                                                vec, mat);
            structural_elem->side_external_force(J!=NULL?true:false,
                                                 vec, mat,
                                                 _side_bc_map);
            structural_elem->volume_external_force(J!=NULL?true:false,
                                                   vec, mat,
                                                   _vol_bc_map);
        }
        else {
            structural_elem->sensitivity_param = param;
            structural_elem->internal_force_sensitivity(J!=NULL?true:false,
                                                        vec, mat, false);
            structural_elem->prestress_force_sensitivity(J!=NULL?true:false,
                                                         vec, mat);
            if (_analysis_type == MAST::DYNAMIC)
                structural_elem->inertial_force_sensitivity(J!=NULL?true:false,
                                                            vec, mat);
            structural_elem->side_external_force_sensitivity(J!=NULL?true:false,
                                                             vec, mat,
                                                             _side_bc_map);
            structural_elem->volume_external_force_sensitivity(J!=NULL?true:false,
                                                               vec, mat,
                                                               _vol_bc_map);
            // scale vector by -1 to account for the fact that the sensitivity
            // appears on RHS of the system with a -ve sign on it
            vec.scale(-1.);
        }
        
        if (R && J)
            _system.get_dof_map().constrain_element_matrix_and_vector(mat, vec, dof_indices);
        else if (R)
            _system.get_dof_map().constrain_element_vector(vec, dof_indices);
        else
            _system.get_dof_map().constrain_element_matrix(mat, dof_indices);
        
        // add to the global matrices
        if (R) R->add_vector(vec, dof_indices);
        if (J) J->add_matrix(mat, dof_indices);
    }
#endif // LIBMESH_USE_COMPLEX_NUMBERS
}



void
MAST::StructuralSystemAssembly::assemble_small_disturbance_aerodynamic_force (const NumericVector<Number>& X,
                                                                              NumericVector<Number>& F) {
    F.zero();
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    DenseVector<Number> vec, sol;
    DenseMatrix<Number> mat;
    std::vector<dof_id_type> dof_indices;
    const DofMap& dof_map = _system.get_dof_map();
    std::auto_ptr<MAST::StructuralElementBase> structural_elem;
    
    AutoPtr<NumericVector<Number> > localized_solution =
    NumericVector<Number>::build(_system.comm());
    localized_solution->init(_system.n_dofs(), _system.n_local_dofs(),
                             _system.get_dof_map().get_send_list(),
                             false, GHOSTED);
    X.localize(*localized_solution, _system.get_dof_map().get_send_list());
    
    MeshBase::const_element_iterator       el     = _system.get_mesh().active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = _system.get_mesh().active_local_elements_end();
    
    std::multimap<boundary_id_type, MAST::BoundaryCondition*> local_side_bc_map;
    std::multimap<subdomain_id_type, MAST::BoundaryCondition*> local_vol_bc_map;
    
    {
        // create a map of only the specified type of load
        std::multimap<boundary_id_type, MAST::BoundaryCondition*>::const_iterator
        it = _side_bc_map.begin(), end = _side_bc_map.end();
        for ( ; it!= end; it++)
            if (it->second->type() == MAST::SMALL_DISTURBANCE_MOTION)
                local_side_bc_map.insert(*it);
    }
    
    {
        // create a map of only the specified type of load
        std::multimap<subdomain_id_type, MAST::BoundaryCondition*>::const_iterator
        it = _vol_bc_map.begin(), end = _vol_bc_map.end();
        for ( ; it!= end; it++)
            if (it->second->type() == MAST::SMALL_DISTURBANCE_MOTION)
                local_vol_bc_map.insert(*it);
    }
    
    for ( ; el != end_el; ++el) {
        
        const Elem* elem = *el;
        
        dof_map.dof_indices (elem, dof_indices);
        
        const MAST::ElementPropertyCardBase& p_card = this->get_property_card(*elem);
        
        // create the structural element for analysis
        structural_elem.reset(MAST::build_structural_element
                              (_system, *elem, p_card).release());
        
        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        sol.resize(ndofs);
        vec.resize(ndofs);
        mat.resize(ndofs, ndofs);
        
        //for (unsigned int i=0; i<dof_indices.size(); i++)
        //    sol(i) = (*localized_solution)(dof_indices[i]);
        structural_elem->local_solution.resize(sol.size());
        //structural_elem->transform_to_local_system(sol, structural_elem->local_solution);
        
        // now get the vector values
        structural_elem->side_external_force(false, vec, mat,
                                             local_side_bc_map);
        structural_elem->volume_external_force(false, vec, mat,
                                               local_vol_bc_map);
        
        // constrain the vector
        _system.get_dof_map().constrain_element_vector(vec, dof_indices);
        
        // add to the global vector
        F.add_vector(vec, dof_indices);
    }
    
    F.close();
}





void
MAST::StructuralSystemAssembly::_assemble_matrices_for_modal_analysis(SparseMatrix<Number>&  matrix_A,
                                                                      SparseMatrix<Number>&  matrix_B,
                                                                      const MAST::FieldFunctionBase* param,
                                                                      const NumericVector<Number>* static_sol,
                                                                      const NumericVector<Number>* static_sol_sens) {
#ifndef LIBMESH_USE_COMPLEX_NUMBERS
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    DenseVector<Real> vec, sol;
    DenseMatrix<Real> mat1, mat2;
    std::vector<dof_id_type> dof_indices;
    const DofMap& dof_map = _system.get_dof_map();
    std::auto_ptr<MAST::StructuralElementBase> structural_elem;
    
    const bool if_exchange_AB_matrices =
    _system.get_equation_systems().parameters.get<bool>("if_exchange_AB_matrices");
    
    AutoPtr<NumericVector<Number> > localized_solution, localized_solution_sens;
    if (static_sol) { // static deformation is provided
        localized_solution.reset(NumericVector<Number>::build(_system.comm()).release());
        localized_solution->init(_system.n_dofs(), _system.n_local_dofs(),
                                 _system.get_dof_map().get_send_list(),
                                 false, GHOSTED);
        static_sol->localize(*localized_solution, _system.get_dof_map().get_send_list());
        // do the same for sensitivity if the parameters were provided
        if (param) {
            // make sure that the sensitivity was also provided
            libmesh_assert(static_sol_sens);
            
            localized_solution_sens.reset(NumericVector<Number>::build(_system.comm()).release());
            localized_solution_sens->init(_system.n_dofs(), _system.n_local_dofs(),
                                     _system.get_dof_map().get_send_list(),
                                     false, GHOSTED);
            static_sol_sens->localize(*localized_solution_sens, _system.get_dof_map().get_send_list());
        }
    }
    
    MeshBase::const_element_iterator       el     = _system.get_mesh().active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = _system.get_mesh().active_local_elements_end();
    
    for ( ; el != end_el; ++el) {
        
        const Elem* elem = *el;
        
        dof_map.dof_indices (elem, dof_indices);
        
        const MAST::ElementPropertyCardBase& p_card = this->get_property_card(*elem);
        
        // create the structural element for analysis
        structural_elem.reset(MAST::build_structural_element
                              (_system, *elem, p_card).release());
        
        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        sol.resize(ndofs);
        vec.resize(ndofs);
        mat1.resize(ndofs, ndofs);
        mat2.resize(ndofs, ndofs);

        // resize the solution to the correct size
        structural_elem->local_solution.resize(sol.size());
        structural_elem->local_velocity.resize(sol.size());
        structural_elem->local_acceleration.resize(sol.size());

        // if the static solution is provided, initialize the element solution
        if (static_sol) {
            for (unsigned int i=0; i<dof_indices.size(); i++)
                sol(i) = (*localized_solution)(dof_indices[i]);
            
            structural_elem->transform_to_local_system(sol, structural_elem->local_solution);
        }
        
        
        // now get the matrices
        if (!param) {
            structural_elem->internal_force(true, vec, mat1, true);
            structural_elem->prestress_force(true, vec, mat1);
            if (static_sol) {
                // get the Jacobian due to external loads, including thermal stresses
                structural_elem->side_external_force(true, vec, mat1,
                                                     _side_bc_map);
                structural_elem->volume_external_force(true, vec, mat1,
                                                       _vol_bc_map);
            }
            mat1.scale(-1.);
            structural_elem->inertial_force(true, vec, mat2);
        }
        else {
            // get the solution sensitivity if the static solution was provided
            structural_elem->local_solution_sens.resize(sol.size());
            structural_elem->local_velocity_sens.resize(sol.size());
            structural_elem->local_acceleration_sens.resize(sol.size());

            if (static_sol) {
                for (unsigned int i=0; i<dof_indices.size(); i++)
                    sol(i) = (*localized_solution_sens)(dof_indices[i]);
                
                structural_elem->transform_to_local_system(sol, structural_elem->local_solution_sens);
                structural_elem->side_external_force_sensitivity(true, vec, mat2,
                                                                 _side_bc_map);
                structural_elem->volume_external_force_sensitivity(true, vec, mat2,
                                                                   _vol_bc_map);
            }

            structural_elem->sensitivity_param = param;
            structural_elem->internal_force_sensitivity(true, vec, mat1, true);
            structural_elem->prestress_force_sensitivity(true, vec, mat1); mat1.scale(-1.);
            structural_elem->inertial_force_sensitivity(true, vec, mat2);
        }
        
        // constrain the element matrices.
        _system.get_dof_map().constrain_element_matrix(mat1, dof_indices);
        _system.get_dof_map().constrain_element_matrix(mat2, dof_indices);
        
        // add to the global matrices
        if (if_exchange_AB_matrices)
        {
            matrix_A.add_matrix (mat2, dof_indices); // mass
            matrix_B.add_matrix (mat1, dof_indices); // stiffness
        }
        else
        {
            matrix_A.add_matrix (mat1, dof_indices); // stiffness
            matrix_B.add_matrix (mat2, dof_indices); // mass
        }
    }
#endif // LIBMESH_USE_COMPLEX_NUMBERS
}



void
MAST::StructuralSystemAssembly::_assemble_matrices_for_buckling_analysis(SparseMatrix<Number>&  matrix_A,
                                                                         SparseMatrix<Number>&  matrix_B,
                                                                         const MAST::FieldFunctionBase* param,
                                                                         const NumericVector<Number>* static_sol,
                                                                         const NumericVector<Number>* static_sol_sens) {
#ifndef LIBMESH_USE_COMPLEX_NUMBERS
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    DenseVector<Real> vec, sol;
    DenseMatrix<Real> mat1, mat2, mat3;
    std::vector<dof_id_type> dof_indices;
    const DofMap& dof_map = _system.get_dof_map();
    std::auto_ptr<MAST::StructuralElementBase> structural_elem;
    
    const bool if_exchange_AB_matrices =
    _system.get_equation_systems().parameters.get<bool>("if_exchange_AB_matrices");
    
    AutoPtr<NumericVector<Number> > localized_solution, localized_solution_sens;
    if (static_sol) { // static deformation is provided
        localized_solution.reset(NumericVector<Number>::build(_system.comm()).release());
        localized_solution->init(_system.n_dofs(), _system.n_local_dofs(),
                                 _system.get_dof_map().get_send_list(),
                                 false, GHOSTED);
        static_sol->localize(*localized_solution, _system.get_dof_map().get_send_list());
        // do the same for sensitivity if the parameters were provided
        if (param) {
            // make sure that the sensitivity was also provided
            libmesh_assert(static_sol_sens);
            
            localized_solution_sens.reset(NumericVector<Number>::build(_system.comm()).release());
            localized_solution_sens->init(_system.n_dofs(), _system.n_local_dofs(),
                                          _system.get_dof_map().get_send_list(),
                                          false, GHOSTED);
            static_sol_sens->localize(*localized_solution_sens, _system.get_dof_map().get_send_list());
        }
    }
    
    MeshBase::const_element_iterator       el     = _system.get_mesh().active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = _system.get_mesh().active_local_elements_end();
    
    for ( ; el != end_el; ++el) {
        
        const Elem* elem = *el;
        
        dof_map.dof_indices (elem, dof_indices);
        
        const MAST::ElementPropertyCardBase& p_card = this->get_property_card(*elem);
        
        // create the structural element for analysis
        structural_elem.reset(MAST::build_structural_element
                              (_system, *elem, p_card).release());
        
        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        sol.resize(ndofs);
        vec.resize(ndofs);
        mat1.resize(ndofs, ndofs);
        mat2.resize(ndofs, ndofs);
        
        // resize the solution to the correct size
        structural_elem->local_solution.resize(sol.size());
        
        if (!param) {
            // set the local solution to zero for the load INdependent stiffness matrix
            structural_elem->local_solution.resize(sol.size());
            structural_elem->internal_force(true, vec, mat1, true); mat1.scale(-1.);

            // if the static solution is provided, initialize the element solution
            if (static_sol) {
                for (unsigned int i=0; i<dof_indices.size(); i++)
                    sol(i) = (*localized_solution)(dof_indices[i]);
                
                structural_elem->transform_to_local_system(sol, structural_elem->local_solution);
                
                // if displacement is zero, mat1 = mat2
                structural_elem->internal_force(true, vec, mat2, true);
                mat2.add(1., mat1); // subtract to get the purely load dependent part
                
                structural_elem->side_external_force(true, vec, mat2,
                                                     _side_bc_map);
                structural_elem->volume_external_force(true, vec, mat2,
                                                       _vol_bc_map);
            }
            
            structural_elem->prestress_force(true, vec, mat2);
        }
        else {
            structural_elem->sensitivity_param = param;
            // set the local solution to zero for the load INdependent stiffness matrix
            structural_elem->local_solution.resize(sol.size());
            structural_elem->local_solution_sens.resize(sol.size());
            structural_elem->internal_force_sensitivity(true, vec, mat1, true); mat1.scale(-1.);

            // now use the solution to get the load dependent stiffness matrix
            if (static_sol) {
                // displacement
                for (unsigned int i=0; i<dof_indices.size(); i++)
                    sol(i) = (*localized_solution)(dof_indices[i]);
                structural_elem->transform_to_local_system(sol, structural_elem->local_solution);
                
                // displacement sensitivity
                for (unsigned int i=0; i<dof_indices.size(); i++)
                    sol(i) = (*localized_solution_sens)(dof_indices[i]);
                structural_elem->transform_to_local_system(sol, structural_elem->local_solution_sens);
                
                // if displacement is zero, mat1 = mat2
                structural_elem->internal_force_sensitivity(true, vec, mat2, true);
                mat2.add(1., mat1); // subtract to get the purely load dependent part
                structural_elem->side_external_force_sensitivity(true, vec, mat2,
                                                                 _side_bc_map);
                structural_elem->volume_external_force_sensitivity(true, vec, mat2,
                                                                   _vol_bc_map);
            }
            
            structural_elem->prestress_force_sensitivity(true, vec, mat2);
        }
        
        // constrain the element matrices.
        _system.get_dof_map().constrain_element_matrix(mat1, dof_indices);
        _system.get_dof_map().constrain_element_matrix(mat2, dof_indices);
        
        // add to the global matrices
        if (if_exchange_AB_matrices)
        {
            matrix_A.add_matrix (mat2, dof_indices); // load dependent
            matrix_B.add_matrix (mat1, dof_indices); // load independent
        }
        else
        {
            matrix_A.add_matrix (mat1, dof_indices); // load independent
            matrix_B.add_matrix (mat2, dof_indices); // load dependent
        }
    }
    
#endif // LIBMESH_USE_COMPLEX_NUMBERS
}



void
MAST::StructuralSystemAssembly::calculate_max_elem_stress(const NumericVector<Number>& X,
                                                          std::vector<Real>& stress,
                                                          const MAST::FieldFunctionBase* param) {
#ifndef LIBMESH_USE_COMPLEX_NUMBERS

    // resize the stress vector to store values for each element
    stress.resize(_system.get_mesh().n_active_local_elem());
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    DenseVector<Real> sol;
    DenseMatrix<Real> mat;
    std::vector<dof_id_type> dof_indices;
    const DofMap& dof_map = _system.get_dof_map();
    std::auto_ptr<MAST::StructuralElementBase> structural_elem;
    
    AutoPtr<NumericVector<Number> > localized_solution =
    NumericVector<Number>::build(_system.comm());
    localized_solution->init(_system.n_dofs(), _system.n_local_dofs(),
                             _system.get_dof_map().get_send_list(),
                             false, GHOSTED);
    X.localize(*localized_solution, _system.get_dof_map().get_send_list());
    
    MeshBase::const_element_iterator       el     = _system.get_mesh().active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = _system.get_mesh().active_local_elements_end();
    
    unsigned int counter = 0;
    
    for ( ; el != end_el; ++el) {
        
        const Elem* elem = *el;
        
        dof_map.dof_indices (elem, dof_indices);
        
        const MAST::ElementPropertyCardBase& p_card = this->get_property_card(*elem);
        
        // create the structural element for analysis
        structural_elem.reset(MAST::build_structural_element
                              (_system, *elem, p_card).release());
        
        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        sol.resize(ndofs);
        mat.resize(ndofs, ndofs);
        
        for (unsigned int i=0; i<dof_indices.size(); i++)
            sol(i) = (*localized_solution)(dof_indices[i]);
        structural_elem->local_solution.resize(sol.size());
        structural_elem->transform_to_local_system(sol, structural_elem->local_solution);
        
        // now get the vector values
        if (!param) {
            stress[counter] = structural_elem->max_von_mises_stress();
        }
        else {
            structural_elem->sensitivity_param = param;
            stress[counter] = structural_elem->max_von_mises_stress_sensitivity();
        }
        counter++;
    }
    
#endif // LIBMESH_USE_COMPLEX_NUMBERS
}



void
MAST::StructuralSystemAssembly::get_dirichlet_dofs(std::set<unsigned int>& dof_ids) const {
    
    dof_ids.clear();
    
    // first prepare a map of boundary ids and the constrained vars on that
    // boundary
    std::map<boundary_id_type, std::vector<unsigned int> >  constrained_vars_map;
    
    // now populate the map for all the give boundaries
    std::multimap<boundary_id_type, MAST::BoundaryCondition*>::const_iterator
    it = _side_bc_map.begin(), end = _side_bc_map.end();
    
    for ( ; it != end; it++)
        if (it->second->type() == MAST::DISPLACEMENT_DIRICHLET) {
            // get the displacement dirichlet condition
            DirichletBoundary& dirichlet_b =
            (dynamic_cast<MAST::DisplacementDirichletBoundaryCondition*>(it->second))->dirichlet_boundary();
            
            constrained_vars_map[it->first] = dirichlet_b.variables;
        }
    
    
    // now collect the ids that correspond to the specified boundary conditions
    
    // Get a constant reference to the mesh object
    const MeshBase& mesh = _system.get_mesh();
    
    // The dimension that we are running.
    const unsigned int dim = mesh.mesh_dimension();
    
    // Get a constant reference to the Finite Element type
    // for the first (and only) variable in the system.
    FEType fe_type = _system.get_dof_map().variable_type(0);
    
    const DofMap& dof_map = _system.get_dof_map();
    
    std::vector<dof_id_type> dof_indices;
    
    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    
    for ( ; el != end_el; ++el) {
        const Elem* elem = *el;
        
        // boundary condition is applied only on sides with no neighbors
        // and if the side's boundary id has a boundary condition tag on it
        for (unsigned int s=0; s<elem->n_sides(); s++)
            if ((*el)->neighbor(s) == NULL &&
                mesh.boundary_info->n_boundary_ids(elem, s)) {
                
                std::vector<boundary_id_type> bc_ids = mesh.boundary_info->boundary_ids(elem, s);
                
                for (unsigned int i_bid=0; i_bid<bc_ids.size(); i_bid++)
                    if (constrained_vars_map.count(bc_ids[i_bid])) {
                        
                        const std::vector<unsigned int>& vars = constrained_vars_map[bc_ids[i_bid]];
                        // now iterate over each constrained variable for this boundary
                        // and collect its dofs
                        for (unsigned int i_var=0; i_var<vars.size(); i_var++) {
                            
                            dof_indices.clear();
                            dof_map.dof_indices (*el, dof_indices, vars[i_var]);
                            
                            // All boundary dofs are Dirichlet dofs in this case
                            std::vector<unsigned int> side_dofs;
                            FEInterface::dofs_on_side(*el, dim, fe_type,
                                                      s, side_dofs);
                            
                            for(unsigned int ii=0; ii<side_dofs.size(); ii++)
                                dof_ids.insert(dof_indices[side_dofs[ii]]);
                        }
                    } // end of boundary loop
            } // end of side loop
    }// end of element loop
    return;
    
}

