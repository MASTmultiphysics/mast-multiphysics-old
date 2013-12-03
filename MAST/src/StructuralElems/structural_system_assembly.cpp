
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


MAST::StructuralSystemAssembly::StructuralSystemAssembly(System& sys,
                                                         MAST::StructuralAnalysisType t,
                                                         GetPot& infile):

_system(sys),
_analysis_type(t),
_infile(infile),
_if_same_property_for_all_elems(false),
_property(NULL)
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
}


void
MAST::StructuralSystemAssembly::add_side_load(boundary_id_type bid,
                                              MAST::BoundaryCondition& load) {
    // make sure that this boundary and load haven't already been applied
    std::pair<std::multimap<boundary_id_type, MAST::BoundaryCondition*>::const_iterator,
    std::multimap<boundary_id_type, MAST::BoundaryCondition*>::const_iterator> it =
    _side_bc_map.equal_range(bid);
    
    for ( ; it.first != it.second; it.first++)
        libmesh_assert(it.first->second != &load);
    
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
MAST::StructuralSystemAssembly::set_property_for_all_elems(const MAST::ElementPropertyCardBase& prop) {
    _if_same_property_for_all_elems = true;
    _element_property.clear();
    _property = &prop;
}


void
MAST::StructuralSystemAssembly::set_property_for_elem(const Elem& e,
                                                      const MAST::ElementPropertyCardBase& prop) {
    _if_same_property_for_all_elems = false;
    _element_property[&e] = &prop;
}



const MAST::ElementPropertyCardBase&
MAST::StructuralSystemAssembly::get_property_card(const Elem& elem) const {
    if (_if_same_property_for_all_elems)
        return *_property;
    else {
        std::map<const Elem*, const MAST::ElementPropertyCardBase*>::const_iterator
        elem_p_it = _element_property.find(&elem);
        libmesh_assert(elem_p_it != _element_property.end());
        
        return *elem_p_it->second;
    }
}

void
MAST::StructuralSystemAssembly::add_parameter(Real* par, MAST::FunctionBase* f) {
    // make sure valid values are given
    libmesh_assert(par);
    libmesh_assert(f);
    // make sure that the function is dependent on the parameters
    libmesh_assert(f->depends_on(par));
    // make sure it does not already exist in the map
    libmesh_assert(!_parameter_map.count(par));
    
    // now add this to the map
    bool insert_success = _parameter_map.insert
    (std::map<Real*, MAST::FunctionBase*>::value_type(par, f)).second;
    
    libmesh_assert(insert_success);
}



const MAST::FunctionBase*
MAST::StructuralSystemAssembly::get_parameter(Real* par) const {
    // make sure valid values are given
    libmesh_assert(par);
    
    std::map<Real*, const MAST::FunctionBase*>::const_iterator
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
    
    SensitivityParameters sens_params;
    sens_params.add_parameter(this->get_parameter(params[i]), 1);

    sensitivity_rhs.zero();

    switch (_analysis_type) {
        case MAST::STATIC:
        case MAST::DYNAMIC:
            libmesh_assert(_system.system_type() == "NonlinearImplicit");
            _assemble_residual_and_jacobian(*_system.solution, &sensitivity_rhs,
                                            NULL,
                                            dynamic_cast<NonlinearImplicitSystem&>(_system),
                                            &sens_params);
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
    
    switch (_analysis_type) {
        case MAST::MODAL:
            _assemble_matrices_for_modal_analysis(*_system.solution,
                                                  matrix_A,
                                                  matrix_B,
                                                  NULL);
            break;

        case MAST::BUCKLING:
            _assemble_matrices_for_buckling_analysis(*_system.solution,
                                                     matrix_A,
                                                     matrix_B,
                                                     NULL);
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
    
    SensitivityParameters sens_params;
    sens_params.add_parameter(this->get_parameter(params[i]), 1);
    
    sensitivity_A->zero();
    sensitivity_B->zero();
    
    switch (_analysis_type) {
        case MAST::MODAL:
            _assemble_matrices_for_modal_analysis(*_system.solution,
                                                  *sensitivity_A,
                                                  *sensitivity_B,
                                                  &sens_params);
            break;
            
        case MAST::BUCKLING:
            _assemble_matrices_for_buckling_analysis(*_system.solution,
                                                     *sensitivity_A,
                                                     *sensitivity_B,
                                                     &sens_params);
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
                                                                 const MAST::SensitivityParameters* params) {
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
        if (!params) {
            structural_elem->internal_force(J!=NULL?true:false,
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
            structural_elem->sensitivity_params = params;
            structural_elem->internal_force_sensitivity(J!=NULL?true:false,
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
MAST::StructuralSystemAssembly::_assemble_matrices_for_modal_analysis(const NumericVector<Number>& X,
                                                                      SparseMatrix<Number>&  matrix_A,
                                                                      SparseMatrix<Number>&  matrix_B,
                                                                      const MAST::SensitivityParameters* params) {
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
        mat1.resize(ndofs, ndofs);
        mat2.resize(ndofs, ndofs);
        
        for (unsigned int i=0; i<dof_indices.size(); i++)
            sol(i) = (*localized_solution)(dof_indices[i]);
        
        structural_elem->local_solution.resize(sol.size());
        structural_elem->transform_to_local_system(sol, structural_elem->local_solution);
        sol.zero();
        structural_elem->local_acceleration = sol;
        
        
        // now get the matrices
        if (!params) {
            structural_elem->internal_force(true, vec, mat1); mat1.scale(-1.);
            structural_elem->inertial_force(true, vec, mat2);
        }
        else {
            structural_elem->sensitivity_params = params;
            structural_elem->internal_force_sensitivity(true, vec, mat1); mat1.scale(-1.);
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
MAST::StructuralSystemAssembly::_assemble_matrices_for_buckling_analysis(const NumericVector<Number>& X,
                                                                         SparseMatrix<Number>&  matrix_A,
                                                                         SparseMatrix<Number>&  matrix_B,
                                                                         const MAST::SensitivityParameters* params) {
#ifndef LIBMESH_USE_COMPLEX_NUMBERS
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    DenseVector<Real> vec, sol;
    DenseMatrix<Real> mat1, mat2, mat3;
    std::vector<dof_id_type> dof_indices;
    const DofMap& dof_map = _system.get_dof_map();
    std::auto_ptr<MAST::StructuralElementBase> structural_elem;
    
    AutoPtr<NumericVector<Number> > localized_solution =
    NumericVector<Number>::build(_system.comm());
    localized_solution->init(_system.n_dofs(), _system.n_local_dofs(),
                             _system.get_dof_map().get_send_list(),
                             false, GHOSTED);
    X.localize(*localized_solution, _system.get_dof_map().get_send_list());

    const bool if_exchange_AB_matrices =
    _system.get_equation_systems().parameters.get<bool>("if_exchange_AB_matrices");
    
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
        mat3.resize(ndofs, ndofs);
        
        for (unsigned int i=0; i<dof_indices.size(); i++)
            sol(i) = (*localized_solution)(dof_indices[i]);
        
        if (!params) {
            // set the local solution to zero for the load INdependent stiffness matrix
            structural_elem->local_solution.resize(sol.size());
            structural_elem->internal_force(true, vec, mat1); mat1.scale(-1.);
            
            // now use the solution to get the load dependent stiffness matrix
            structural_elem->transform_to_local_system(sol, structural_elem->local_solution);
            if (sol.l2_norm() > 0.) { // if displacement is zero, mat1 = mat2
                structural_elem->internal_force(true, vec, mat2);
                mat2.add(1., mat1); // subtract to get the purely load dependent part
            }
            structural_elem->prestress_force(true, vec, mat3);
            mat2.add(1., mat3);
        }
        else {
            structural_elem->sensitivity_params = params;
            // set the local solution to zero for the load INdependent stiffness matrix
            structural_elem->local_solution.resize(sol.size());
            structural_elem->internal_force_sensitivity(true, vec, mat1); mat1.scale(-1.);
            
            // now use the solution to get the load dependent stiffness matrix
            structural_elem->transform_to_local_system(sol, structural_elem->local_solution);
            if (sol.l2_norm() > 0.) { // if displacement is zero, mat1 = mat2
                structural_elem->internal_force_sensitivity(true, vec, mat2);
                mat2.add(1., mat1); // subtract to get the purely load dependent part
            }
            structural_elem->prestress_force_sensitivity(true, vec, mat3);
            mat2.add(1., mat3);
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
MAST::StructuralSystemAssembly::get_dirichlet_dofs(std::set<unsigned int>& dof_ids) const {
    dof_ids.clear();
    
    // Get a constant reference to the mesh object.
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
    
    for ( ; el != end_el; ++el)
    {
        dof_map.dof_indices (*el, dof_indices, 1); // uy
        
        // All boundary dofs are Dirichlet dofs in this case
        for (unsigned int s=0; s<(*el)->n_sides(); s++)
            if ((*el)->neighbor(s) == NULL)
            {
                std::vector<unsigned int> side_dofs;
                FEInterface::dofs_on_side(*el, dim, fe_type,
                                          s, side_dofs);
                
                for(unsigned int ii=0; ii<side_dofs.size(); ii++)
                    dof_ids.insert(dof_indices[side_dofs[ii]]);
            }
        
        // also add the dofs for variable u, v and tz
        dof_indices.clear();
        dof_map.dof_indices(*el, dof_indices, 0); // ux
        for (unsigned int i=0; i<dof_indices.size(); i++)
            dof_ids.insert(dof_indices[i]);
        
        dof_indices.clear();
        dof_map.dof_indices(*el, dof_indices, 2); // uz
        for (unsigned int i=0; i<dof_indices.size(); i++)
            dof_ids.insert(dof_indices[i]);
        
        dof_indices.clear();
        dof_map.dof_indices(*el, dof_indices, 3); // tx
        for (unsigned int i=0; i<dof_indices.size(); i++)
            dof_ids.insert(dof_indices[i]);

        dof_indices.clear();
        dof_map.dof_indices(*el, dof_indices, 4); // ty
        for (unsigned int i=0; i<dof_indices.size(); i++)
            dof_ids.insert(dof_indices[i]);
    } // end of element loop
    
    /**
     * All done!
     */
    return;
    
}

