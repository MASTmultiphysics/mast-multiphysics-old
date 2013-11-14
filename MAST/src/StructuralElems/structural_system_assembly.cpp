
// libMesh include files
#include "libmesh/nonlinear_solver.h"
#include "libmesh/condensed_eigen_system.h"
#include "libmesh/fe_interface.h"
#include "libmesh/parameter_vector.h"

// MAST includes
#include "StructuralElems/structural_system_assembly.h"
#include "StructuralElems/surface_pressure_load.h"
#include "StructuralElems/structural_elem_base.h"
#include "PropertyCards/element_property_card_base.h"


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
MAST::StructuralSystemAssembly::residual_and_jacobian (const NumericVector<Number>& X,
                                                       NumericVector<Number>* R,
                                                       SparseMatrix<Number>*  J,
                                                       NonlinearImplicitSystem& S) {
    
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
}




bool
MAST::StructuralSystemAssembly::sensitivity_assemble (const ParameterVector& params,
                                                      const unsigned int i,
                                                      NumericVector<Number>& sensitivity_rhs) {
    SensitivityParameters sens_params;
    sens_params.add_parameter(this->get_parameter(params[i]), 1);
    
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
    
    // currently, all relevant parameter sensitivities are calculated
    return true;
}



void
MAST::StructuralSystemAssembly::assemble() {
    
    SparseMatrix<Number>&  matrix_A = *(dynamic_cast<EigenSystem&>(_system).matrix_A);
    SparseMatrix<Number>&  matrix_B = *(dynamic_cast<EigenSystem&>(_system).matrix_B);

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
}




bool
MAST::StructuralSystemAssembly::sensitivity_assemble (const ParameterVector& params,
                                                      const unsigned int i,
                                                      SparseMatrix<Number>* sensitivity_A,
                                                      SparseMatrix<Number>* sensitivity_B) {
    SensitivityParameters sens_params;
    sens_params.add_parameter(this->get_parameter(params[i]), 1);
    
    
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
    // currently, all relevant parameter sensitivities are calculated
    return true;
}




void
MAST::StructuralSystemAssembly::_assemble_residual_and_jacobian (const NumericVector<Number>& X,
                                                                 NumericVector<Number>* R,
                                                                 SparseMatrix<Number>*  J,
                                                                 NonlinearImplicitSystem& S,
                                                                 const MAST::SensitivityParameters* params) {
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    DenseVector<Real> vec, sol;
    DenseMatrix<Real> mat;
    std::vector<dof_id_type> dof_indices;
    const DofMap& dof_map = _system.get_dof_map();
    std::auto_ptr<MAST::StructuralElementBase> structural_elem;
    
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
            sol(i) = X(dof_indices[i]);
        
        structural_elem->local_solution = sol;
        
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
}



void
MAST::StructuralSystemAssembly::_assemble_matrices_for_modal_analysis(const NumericVector<Number>& X,
                                                                      SparseMatrix<Number>&  matrix_A,
                                                                      SparseMatrix<Number>&  matrix_B,
                                                                      const MAST::SensitivityParameters* params) {
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    DenseVector<Real> vec, sol;
    DenseMatrix<Real> mat1, mat2;
    std::vector<dof_id_type> dof_indices;
    const DofMap& dof_map = _system.get_dof_map();
    std::auto_ptr<MAST::StructuralElementBase> structural_elem;

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
        
        for (unsigned int i=0; i<dof_indices.size(); i++)
            sol(i) = X(dof_indices[i]);
        
        structural_elem->local_solution = sol;
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
}



void
MAST::StructuralSystemAssembly::_assemble_matrices_for_buckling_analysis(const NumericVector<Number>& X,
                                                                         SparseMatrix<Number>&  matrix_A,
                                                                         SparseMatrix<Number>&  matrix_B,
                                                                         const MAST::SensitivityParameters* params) {
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    DenseVector<Real> vec, sol;
    DenseMatrix<Real> mat1, mat2, mat3;
    std::vector<dof_id_type> dof_indices;
    const DofMap& dof_map = _system.get_dof_map();
    std::auto_ptr<MAST::StructuralElementBase> structural_elem;
    
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
            sol(i) = X(dof_indices[i]);
        
        if (!params) {
            // set the local solution to zero for the load INdependent stiffness matrix
            structural_elem->local_solution.resize(sol.size());
            structural_elem->internal_force(true, vec, mat1); mat1.scale(-1.);
            
            // now use the solution to get the load dependent stiffness matrix
            structural_elem->local_solution = sol;
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
            structural_elem->local_solution = sol;
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
        dof_map.dof_indices (*el, dof_indices, 2); // uz
        
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
        dof_map.dof_indices(*el, dof_indices, 1); // uy
        for (unsigned int i=0; i<dof_indices.size(); i++)
            dof_ids.insert(dof_indices[i]);
        
        dof_indices.clear();
        dof_map.dof_indices(*el, dof_indices, 5); // tz
        for (unsigned int i=0; i<dof_indices.size(); i++)
            dof_ids.insert(dof_indices[i]);
    } // end of element loop
    
    /**
     * All done!
     */
    return;
    
}



void assemble_force_vec(System& sys,
                             SurfacePressureLoad& surf_press,
                             SurfaceMotionBase& surf_motion,
                             NumericVector<Number>& fvec)
{
//#ifdef LIBMESH_USE_COMPLEX_NUMBERS
//    // Get a constant reference to the mesh object.
//    const MeshBase& mesh = sys.get_mesh();
//    
//    // The dimension that we are running.
//    const unsigned int dim = mesh.mesh_dimension();
//    
//    // A reference to the \p DofMap object for this system.  The \p DofMap
//    // object handles the index translation from node and element numbers
//    // to degree of freedom numbers.
//    const DofMap& dof_map = sys.get_dof_map();
//    std::vector<dof_id_type> dof_indices;
//
//    // The element mass and stiffness matrices.
//    DenseVector<Number>   fvec_e;
//    
//    FEType fe_type(FIRST, LAGRANGE);
//    AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));
//    QGauss qrule (dim, FIFTH);
//    fe->attach_quadrature_rule (&qrule);
//    const std::vector<Real>& JxW = fe->get_JxW();
//    const std::vector<Point>& q_point = fe->get_xyz();
//    const std::vector<std::vector<Real> >& phi = fe->get_phi();
//    
//    
//    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
//    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
//    Number press, dpress;
//    DenseVector<Number> utrans, dn_rot; utrans.resize(3); dn_rot.resize(3);
//    Point normal; normal.zero(); normal(dim) = -1.;
//    
//    fvec.zero(); fvec.close();
//    
//    for ( ; el != end_el; ++el)
//    {
//        dof_map.dof_indices (*el, dof_indices);
//        fvec_e.resize (dof_indices.size());
//        fe->reinit (*el);
//        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
//        {
//            surf_press.surface_pressure(q_point[qp], press, dpress);
//            surf_motion.surface_velocity_frequency_domain(q_point[qp], normal,
//                                                          utrans, dn_rot);
////            press = 0.;
////            dpress = Complex(2./4.*std::real(dn_rot(0)),  2./4./.1*std::imag(utrans(1)));
////            std::cout << q_point[qp](0)
////            << "  " << std::real(utrans(1))
////            << "  " << std::imag(utrans(1))
////            << "  " << std::real(dn_rot(0))
////            << "  " << std::imag(dn_rot(0))
////            << "  " << std::real(press)
////            << "  " << std::imag(press)
////            << "  " << std::real(dpress)
////            << "  " << std::imag(dpress) << std::endl;
//            
//            for (unsigned int i=0; i<3; i++)
//                for (unsigned int iphi=0; iphi<phi.size(); iphi++)
//                    fvec_e(i*phi.size()+iphi) += JxW[qp] * phi[iphi][qp] *
//                    ( press * dn_rot(i) + // steady pressure
//                     dpress * normal(i)); // unsteady pressure
//        }
//        
//        fvec.add_vector (fvec_e, dof_indices);
//    }
//    fvec.close();
//#endif // LIBMESH_USE_COMPLEX_NUMBERS
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







