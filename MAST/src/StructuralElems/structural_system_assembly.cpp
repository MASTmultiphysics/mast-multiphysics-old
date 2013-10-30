
// libMesh include files
#include "libmesh/nonlinear_solver.h"
#include "libmesh/condensed_eigen_system.h"
#include "libmesh/fe_interface.h"

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
MAST::StructuralSystemAssembly::assemble() {
    
    switch (_analysis_type) {
        case MAST::MODAL:
            _assemble_matrices_for_modal_analysis();
            break;
            
        default:
            libmesh_error();
            break;
    }
}




void
MAST::StructuralSystemAssembly::residual_and_jacobian (const NumericVector<Number>& X,
                                                       NumericVector<Number>* R,
                                                       SparseMatrix<Number>*  J,
                                                       NonlinearImplicitSystem& S) {
    
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
        sol.resize(dof_indices.size());
        vec.resize(dof_indices.size());
        mat.resize(dof_indices.size(), dof_indices.size());
        
        for (unsigned int i=0; i<dof_indices.size(); i++)
            sol(i) = X(dof_indices[i]);
        
        structural_elem->local_solution = sol;
        
        // now get the vector values
        structural_elem->internal_force(J!=NULL?true:false, vec, mat);
        
        // add to the global matrices
        if (R) R->add_vector(vec, dof_indices);
        if (J) J->add_matrix(mat, dof_indices);
    }
}



void
MAST::StructuralSystemAssembly::_assemble_matrices_for_modal_analysis() {
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    DenseVector<Real> vec, sol;
    DenseMatrix<Real> mat1, mat2;
    std::vector<dof_id_type> dof_indices;
    const DofMap& dof_map = _system.get_dof_map();
    std::auto_ptr<MAST::StructuralElementBase> structural_elem;

    const bool if_exchange_AB_matrices =
    _system.get_equation_systems().parameters.get<bool>("if_exchange_AB_matrices");
    
    SparseMatrix<Number>&  matrix_A = *(dynamic_cast<EigenSystem&>(_system).matrix_A);
    SparseMatrix<Number>&  matrix_B = *(dynamic_cast<EigenSystem&>(_system).matrix_B);

    
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
        sol.resize(dof_indices.size());
        vec.resize(dof_indices.size());
        mat1.resize(dof_indices.size(), dof_indices.size());
        mat2.resize(dof_indices.size(), dof_indices.size());
        
        for (unsigned int i=0; i<dof_indices.size(); i++)
            sol(i) = (*_system.solution)(dof_indices[i]);
        
        structural_elem->local_solution = sol;
        sol.zero();
        structural_elem->local_acceleration = sol;
        
        
        // now get the vector values
        structural_elem->internal_force(true, vec, mat1); mat1.scale(-1.);
        structural_elem->inertial_force(true, vec, mat2);
        
        // add to the global matrices
        if (if_exchange_AB_matrices)
        {
            matrix_A.add_matrix (mat1, dof_indices);
            matrix_B.add_matrix (mat2, dof_indices);
        }
        else
        {
            matrix_A.add_matrix (mat2, dof_indices);
            matrix_B.add_matrix (mat1, dof_indices);
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






