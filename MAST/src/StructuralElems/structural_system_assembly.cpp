
// libMesh include files
#include "libmesh/nonlinear_solver.h"

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
            sol(i) = (*_system.solution)(dof_indices[i]);
        
        structural_elem->local_solution = sol;
        
        // now get the vector values
        structural_elem->internal_force(J!=NULL?true:false, vec, mat);
        
        // add to the global matrices
        if (R) R->add_vector(vec, dof_indices);
        if (J) J->add_matrix(mat, dof_indices);
    }
}




void assemble_plate_matrices(EquationSystems& es,
                             const std::string& system_name)
{
//    // Get a constant reference to the mesh object.
//    const MeshBase& mesh = es.get_mesh();
//    
//    // The dimension that we are running.
//    const unsigned int dim = mesh.mesh_dimension();
//    
//    // Get a reference to our system.
//    EigenSystem & eigen_system = es.get_system<EigenSystem> (system_name);
//    
//    // A reference to the two system matrices
//    SparseMatrix<Number>&  matrix_A = *eigen_system.matrix_A;
//    SparseMatrix<Number>&  matrix_B = *eigen_system.matrix_B;
//    
//    // A reference to the \p DofMap object for this system.  The \p DofMap
//    // object handles the index translation from node and element numbers
//    // to degree of freedom numbers.
//    const DofMap& dof_map = eigen_system.get_dof_map();
//    
//    // The element mass and stiffness matrices.
//    DenseMatrix<Number>   Me;
//    DenseMatrix<Number>   Ke;
//    
//    // This vector will hold the degree of freedom indices for
//    // the element.  These define where in the global system
//    // the element degrees of freedom get mapped.
//    std::vector<dof_id_type> dof_indices;
//    
//    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
//    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
//    
//    for ( ; el != end_el; ++el)
//    {
//        dof_map.dof_indices (*el, dof_indices);
//        Ke.resize (dof_indices.size(), dof_indices.size());
//        Me.resize (dof_indices.size(), dof_indices.size());
//        
//        FESystem::Quadrature::TrapezoidQuadrature q_rule_shear, q_rule_bending;
//        FESystem::FiniteElement::FELagrange fe, fe_tri6;
//        FESystem::Structures::DKTPlate dkt_plate;
//        q_rule_bending.init(2, 9);
//        
//        // initialize the geometric element
//        std::auto_ptr<FESystem::Mesh::ElemBase> elem(new FESystem::Mesh::Tri3(false));
//        
//        FESystem::Numerics::DenseMatrix<Real> basis; basis.resize(3, 3); basis.setToIdentity();
//        FESystem::Geometry::Point origin(3);
//        FESystem::Geometry::RectangularCoordinateSystem cs(origin, basis);
//        std::vector<FESystem::Mesh::Node*> nodes(3);
//        
//        for (unsigned int i=0; i<3; i++)
//        {
//            nodes[i] = new FESystem::Mesh::Node(cs);
//            for (unsigned int j=0; j<3; j++)
//                nodes[i]->setVal(j, (*el)->point(i)(j));
//            elem->setNode(i, *nodes[i]);
//        }
//        fe.reinit(*elem);
//        
//        dkt_plate.initialize(*elem, fe, fe_tri6, q_rule_bending, q_rule_bending,
//                             72.0e9, 0.33, 2700., 0.01);
//        
//        FESystem::Numerics::DenseMatrix<Real> plate_elem_mat, elem_mat;
//        FESystem::Numerics::LocalVector<Real> plate_elem_vec, elem_vec;
//        plate_elem_mat.resize(dkt_plate.getNElemDofs(), dkt_plate.getNElemDofs());
//        elem_mat.resize(18, 18);
//        plate_elem_vec.resize(dkt_plate.getNElemDofs());
//        elem_vec.resize(18);
//        
//        dkt_plate.calculateStiffnessMatrix(plate_elem_mat);
//        dkt_plate.transformMatrixToGlobalSystem(plate_elem_mat, elem_mat);
//        for (unsigned int i=0; i<18; i++)
//            for (unsigned int j=0; j<18; j++)
//                Ke(i, j) = elem_mat.getVal(i, j);
//        
//        dkt_plate.calculateDiagonalMassMatrix(plate_elem_vec);
//        dkt_plate.transformVectorToGlobalSystem(plate_elem_vec, elem_vec);
//        for (unsigned int j=0; j<18; j++)
//            Me(j, j) = elem_vec.getVal(j);
//        
//        // clear the pointers
//        elem.reset();
//        for (unsigned int i=0; i<3; i++)
//            delete nodes[i];
//        
//        if (es.parameters.get<bool>("if_exchange_AB_matrices"))
//        {
//            matrix_A.add_matrix (Me, dof_indices);
//            matrix_B.add_matrix (Ke, dof_indices);
//        }
//        else
//        {
//            matrix_A.add_matrix (Ke, dof_indices);
//            matrix_B.add_matrix (Me, dof_indices);
//        }
//    }
}



void assemble_beam_matrices(EquationSystems& es,
                            const std::string& system_name)
{
//    // Get a constant reference to the mesh object.
//    const MeshBase& mesh = es.get_mesh();
//    
//    // The dimension that we are running.
//    const unsigned int dim = mesh.mesh_dimension();
//    
//    // Get a reference to our system.
//    EigenSystem & eigen_system = es.get_system<EigenSystem> (system_name);
//    
//    // A reference to the two system matrices
//    SparseMatrix<Number>&  matrix_A = *eigen_system.matrix_A;
//    SparseMatrix<Number>&  matrix_B = *eigen_system.matrix_B;
//    
//    // A reference to the \p DofMap object for this system.  The \p DofMap
//    // object handles the index translation from node and element numbers
//    // to degree of freedom numbers.
//    const DofMap& dof_map = eigen_system.get_dof_map();
//    
//    // The element mass and stiffness matrices.
//    DenseMatrix<Number>   Me;
//    DenseMatrix<Number>   Ke;
//    
//    // This vector will hold the degree of freedom indices for
//    // the element.  These define where in the global system
//    // the element degrees of freedom get mapped.
//    std::vector<dof_id_type> dof_indices;
//    
//    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
//    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
//    
//    for ( ; el != end_el; ++el)
//    {
//        dof_map.dof_indices (*el, dof_indices);
//        Ke.resize (dof_indices.size(), dof_indices.size());
//        Me.resize (dof_indices.size(), dof_indices.size());
//        
//        FESystem::Quadrature::TrapezoidQuadrature q_rule_shear, q_rule_bending;
//        FESystem::FiniteElement::FELagrange fe, fe_tri6;
//        FESystem::Structures::EulerBernoulliBeam beam;
//        q_rule_bending.init(1, 9);
//        
//        // initialize the geometric element
//        std::auto_ptr<FESystem::Mesh::Edge2> elem(new FESystem::Mesh::Edge2(false));
//        
//        FESystem::Numerics::DenseMatrix<Real> basis; basis.resize(3, 3); basis.setToIdentity();
//        FESystem::Geometry::Point origin(3);
//        FESystem::Geometry::RectangularCoordinateSystem cs(origin, basis);
//        FESystem::Numerics::LocalVector<Real> yvec; yvec.resize(3); yvec.setVal(1, 1.);
//        std::vector<FESystem::Mesh::Node*> nodes(2);
//        
//        elem->setVectorForXYPlane(yvec);
//        for (unsigned int i=0; i<2; i++)
//        {
//            nodes[i] = new FESystem::Mesh::Node(cs);
//            for (unsigned int j=0; j<2; j++)
//                nodes[i]->setVal(j, (*el)->point(i)(j));
//            elem->setNode(i, *nodes[i]);
//        }
//        fe.reinit(*elem);
//        
//        beam.initialize(*elem, fe, q_rule_bending,
//                        72.0e9, 0.33, 2700., 8.33333E-08, 8.33333E-08, 0.01);
//        
//        FESystem::Numerics::DenseMatrix<Real> beam_elem_mat, elem_mat;
//        FESystem::Numerics::LocalVector<Real> beam_elem_vec, elem_vec;
//        beam_elem_mat.resize(beam.getNElemDofs(), beam.getNElemDofs());
//        elem_mat.resize(12, 12);
//        beam_elem_vec.resize(beam.getNElemDofs());
//        elem_vec.resize(12);
//        
//        beam.calculateStiffnessMatrix(beam_elem_mat);
//        beam.transformMatrixToGlobalSystem(beam_elem_mat, elem_mat);
//        for (unsigned int i=0; i<12; i++)
//            for (unsigned int j=0; j<12; j++)
//                Ke(i, j) = elem_mat.getVal(i, j);
//        
//        beam.calculateDiagonalMassMatrix(beam_elem_vec);
//        beam.transformVectorToGlobalSystem(beam_elem_vec, elem_vec);
//        for (unsigned int j=0; j<12; j++)
//            Me(j, j) = elem_vec.getVal(j);
//        
//        // clear the pointers
//        elem.reset();
//        for (unsigned int i=0; i<2; i++)
//            delete nodes[i];
//        
//        if (es.parameters.get<bool>("if_exchange_AB_matrices"))
//        {
//            matrix_A.add_matrix (Me, dof_indices);
//            matrix_B.add_matrix (Ke, dof_indices);
//        }
//        else
//        {
//            matrix_A.add_matrix (Ke, dof_indices);
//            matrix_B.add_matrix (Me, dof_indices);
//        }
//    }
}




void get_plate_dirichlet_dofs(EquationSystems& es,
                              const std::string& system_name,
                              std::set<unsigned int>& dirichlet_dof_ids)
{
//    dirichlet_dof_ids.clear();
//    
//    // Get a constant reference to the mesh object.
//    const MeshBase& mesh = es.get_mesh();
//    
//    // The dimension that we are running.
//    const unsigned int dim = mesh.mesh_dimension();
//    
//    // Get a reference to our system.
//    EigenSystem & eigen_system = es.get_system<EigenSystem> (system_name);
//    
//    // Get a constant reference to the Finite Element type
//    // for the first (and only) variable in the system.
//    FEType fe_type = eigen_system.get_dof_map().variable_type(0);
//    
//    const DofMap& dof_map = eigen_system.get_dof_map();
//    
//    std::vector<dof_id_type> dof_indices;
//    
//    
//    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
//    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
//    
//    for ( ; el != end_el; ++el)
//    {
//        dof_map.dof_indices (*el, dof_indices, 2); // uz
//
//        // All boundary dofs are Dirichlet dofs in this case
//        for (unsigned int s=0; s<(*el)->n_sides(); s++)
//            if ((*el)->neighbor(s) == NULL)
//            {
//                std::vector<unsigned int> side_dofs;
//                FEInterface::dofs_on_side(*el, dim, fe_type,
//                                          s, side_dofs);
//                
//                for(unsigned int ii=0; ii<side_dofs.size(); ii++)
//                    dirichlet_dof_ids.insert(dof_indices[side_dofs[ii]]);
//            }
//        
//        // also add the dofs for variable u, v and tz
//        dof_indices.clear();
//        dof_map.dof_indices(*el, dof_indices, 0); // ux
//        for (unsigned int i=0; i<dof_indices.size(); i++)
//            dirichlet_dof_ids.insert(dof_indices[i]);
//        
//        dof_indices.clear();
//        dof_map.dof_indices(*el, dof_indices, 1); // uy
//        for (unsigned int i=0; i<dof_indices.size(); i++)
//            dirichlet_dof_ids.insert(dof_indices[i]);
//        
//        dof_indices.clear();
//        dof_map.dof_indices(*el, dof_indices, 5); // tz
//        for (unsigned int i=0; i<dof_indices.size(); i++)
//            dirichlet_dof_ids.insert(dof_indices[i]);
//    } // end of element loop
//
//    /**
//     * All done!
//     */
//    return;
//    
}



void get_beam_dirichlet_dofs(EquationSystems& es,
                             const std::string& system_name,
                             std::set<unsigned int>& dirichlet_dof_ids)
{
//    dirichlet_dof_ids.clear();
//    
//    // Get a constant reference to the mesh object.
//    const MeshBase& mesh = es.get_mesh();
//    
//    // The dimension that we are running.
//    const unsigned int dim = mesh.mesh_dimension();
//    
//    // Get a reference to our system.
//    EigenSystem & eigen_system = es.get_system<EigenSystem> (system_name);
//    
//    // Get a constant reference to the Finite Element type
//    // for the first (and only) variable in the system.
//    FEType fe_type = eigen_system.get_dof_map().variable_type(0);
//    
//    const DofMap& dof_map = eigen_system.get_dof_map();
//    
//    std::vector<dof_id_type> dof_indices;
//    
//    
//    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
//    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
//    
//    for ( ; el != end_el; ++el)
//    {
//        dof_map.dof_indices (*el, dof_indices, 1); // uy
//        
//        // All boundary dofs are Dirichlet dofs in this case
//        for (unsigned int s=0; s<(*el)->n_sides(); s++)
//            if ((*el)->neighbor(s) == NULL)
//            {
//                std::vector<unsigned int> side_dofs;
//                FEInterface::dofs_on_side(*el, dim, fe_type,
//                                          s, side_dofs);
//                
//                for(unsigned int ii=0; ii<side_dofs.size(); ii++)
//                    dirichlet_dof_ids.insert(dof_indices[side_dofs[ii]]);
//            }
//        
//        // also add the dofs for variable u, w, tx, ty
//        dof_indices.clear();
//        dof_map.dof_indices(*el, dof_indices, 0); // ux
//        for (unsigned int i=0; i<dof_indices.size(); i++)
//            dirichlet_dof_ids.insert(dof_indices[i]);
//        
//        dof_indices.clear();
//        dof_map.dof_indices(*el, dof_indices, 2); // uw
//        for (unsigned int i=0; i<dof_indices.size(); i++)
//            dirichlet_dof_ids.insert(dof_indices[i]);
//        
//        dof_indices.clear();
//        dof_map.dof_indices(*el, dof_indices, 3); // tx
//        for (unsigned int i=0; i<dof_indices.size(); i++)
//            dirichlet_dof_ids.insert(dof_indices[i]);
//
//        dof_indices.clear();
//        dof_map.dof_indices(*el, dof_indices, 4); // ty
//        for (unsigned int i=0; i<dof_indices.size(); i++)
//            dirichlet_dof_ids.insert(dof_indices[i]);
//    } // end of element loop
//    
//    /**
//     * All done!
//     */
//    return;
//    
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






