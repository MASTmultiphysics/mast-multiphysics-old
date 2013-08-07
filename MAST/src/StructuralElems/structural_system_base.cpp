
// MAST includes
#include "StructuralElems/structural_system_base.h"
#include "StructuralElems/surface_pressure_load.h"

// libMesh include files
#include "libmesh/equation_systems.h"
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/parameters.h"
#include "libmesh/quadrature.h"
#include "libmesh/dof_map.h"
#include "libmesh/zero_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/mesh_function.h"
#include "libmesh/fem_function_base.h"
#include "libmesh/eigen_system.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/fe_interface.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/fe.h"

// FESystem includes
#include "Quadrature/TrapezoidQuadrature.h"
#include "Mesh/Tri3.h"
#include "Mesh/Edge2.h"
#include "Mesh/Node.h"
#include "Geom/RectangularCoordinateSystem.h"
#include "FiniteElems/FELagrange.h"
#include "Numerics/DenseMatrix.h"
#include "Numerics/LocalVector.h"
#include "Disciplines/Structure/DKTPlate.h"
#include "Disciplines/Structure/EulerBernoulliBeam.h"



void StructuralSystemBase::init_data ()
{
    this->use_fixed_solution = true;
    
    unsigned int dim = this->get_mesh().mesh_dimension();
    
    vars.resize(6);

    
    GetPot infile("structural.in");
    
    unsigned int o = infile("fe_order", 1);
    std::string fe_family = infile("fe_family", std::string("LAGRANGE"));
    FEFamily fefamily = Utility::string_to_enum<FEFamily>(fe_family);
    
    vars.push_back(this->add_variable ( "ux", static_cast<Order>(o), fefamily));
    vars.push_back(this->add_variable ( "uy", static_cast<Order>(o), fefamily));
    vars.push_back(this->add_variable ( "uz", static_cast<Order>(o), fefamily));
    vars.push_back(this->add_variable ( "tx", static_cast<Order>(o), fefamily));
    vars.push_back(this->add_variable ( "ty", static_cast<Order>(o), fefamily));
    vars.push_back(this->add_variable ( "tz", static_cast<Order>(o), fefamily));
    
    for (unsigned int i=0; i<vars.size(); i++)
        this->time_evolving(vars[i]);
    
    // Useful debugging options
    this->verify_analytic_jacobians = infile("verify_analytic_jacobians", 0.);
    this->print_jacobians = infile("print_jacobians", false);
    this->print_element_jacobians = infile("print_element_jacobians", false);
    
    // read in the boudnary conditions
    unsigned int n_bc, bc_id;
    // first the inviscid conditions
    // slip wall bc
    n_bc = infile("n_dirichlet_bc", 0);
    if (n_bc > 0)
    {
        for (unsigned int i_bc=0; i_bc<n_bc; i_bc++)
        {
            bc_id = infile("dirichlet_bc", 0, i_bc);
            _boundary_condition.insert
            (std::multimap<unsigned int, StructuralBoundaryConditionType>::value_type
             (bc_id, DIRICHLET));
        }
    }
    
    
    // initialize the Dirichlet boundary conditions
    std::multimap<unsigned, StructuralBoundaryConditionType>::const_iterator
    bc_it = _boundary_condition.begin(), bc_end = _boundary_condition.end();
    
    std::set<boundary_id_type> dirichlet_boundary;
    
    for ( ; bc_it!=bc_end; bc_it++)
    {
        switch (bc_it->second)
        {
            case DIRICHLET: // u_i = 0, tx_i = 0;
                dirichlet_boundary.insert(bc_it->first);
                break;
                
            default:
                break;
        }
    }
    
    if (dirichlet_boundary.size() > 0)
    {
        // Dirichlet Boundary condition
        ZeroFunction<Number> zero_function;
        this->get_dof_map().add_dirichlet_boundary(DirichletBoundary(dirichlet_boundary,
                                                                     vars,
                                                                     &zero_function));
    }
    
    // ask the parent class to initilize the data
    FEMSystem::init_data();
}



void StructuralSystemBase::init_context(DiffContext &context)
{
    FEMContext &c = libmesh_cast_ref<FEMContext&>(context);
    
    std::vector<FEBase*> elem_fe(dim+2);
    
    for (unsigned int i=0; i<dim+2; i++)
    {
        c.get_element_fe( vars[i], elem_fe[i]);
        elem_fe[i]->get_JxW();
        elem_fe[i]->get_phi();
        elem_fe[i]->get_dphi();
        elem_fe[i]->get_xyz();
    }
    
    std::vector<FEBase*> elem_side_fe(dim+2);
    
    for (unsigned int i=0; i<dim+2; i++)
    {
        c.get_side_fe( vars[i], elem_side_fe[i]);
        elem_side_fe[i]->get_JxW();
        elem_side_fe[i]->get_phi();
        elem_side_fe[i]->get_xyz();
        elem_side_fe[i]->get_dphi();
    }
}



bool StructuralSystemBase::element_time_derivative (bool request_jacobian,
                                                    DiffContext &context)
{
    FEMContext &c = libmesh_cast_ref<FEMContext&>(context);
    DenseMatrix<Number>& Kmat = c.get_elem_jacobian();
    DenseVector<Number>& Fvec = c.get_elem_residual();

    FESystem::Quadrature::TrapezoidQuadrature q_rule_shear, q_rule_bending;
    FESystem::FiniteElement::FELagrange fe, fe_tri6;
    FESystem::Structures::DKTPlate dkt_plate;
    q_rule_bending.init(2, 9);
    
    // initialize the geometric element
    std::auto_ptr<FESystem::Mesh::ElemBase> elem(new FESystem::Mesh::Tri3(false));

    FESystem::Numerics::DenseMatrix<Real> basis; basis.resize(3, 3); basis.setToIdentity();
    FESystem::Geometry::Point origin(3);
    FESystem::Geometry::RectangularCoordinateSystem cs(origin, basis);
    std::vector<FESystem::Mesh::Node*> nodes(3);
    
    for (unsigned int i=0; i<3; i++)
    {
        nodes[i] = new FESystem::Mesh::Node(cs);
        for (unsigned int j=0; j<3; j++)
            nodes[i]->setVal(j, c.get_elem().point(i)(j));
        elem->setNode(i, *nodes[i]);
    }
    fe.reinit(*elem);
    
    dkt_plate.initialize(*elem, fe, fe_tri6, q_rule_bending, q_rule_bending,
                         72.0e9, 0.33, 2700., 0.002);

    FESystem::Numerics::DenseMatrix<Real> plate_elem_mat, elem_mat;
    FESystem::Numerics::LocalVector<Real> plate_elem_vec, elem_vec;
    plate_elem_mat.resize(dkt_plate.getNElemDofs(), dkt_plate.getNElemDofs());
    elem_mat.resize(18, 18);
    plate_elem_vec.resize(dkt_plate.getNElemDofs());
    elem_vec.resize(18);
    
    // stiffness matrix
    dkt_plate.calculateStiffnessMatrix(plate_elem_mat);
    dkt_plate.transformMatrixToGlobalSystem(plate_elem_mat, elem_mat);

    // set small values for the u, v, and tz dof stiffness
    for (unsigned int i=0; i<3; i++)
    {
        elem_mat.setVal(   i,     i, 1.0e16);
        elem_mat.setVal( i+3,   i+3, 1.0e16);
        elem_mat.setVal(i+15,  i+15, 1.0e16);
    }
    
    
    // force
    dkt_plate.calculateDistributedLoad(1.0e6, plate_elem_vec);
    dkt_plate.transformVectorToGlobalSystem(plate_elem_vec, elem_vec);
    
    
    //    for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
        DenseMatrix<Number> tmpmat; tmpmat.resize(18, 18);
        DenseVector<Number> tmp; tmp.resize(18);
        
        for (unsigned int i=0; i<18; i++)
            for (unsigned int j=0; j<18; j++)
                tmpmat(i, j) = elem_mat.getVal(i, j);

        tmpmat.vector_mult(tmp, c.get_elem_solution());
        Fvec.add(1.0, tmp);
        
        for (unsigned int i=0; i<18; i++)
            Fvec(i) -= elem_vec.getVal(i);
        
        if (request_jacobian && c.elem_solution_derivative)
        {
            // now copy the matrix values to the stiffness mat
            for (unsigned int i=0; i<18; i++)
                for (unsigned int j=0; j<18; j++)
                    Kmat(i, j) += elem_mat.getVal(i, j);
        }
    } // end of the quadrature point qp-loop
        
    //    std::cout << "inside element time derivative " << std::endl;
    //    c.elem->print_info();
    //    std::cout << "sol: " << std::endl; c.elem_solution.print(std::cout);
    //    std::cout << "res: " << std::endl; Fvec.print(std::cout);
    //    if (request_jacobian && c.elem_solution_derivative)
    //        Kmat.print(std::cout);
    
    // clear the pointers
    elem.reset();
    for (unsigned int i=0; i<3; i++)
        delete nodes[i];
        
    return request_jacobian;
}



bool StructuralSystemBase::side_time_derivative (bool request_jacobian,
                                        DiffContext &context)
{
    
    return request_jacobian;
}





bool StructuralSystemBase::mass_residual (bool request_jacobian,
                                          DiffContext &context)
{
    FEMContext &c = libmesh_cast_ref<FEMContext&>(context);
    DenseMatrix<Number>& Kmat = c.get_elem_jacobian();
    DenseVector<Number>& Fvec = c.get_elem_residual();
    
    FESystem::Quadrature::TrapezoidQuadrature q_rule_shear, q_rule_bending;
    FESystem::FiniteElement::FELagrange fe, fe_tri6;
    FESystem::Structures::DKTPlate dkt_plate;
    q_rule_bending.init(2, 9);
    
    // initialize the geometric element
    std::auto_ptr<FESystem::Mesh::ElemBase> elem(new FESystem::Mesh::Tri3(false));
    
    FESystem::Numerics::DenseMatrix<Real> basis; basis.resize(3, 3); basis.setToIdentity();
    FESystem::Geometry::Point origin(3);
    FESystem::Geometry::RectangularCoordinateSystem cs(origin, basis);
    std::vector<FESystem::Mesh::Node*> nodes(3);
    
    for (unsigned int i=0; i<3; i++)
    {
        nodes[i] = new FESystem::Mesh::Node(cs);
        for (unsigned int j=0; j<3; j++)
            nodes[i]->setVal(j, c.get_elem().point(i)(j));
        elem->setNode(i, *nodes[i]);
    }
    fe.reinit(*elem);
    
    dkt_plate.initialize(*elem, fe, fe_tri6, q_rule_bending, q_rule_bending, 72.0e9, 0.33, 2700., 0.002);
    
    FESystem::Numerics::DenseMatrix<Real> plate_elem_mat, elem_mat;
    FESystem::Numerics::LocalVector<Real> plate_elem_vec, elem_vec;
    plate_elem_mat.resize(dkt_plate.getNElemDofs(), dkt_plate.getNElemDofs());
    elem_mat.resize(18, 18);
    plate_elem_vec.resize(dkt_plate.getNElemDofs());
    elem_vec.resize(18);
    
    // stiffness matrix
    dkt_plate.calculateConsistentMassMatrix(plate_elem_mat);
    dkt_plate.transformMatrixToGlobalSystem(plate_elem_mat, elem_mat);
    
    // set small values for the u, v, and tz dof mass
    for (unsigned int i=0; i<3; i++)
    {
        elem_mat.setVal(   i,     i, 1.0e-6);
        elem_mat.setVal( i+3,   i+3, 1.0e-6);
        elem_mat.setVal(i+15,  i+15, 1.0e-6);
    }
    
    //    for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
        DenseMatrix<Number> tmpmat; tmpmat.resize(18, 18);
        DenseVector<Number> tmp; tmp.resize(18);
        
        for (unsigned int i=0; i<18; i++)
            for (unsigned int j=0; j<18; j++)
                tmpmat(i, j) = elem_mat.getVal(i, j);
        
        tmpmat.vector_mult(tmp, c.get_elem_solution());
        Fvec.add(1.0, tmp);
        
        if (request_jacobian && c.elem_solution_derivative)
        {
            // now copy the matrix values to the stiffness mat
            for (unsigned int i=0; i<18; i++)
                for (unsigned int j=0; j<18; j++)
                    Kmat(i, j) += elem_mat.getVal(i, j);
        }
    } // end of the quadrature point qp-loop
    
    // clear the pointers
    elem.reset();
    for (unsigned int i=0; i<3; i++)
        delete nodes[i];
    
    
    //    std::cout << "inside mass residual " << std::endl;
    //    std::cout << "elem velocity" << std::endl; c.elem_solution.print(std::cout);
    //    std::cout << "elem solution" << std::endl; c.elem_fixed_solution.print(std::cout);
    //    std::cout << "mass vec: " << std::endl; Fvec.print(std::cout);
    //    if (request_jacobian && c.elem_solution_derivative)
    //        Kmat.print(std::cout);
    
    return request_jacobian;
}



void assemble_plate_matrices(EquationSystems& es,
                             const std::string& system_name)
{
    // Get a constant reference to the mesh object.
    const MeshBase& mesh = es.get_mesh();
    
    // The dimension that we are running.
    const unsigned int dim = mesh.mesh_dimension();
    
    // Get a reference to our system.
    EigenSystem & eigen_system = es.get_system<EigenSystem> (system_name);
    
    // A reference to the two system matrices
    SparseMatrix<Number>&  matrix_A = *eigen_system.matrix_A;
    SparseMatrix<Number>&  matrix_B = *eigen_system.matrix_B;
    
    // A reference to the \p DofMap object for this system.  The \p DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.
    const DofMap& dof_map = eigen_system.get_dof_map();
    
    // The element mass and stiffness matrices.
    DenseMatrix<Number>   Me;
    DenseMatrix<Number>   Ke;
    
    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector<dof_id_type> dof_indices;
    
    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    
    for ( ; el != end_el; ++el)
    {
        dof_map.dof_indices (*el, dof_indices);
        Ke.resize (dof_indices.size(), dof_indices.size());
        Me.resize (dof_indices.size(), dof_indices.size());
        
        FESystem::Quadrature::TrapezoidQuadrature q_rule_shear, q_rule_bending;
        FESystem::FiniteElement::FELagrange fe, fe_tri6;
        FESystem::Structures::DKTPlate dkt_plate;
        q_rule_bending.init(2, 9);
        
        // initialize the geometric element
        std::auto_ptr<FESystem::Mesh::ElemBase> elem(new FESystem::Mesh::Tri3(false));
        
        FESystem::Numerics::DenseMatrix<Real> basis; basis.resize(3, 3); basis.setToIdentity();
        FESystem::Geometry::Point origin(3);
        FESystem::Geometry::RectangularCoordinateSystem cs(origin, basis);
        std::vector<FESystem::Mesh::Node*> nodes(3);
        
        for (unsigned int i=0; i<3; i++)
        {
            nodes[i] = new FESystem::Mesh::Node(cs);
            for (unsigned int j=0; j<3; j++)
                nodes[i]->setVal(j, (*el)->point(i)(j));
            elem->setNode(i, *nodes[i]);
        }
        fe.reinit(*elem);
        
        dkt_plate.initialize(*elem, fe, fe_tri6, q_rule_bending, q_rule_bending,
                             72.0e9, 0.33, 2700., 0.002);
        
        FESystem::Numerics::DenseMatrix<Real> plate_elem_mat, elem_mat;
        FESystem::Numerics::LocalVector<Real> plate_elem_vec, elem_vec;
        plate_elem_mat.resize(dkt_plate.getNElemDofs(), dkt_plate.getNElemDofs());
        elem_mat.resize(18, 18);
        plate_elem_vec.resize(dkt_plate.getNElemDofs());
        elem_vec.resize(18);
        
        dkt_plate.calculateStiffnessMatrix(plate_elem_mat);
        dkt_plate.transformMatrixToGlobalSystem(plate_elem_mat, elem_mat);
        for (unsigned int i=0; i<18; i++)
            for (unsigned int j=0; j<18; j++)
                Ke(i, j) = elem_mat.getVal(i, j);
        
        dkt_plate.calculateConsistentMassMatrix(plate_elem_mat);
        // put small values on the diagonal for rotation dofs
        for (unsigned int i=0; i<6; i++)
            plate_elem_mat.setVal(3+i, 3+i, 1.0);
        dkt_plate.transformMatrixToGlobalSystem(plate_elem_mat, elem_mat);
        for (unsigned int i=0; i<18; i++)
            for (unsigned int j=0; j<18; j++)
                Me(i, j) = elem_mat.getVal(i, j);
        
        // clear the pointers
        elem.reset();
        for (unsigned int i=0; i<3; i++)
            delete nodes[i];
        
        matrix_A.add_matrix (Me, dof_indices);
        matrix_B.add_matrix (Ke, dof_indices);
    }
}



void assemble_beam_matrices(EquationSystems& es,
                            const std::string& system_name)
{
    // Get a constant reference to the mesh object.
    const MeshBase& mesh = es.get_mesh();
    
    // The dimension that we are running.
    const unsigned int dim = mesh.mesh_dimension();
    
    // Get a reference to our system.
    EigenSystem & eigen_system = es.get_system<EigenSystem> (system_name);
    
    // A reference to the two system matrices
    SparseMatrix<Number>&  matrix_A = *eigen_system.matrix_A;
    SparseMatrix<Number>&  matrix_B = *eigen_system.matrix_B;
    
    // A reference to the \p DofMap object for this system.  The \p DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.
    const DofMap& dof_map = eigen_system.get_dof_map();
    
    // The element mass and stiffness matrices.
    DenseMatrix<Number>   Me;
    DenseMatrix<Number>   Ke;
    
    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector<dof_id_type> dof_indices;
    
    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    
    for ( ; el != end_el; ++el)
    {
        dof_map.dof_indices (*el, dof_indices);
        Ke.resize (dof_indices.size(), dof_indices.size());
        Me.resize (dof_indices.size(), dof_indices.size());
        
        FESystem::Quadrature::TrapezoidQuadrature q_rule_shear, q_rule_bending;
        FESystem::FiniteElement::FELagrange fe, fe_tri6;
        FESystem::Structures::EulerBernoulliBeam beam;
        q_rule_bending.init(1, 9);
        
        // initialize the geometric element
        std::auto_ptr<FESystem::Mesh::Edge2> elem(new FESystem::Mesh::Edge2(false));
        
        FESystem::Numerics::DenseMatrix<Real> basis; basis.resize(3, 3); basis.setToIdentity();
        FESystem::Geometry::Point origin(3);
        FESystem::Geometry::RectangularCoordinateSystem cs(origin, basis);
        FESystem::Numerics::LocalVector<Real> yvec; yvec.resize(3); yvec.setVal(1, 1.);
        std::vector<FESystem::Mesh::Node*> nodes(2);
        
        elem->setVectorForXYPlane(yvec);
        for (unsigned int i=0; i<2; i++)
        {
            nodes[i] = new FESystem::Mesh::Node(cs);
            for (unsigned int j=0; j<2; j++)
                nodes[i]->setVal(j, (*el)->point(i)(j));
            elem->setNode(i, *nodes[i]);
        }
        fe.reinit(*elem);
        
        beam.initialize(*elem, fe, q_rule_bending,
                        72.0e9, 0.33, 2700., 6.667e-9, 1.6667e-9, 2.0e-4);
        
        FESystem::Numerics::DenseMatrix<Real> beam_elem_mat, elem_mat;
        FESystem::Numerics::LocalVector<Real> beam_elem_vec, elem_vec;
        beam_elem_mat.resize(beam.getNElemDofs(), beam.getNElemDofs());
        elem_mat.resize(12, 12);
        beam_elem_vec.resize(beam.getNElemDofs());
        elem_vec.resize(12);
        
        beam.calculateStiffnessMatrix(beam_elem_mat);
        beam.transformMatrixToGlobalSystem(beam_elem_mat, elem_mat);
        for (unsigned int i=0; i<12; i++)
            for (unsigned int j=0; j<12; j++)
                Ke(i, j) = elem_mat.getVal(i, j);
        
        beam.calculateConsistentMassMatrix(beam_elem_mat);
        beam.transformMatrixToGlobalSystem(beam_elem_mat, elem_mat);
        for (unsigned int i=0; i<12; i++)
            for (unsigned int j=0; j<12; j++)
                Me(i, j) = elem_mat.getVal(i, j);
        
        // clear the pointers
        elem.reset();
        for (unsigned int i=0; i<2; i++)
            delete nodes[i];
        
        matrix_A.add_matrix (Me, dof_indices);
        matrix_B.add_matrix (Ke, dof_indices);
    }
}




void get_plate_dirichlet_dofs(EquationSystems& es,
                              const std::string& system_name,
                              std::set<unsigned int>& dirichlet_dof_ids)
{
    dirichlet_dof_ids.clear();
    
    // Get a constant reference to the mesh object.
    const MeshBase& mesh = es.get_mesh();
    
    // The dimension that we are running.
    const unsigned int dim = mesh.mesh_dimension();
    
    // Get a reference to our system.
    EigenSystem & eigen_system = es.get_system<EigenSystem> (system_name);
    
    // Get a constant reference to the Finite Element type
    // for the first (and only) variable in the system.
    FEType fe_type = eigen_system.get_dof_map().variable_type(0);
    
    const DofMap& dof_map = eigen_system.get_dof_map();
    
    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector<dof_id_type> dof_indices;
    
    
    // Now we will loop over all the elements in the mesh that
    // live on the local processor. We will compute the element
    // matrix and right-hand-side contribution.  In case users
    // later modify this program to include refinement, we will
    // be safe and will only consider the active elements;
    // hence we use a variant of the \p active_elem_iterator.
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
                    dirichlet_dof_ids.insert(dof_indices[side_dofs[ii]]);
            }
        
        // also add the dofs for variable u, v and tz
        dof_indices.clear();
        dof_map.dof_indices(*el, dof_indices, 0); // ux
        for (unsigned int i=0; i<dof_indices.size(); i++)
            dirichlet_dof_ids.insert(dof_indices[i]);
        
        dof_indices.clear();
        dof_map.dof_indices(*el, dof_indices, 1); // uy
        for (unsigned int i=0; i<dof_indices.size(); i++)
            dirichlet_dof_ids.insert(dof_indices[i]);
        
        dof_indices.clear();
        dof_map.dof_indices(*el, dof_indices, 5); // tz
        for (unsigned int i=0; i<dof_indices.size(); i++)
            dirichlet_dof_ids.insert(dof_indices[i]);
    } // end of element loop

    /**
     * All done!
     */
    return;
    
}



void get_beam_dirichlet_dofs(EquationSystems& es,
                             const std::string& system_name,
                             std::set<unsigned int>& dirichlet_dof_ids)
{
    dirichlet_dof_ids.clear();
    
    // Get a constant reference to the mesh object.
    const MeshBase& mesh = es.get_mesh();
    
    // The dimension that we are running.
    const unsigned int dim = mesh.mesh_dimension();
    
    // Get a reference to our system.
    EigenSystem & eigen_system = es.get_system<EigenSystem> (system_name);
    
    // Get a constant reference to the Finite Element type
    // for the first (and only) variable in the system.
    FEType fe_type = eigen_system.get_dof_map().variable_type(0);
    
    const DofMap& dof_map = eigen_system.get_dof_map();
    
    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector<dof_id_type> dof_indices;
    
    
    // Now we will loop over all the elements in the mesh that
    // live on the local processor. We will compute the element
    // matrix and right-hand-side contribution.  In case users
    // later modify this program to include refinement, we will
    // be safe and will only consider the active elements;
    // hence we use a variant of the \p active_elem_iterator.
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
                    dirichlet_dof_ids.insert(dof_indices[side_dofs[ii]]);
            }
        
        // also add the dofs for variable u, w, tx, ty
        dof_indices.clear();
        dof_map.dof_indices(*el, dof_indices, 0); // ux
        for (unsigned int i=0; i<dof_indices.size(); i++)
            dirichlet_dof_ids.insert(dof_indices[i]);
        
        dof_indices.clear();
        dof_map.dof_indices(*el, dof_indices, 2); // uw
        for (unsigned int i=0; i<dof_indices.size(); i++)
            dirichlet_dof_ids.insert(dof_indices[i]);
        
        dof_indices.clear();
        dof_map.dof_indices(*el, dof_indices, 3); // tx
        for (unsigned int i=0; i<dof_indices.size(); i++)
            dirichlet_dof_ids.insert(dof_indices[i]);

        dof_indices.clear();
        dof_map.dof_indices(*el, dof_indices, 4); // ty
        for (unsigned int i=0; i<dof_indices.size(); i++)
            dirichlet_dof_ids.insert(dof_indices[i]);
    } // end of element loop
    
    /**
     * All done!
     */
    return;
    
}



void assemble_beam_force_vec(System& sys,
                             SurfacePressureLoad& surf_press,
                             SurfaceMotionBase& surf_motion,
                             NumericVector<Number>& fvec)
{
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
    // Get a constant reference to the mesh object.
    const MeshBase& mesh = sys.get_mesh();
    
    // The dimension that we are running.
    const unsigned int dim = mesh.mesh_dimension();
    
    // A reference to the \p DofMap object for this system.  The \p DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.
    const DofMap& dof_map = sys.get_dof_map();
    std::vector<dof_id_type> dof_indices;

    // The element mass and stiffness matrices.
    DenseVector<Number>   fvec_e;
    
    FEType fe_type(FIRST, LAGRANGE);
    AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));
    QGauss qrule (dim, FIFTH);
    fe->attach_quadrature_rule (&qrule);
    const std::vector<Real>& JxW = fe->get_JxW();
    const std::vector<Point>& q_point = fe->get_xyz();
    const std::vector<std::vector<Real> >& phi = fe->get_phi();
    
    
    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    Number press, dpress;
    DenseVector<Number> utrans, dn_rot; utrans.resize(3); dn_rot.resize(3);
    Point normal; normal.zero(); normal(1) = 1.; // y_vec
    
    fvec.zero(); fvec.close();
    
    for ( ; el != end_el; ++el)
    {
        dof_map.dof_indices (*el, dof_indices);
        fvec_e.resize (dof_indices.size());
        fe->reinit (*el);
        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
            surf_press.surface_pressure(q_point[qp], press, dpress);
            surf_motion.surface_velocity_frequency_domain(q_point[qp], normal,
                                                          utrans, dn_rot);
            for (unsigned int i=0; i<3; i++)
                for (unsigned int iphi=0; iphi<phi.size(); iphi++)
                    fvec_e(i*phi.size()+iphi) += JxW[qp] * phi[iphi][qp] *
                    ( press * dn_rot(i) + // steady pressure
                     dpress * normal(i)); // unsteady pressure
        }
        
        fvec.add_vector (fvec_e, dof_indices);
    }
    fvec.close();
#endif // LIBMESH_USE_COMPLEX_NUMBERS
}





