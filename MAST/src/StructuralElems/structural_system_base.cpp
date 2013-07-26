
// MAST includes
#include "StructuralElems/structural_system_base.h"

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


// FESystem includes
#include "Quadrature/TrapezoidQuadrature.h"
#include "Mesh/Tri3.h"
#include "Mesh/Node.h"
#include "Geom/RectangularCoordinateSystem.h"
#include "FiniteElems/FELagrange.h"
#include "Numerics/DenseMatrix.h"
#include "Numerics/LocalVector.h"
#include "Disciplines/Structure/DKTPlate.h"

#ifndef LIBMESH_USE_COMPLEX_NUMBERS


void StructuralSystemBase::init_data ()
{
    this->use_fixed_solution = true;
    
    unsigned int dim = this->get_mesh().mesh_dimension();
    
    vars.resize(6);

    
    GetPot infile("structural.in");
    
    unsigned int o = infile("fe_order", 1);
    std::string fe_family = infile("fe_family", std::string("LAGRANGE"));
    FEFamily fefamily = Utility::string_to_enum<FEFamily>(fe_family);
    
    vars[0]  = this->add_variable ( "ux", static_cast<Order>(o), fefamily);
    this->time_evolving(vars[0]);
    vars[1]  = this->add_variable ( "uy", static_cast<Order>(o), fefamily);
    this->time_evolving(vars[1]);
    vars[2]  = this->add_variable ( "uz", static_cast<Order>(o), fefamily);
    this->time_evolving(vars[2]);
    vars[3]  = this->add_variable ( "tx", static_cast<Order>(o), fefamily);
    this->time_evolving(vars[3]);
    vars[4]  = this->add_variable ( "ty", static_cast<Order>(o), fefamily);
    this->time_evolving(vars[4]);
    vars[5]  = this->add_variable ( "tz", static_cast<Order>(o), fefamily);
    this->time_evolving(vars[5]);
    
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
        ZeroFunction<Real> zero_function;
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
            nodes[i]->setVal(j, c.elem->point(i)(j));
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
    dkt_plate.calculateStiffnessMatrix(plate_elem_mat);
    dkt_plate.transformMatrixToGlobalSystem(plate_elem_mat, elem_mat);

    // set small values for the u, v, and tz dof stiffness
    for (unsigned int i=0; i<3; i++)
    {
        elem_mat.setVal(   i,     i, 1.0e-6);
        elem_mat.setVal( i+3,   i+3, 1.0e-6);
        elem_mat.setVal(i+15,  i+15, 1.0e-6);
    }
    
    
    // force
    dkt_plate.calculateDistributedLoad(1.0e6, plate_elem_vec);
    dkt_plate.transformVectorToGlobalSystem(plate_elem_vec, elem_vec);
    
    
    
    // mass
//    dkt_plate.calculateDiagonalMassMatrix(plate_elem_vec);
//    dkt_plate.getActiveElementMatrixIndices(elem_dof_indices);
    

    //    for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
        DenseMatrix<Real> tmpmat; tmpmat.resize(18, 18);
        DenseVector<Real> tmp; tmp.resize(18);
        
        for (unsigned int i=0; i<18; i++)
            for (unsigned int j=0; j<18; j++)
                tmpmat(i, j) = elem_mat.getVal(i, j);

        tmpmat.vector_mult(tmp, c.elem_solution);
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
    
    
//    for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
        if (request_jacobian && c.elem_solution_derivative)
        {
        }
    } // end of the quadrature point qp-loop
    
    //    std::cout << "inside mass residual " << std::endl;
    //    std::cout << "elem velocity" << std::endl; c.elem_solution.print(std::cout);
    //    std::cout << "elem solution" << std::endl; c.elem_fixed_solution.print(std::cout);
    //    std::cout << "mass vec: " << std::endl; Fvec.print(std::cout);
    //    if (request_jacobian && c.elem_solution_derivative)
    //        Kmat.print(std::cout);
    
    return request_jacobian;
}


#endif // LIBMESH_USE_COMPLEX_NUMBERS
