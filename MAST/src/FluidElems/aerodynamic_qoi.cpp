//
//  aerodynamic_qoi.cpp
//  MAST
//
//  Created by Manav Bhatia on 6/11/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//


// MAST includes
#include "FluidElems/aerodynamic_qoi.h"




void AerodynamicQoI::init_qoi( std::vector<Real>& sys_qoi)
{
    //Two qois are calculated: lift and drag
    sys_qoi.resize(4);
    return;
}




void AerodynamicQoI::element_qoi_derivative (libMesh::DiffContext& context,
                                             const libMesh::QoISet& qois)
{
    
}



void AerodynamicQoI::element_qoi (libMesh::DiffContext& context, const libMesh::QoISet& qois)
{

    
    libMesh::FEMContext &c = libMesh::libmesh_cast_ref<libMesh::FEMContext&>(context);

    libMesh::FEBase* elem_fe;
    
    c.get_element_fe(_fluid_vars[0], elem_fe);
    
    // Element Jacobian * quadrature weights for interior integration
    const std::vector<Real> &JxW = elem_fe->get_JxW();
    
    // Physical location of the quadrature points
    const std::vector<libMesh::Point>& qpoint = elem_fe->get_xyz();
    
    DenseRealVector conservative_sol, integrated_force;
    FEMOperatorMatrix  B_mat;
    integrated_force.resize(dim); conservative_sol.resize(dim+2);
    
    std::vector<FEMOperatorMatrix> dB_mat(dim);
    
    PrimitiveSolution p_sol;

    Real total_volume=0., entropy_error=0.;
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++)
    {
        // first update the variables at the current quadrature point
        this->update_solution_at_quadrature_point
        (_fluid_vars, qp, c,
         true, c.get_elem_solution(),
         conservative_sol, p_sol,
         B_mat, dB_mat);
        
        
        total_volume += JxW[qp];
        entropy_error +=
        JxW[qp] * pow(log(p_sol.p/flight_condition->p0() *
                       pow(flight_condition->rho()/p_sol.rho,
                           flight_condition->gas_property.gamma)), 2);
    }
    
    std::vector<Real> vals(2);
    vals[0] = total_volume;
    vals[1] = entropy_error;

    // the third qoi is total volume, and the fourth is the entropy
    for (unsigned int i=2; i<4; i++)
        if (qois.has_index(i))
            c.get_qois()[i] += vals[i-2];

}



void AerodynamicQoI::side_qoi_derivative (libMesh::DiffContext &context,
                                          const libMesh::QoISet& qois)
{
    libmesh_error(); // to be implemented
}




// This function computes the actual QoI
void AerodynamicQoI::side_qoi(libMesh::DiffContext &context, const libMesh::QoISet& qois)
{


    libMesh::FEMContext &c =
    libMesh::libmesh_cast_ref<libMesh::FEMContext&>(context);
    
    // check for the boundary tags to check if the boundary condition needs to be applied to this element
    
    FluidBoundaryConditionType mechanical_bc_type;
    
    unsigned int n_mechanical_bc = 0;
    
    std::multimap<unsigned int, FluidBoundaryConditionType>::const_iterator
    bc_it     = this->_boundary_condition.begin(),
    bc_end    = this->_boundary_condition.end();
    
    
    for ( ; bc_it != bc_end; bc_it++)
    {
        if ( c.has_side_boundary_id(bc_it->first))
        {
            switch (bc_it->second)
            {
                case SLIP_WALL:
                    mechanical_bc_type = SLIP_WALL;
                    n_mechanical_bc++;
                    break;
                    
                case NO_SLIP_WALL:
                    libmesh_assert(_if_viscous);
                    mechanical_bc_type = NO_SLIP_WALL;
                    n_mechanical_bc++;
                    break;
                    
                default:
                    // nothing to be done for other BC kinds
                    break;
            }
        }
    }
    
    // return if no boundary condition is applied
    if (n_mechanical_bc == 0)
        return;
    
    
    libmesh_assert_equal_to(n_mechanical_bc, 1); // make sure that only one is active
    
    libMesh::FEBase * side_fe;
    c.get_side_fe(_fluid_vars[0], side_fe); // assuming all variables have the same FE
    
    // Element Jacobian * quadrature weights for interior integration
    const std::vector<Real> &JxW = side_fe->get_JxW();
    
    // Physical location of the quadrature points
    const std::vector<libMesh::Point>& qpoint = side_fe->get_xyz();
    
    // boundary normals
    const std::vector<libMesh::Point>& face_normals = side_fe->get_normals();
    
    DenseRealVector conservative_sol;
    FEMOperatorMatrix  B_mat;
    conservative_sol.resize(dim+2);
    
    RealVector3 integrated_force;
    
    std::vector<FEMOperatorMatrix> dB_mat(dim);
    
    PrimitiveSolution p_sol;
    
    switch (mechanical_bc_type)
    {
        case NO_SLIP_WALL: // both slip and no-slip wall have the same force integration
        case SLIP_WALL:
        {
            for (unsigned int qp=0; qp<qpoint.size(); qp++)
            {
                // first update the variables at the current quadrature point
                this->update_solution_at_quadrature_point
                (_fluid_vars, qp, c,
                 false, c.get_elem_solution(),
                 conservative_sol, p_sol,
                 B_mat, dB_mat);
                
                // update the force vector
                for (unsigned int i=0; i<dim; i++)
                    integrated_force(i) += face_normals[qp](i)*JxW[qp]*p_sol.p;
            }
        }
            break;
            
        default:
            //nothing to be done for other bc kinds
            break;
    }
    
    std::vector<Real> vals(2);
    vals[0] = integrated_force.dot(flight_condition->lift_normal);
    vals[1] = integrated_force.dot(flight_condition->drag_normal);
    
    // the first qoi is lift, and the second is drag
    for (unsigned int i=0; i<2; i++)
        if (qois.has_index(i))
            c.get_qois()[i] += vals[i];
    

}

