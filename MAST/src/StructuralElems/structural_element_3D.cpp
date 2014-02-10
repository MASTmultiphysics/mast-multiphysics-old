//
//  structural_element_3D.cpp
//  MAST
//
//  Created by Manav Bhatia on 11/14/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

// MAST includes
#include "StructuralElems/structural_element_3D.h"
#include "PropertyCards/element_property_card_base.h"
#include "Numerics/fem_operator_matrix.h"
#include "ThermalElems/temperature_function.h"
#include "BoundaryConditions/boundary_condition.h"


bool
MAST::StructuralElement3D::internal_force (bool request_jacobian,
                                           DenseVector<Real>& f,
                                           DenseMatrix<Real>& jac)
{
    FEMOperatorMatrix Bmat;
    
    const std::vector<Real>& JxW = _fe->get_JxW();
    const std::vector<Point>& xyz = _fe->get_xyz();
    const unsigned int n_phi = (unsigned int)JxW.size();
    const unsigned int n1=6, n2=6*n_phi;
    DenseMatrix<Real> material_mat, tmp_mat1_n1n2, tmp_mat2_n2n2;
    DenseVector<Real>  phi, tmp_vec1_n1, tmp_vec2_n2;
    
    tmp_mat1_n1n2.resize(n1, n2); tmp_mat2_n2n2.resize(n2, n2);
    phi.resize(n_phi); tmp_vec1_n1.resize(n1); tmp_vec2_n2.resize(n2);
    
    
    Bmat.reinit(6, _system.n_vars(), _elem.n_nodes()); // six stress-strain components
    
    std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>> mat_stiff
    (_property.get_property(MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX,
                            *this).release());
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        this->initialize_strain_operator(qp, Bmat);
        
        // get the material matrix
        (*mat_stiff)(xyz[qp], _system.time, material_mat);
        
        // calculate the stress
        Bmat.left_multiply(tmp_mat1_n1n2, material_mat);
        tmp_mat1_n1n2.vector_mult(tmp_vec1_n1, local_solution); // this is stress
        
        // now calculate the internal force vector
        Bmat.vector_mult_transpose(tmp_vec2_n2, tmp_vec1_n1);
        f.add(-JxW[qp], tmp_vec2_n2);
        
        // add the prestress
        
        if (request_jacobian) {
            
            Bmat.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat1_n1n2);
            jac.add(-JxW[qp], tmp_mat2_n2n2);
        }
    }
    
    return request_jacobian;
}



bool
MAST::StructuralElement3D::internal_force_sensitivity (bool request_jacobian,
                                                       DenseVector<Real>& f,
                                                       DenseMatrix<Real>& jac)
{
    // this should be true if the function is called
    libmesh_assert(this->sensitivity_param);
    libmesh_assert(!this->sensitivity_param->is_shape_parameter()); // this is not implemented for now
    
    
    // check if the material property or the provided exterior
    // values, like temperature, are functions of the sensitivity parameter
    bool calculate = false;
    calculate = calculate || _property.depends_on(*(this->sensitivity_param));
    
    // nothing to be calculated if the element does not depend on the
    // sensitivity parameter.
    if (!calculate)
        return false;
    
    FEMOperatorMatrix Bmat;
    
    const std::vector<Real>& JxW = _fe->get_JxW();
    const std::vector<Point>& xyz = _fe->get_xyz();
    const unsigned int n_phi = (unsigned int)JxW.size();
    const unsigned int n1=6, n2=6*n_phi;
    DenseMatrix<Real> material_mat, tmp_mat1_n1n2, tmp_mat2_n2n2;
    DenseVector<Real>  phi, tmp_vec1_n1, tmp_vec2_n2;
    
    tmp_mat1_n1n2.resize(n1, n2); tmp_mat2_n2n2.resize(n2, n2);
    phi.resize(n_phi); tmp_vec1_n1.resize(n1); tmp_vec2_n2.resize(n2);
    
    
    Bmat.reinit(6, _system.n_vars(), _elem.n_nodes()); // six stress-strain components

    std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>> mat_stiff
    (_property.get_property(MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX,
                            *this).release());

    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        this->initialize_strain_operator(qp, Bmat);
        
        // get the material matrix
        mat_stiff->total(*this->sensitivity_param,
                         xyz[qp], _system.time, material_mat);

        // calculate the stress
        Bmat.left_multiply(tmp_mat1_n1n2, material_mat);
        tmp_mat1_n1n2.vector_mult(tmp_vec1_n1, local_solution); // this is stress
        
        // now calculate the internal force vector
        Bmat.vector_mult_transpose(tmp_vec2_n2, tmp_vec1_n1);
        f.add(-JxW[qp], tmp_vec2_n2);
        
        // add the prestress
        
        if (request_jacobian) {
            
            Bmat.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat1_n1n2);
            jac.add(-JxW[qp], tmp_mat2_n2n2);
        }
    }
    
    return request_jacobian;
}




bool
MAST::StructuralElement3D::thermal_force (bool request_jacobian,
                                          DenseVector<Real>& f,
                                          DenseMatrix<Real>& jac,
                                          MAST::BoundaryCondition& p)
{
    FEMOperatorMatrix Bmat;
    
    const std::vector<Real>& JxW = _fe->get_JxW();
    const std::vector<Point>& xyz = _fe->get_xyz();
    const unsigned int n_phi = (unsigned int)_fe->get_phi().size();
    const unsigned int n1= 6, n2=6*n_phi;
    DenseMatrix<Real> material_exp_A_mat, material_exp_B_mat,
    tmp_mat1_n1n2, tmp_mat2_n2n2, stress;
    DenseVector<Real>  tmp_vec1_n1, tmp_vec2_n1, tmp_vec3_n2, delta_t;
    
    tmp_mat1_n1n2.resize(n1, n2); tmp_mat2_n2n2.resize(n2, n2);
    tmp_vec1_n1.resize(n1); tmp_vec2_n1.resize(n1);
    tmp_vec3_n2.resize(n2); delta_t.resize(1);
    
    
    Bmat.reinit(n1, _system.n_vars(), n_phi); // three stress-strain components
    
    std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>> mat
    (_property.get_property(MAST::SECTION_INTEGRATED_MATERIAL_THERMAL_EXPANSION_A_MATRIX,
                            *this).release());
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        // set the temperature vector to the value at this point
        //        _temperature->initialize(xyz[qp]);
        //        delta_t(0) = (*_temperature)() - _temperature->reference();
        
        // this is moved inside the domain since
        (*mat)(xyz[qp], _system.time, material_exp_A_mat);
        
        material_exp_A_mat.vector_mult(tmp_vec1_n1, delta_t); // [C]{alpha (T - T0)}
        
        this->initialize_strain_operator(qp, Bmat);
        
        Bmat.vector_mult_transpose(tmp_vec3_n2, tmp_vec1_n1);
        f.add(JxW[qp], tmp_vec3_n2);
    }
    
    // Jacobian contribution from von Karman strain
    return false;
}




bool
MAST::StructuralElement3D::thermal_force_sensitivity (bool request_jacobian,
                                                      DenseVector<Real>& f,
                                                      DenseMatrix<Real>& jac,
                                                      MAST::BoundaryCondition& p)
{
    
    libmesh_error(); // to be implemented
    
    return false;
}




void
MAST::StructuralElement3D::initialize_strain_operator(const unsigned int qp,
                                                      FEMOperatorMatrix& Bmat) {
    const std::vector<std::vector<RealVectorValue> >& dphi = _fe->get_dphi();
    
    unsigned int n_phi = (unsigned int)dphi.size();
    DenseVector<Real> phi; phi.resize(n_phi);
    
    // now set the shape function values
    // dN/dx
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
        phi(i_nd) = dphi[i_nd][qp](0);
    Bmat.set_shape_function(0, 0, phi); //  epsilon_xx = du/dx
    Bmat.set_shape_function(3, 1, phi); //  gamma_xy = dv/dx + ...
    Bmat.set_shape_function(5, 2, phi); //  gamma_zx = dw/dx + ...
    
    // dN/dy
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
        phi(i_nd) = dphi[i_nd][qp](1);
    Bmat.set_shape_function(1, 1, phi); //  epsilon_yy = dv/dy
    Bmat.set_shape_function(3, 0, phi); //  gamma_xy = du/dy + ...
    Bmat.set_shape_function(4, 2, phi); //  gamma_yz = dw/dy + ...
    
    // dN/dz
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
        phi(i_nd) = dphi[i_nd][qp](2);
    Bmat.set_shape_function(2, 2, phi); //  epsilon_xx = dw/dz
    Bmat.set_shape_function(4, 1, phi); //  gamma_xy = dv/dz + ...
    Bmat.set_shape_function(5, 0, phi); //  gamma_zx = du/dz + ...
    
}


