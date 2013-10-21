//
//  structural_element_3d.h
//  MAST
//
//  Created by Manav Bhatia on 10/21/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_structural_element_3d_h__
#define __MAST_structural_element_3d_h__

// MAST includes
#include "StructuralElems/structural_elem_base.h"
#include "PropertyCards/element_property_card_base.h"
#include "Numerics/fem_operator_matrix.h"

// libMesh includes
#include "libmesh/quadrature.h"



namespace MAST {
    
    class StructuralElement3D: public MAST::StructuralElementBase {
        
    public:
        StructuralElement3D():
        StructuralElementBase()
        { }
        
        virtual bool internal_force(bool request_jacobian,
                                    DenseVector<Real>& f,
                                    DenseMatrix<Real>& jac);
        
    };
}



bool
MAST::StructuralElementBase::internal_force (bool request_jacobian,
                                             DenseVector<Real>& f,
                                             DenseMatrix<Real>& jac)
{
    FEMOperatorMatrix Bmat;
    FEBase* fe = NULL;
    QBase* qrule = NULL;
    get_fe_and_qrule(fe, qrule);
    
    const std::vector<Real>& JxW = fe->get_JxW();
    DenseMatrix<Real> material_mat, tmp_mat1_n1n2, tmp_mat2_n2n2;
    DenseVector<Real>  phi, tmp_vec1_n1, tmp_vec2_n2;
    
    _property->calculate_matrix(*_elem,
                                MAST::MATERIAL_STIFFNESS_MATRIX,
                                material_mat);

    Bmat.reinit(6, _system->n_vars()); // six stress-strain components
    phi.resize(JxW.size());
    const std::vector<std::vector<RealVectorValue> >& dphi = fe->get_dphi();
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {

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

        
        // calculate the stress
        Bmat->left_multiply(tmp_mat1_n1n2, material_mat);
        tmp_mat1_n1n2.vector_mult(tmp_vec1_n1, _local_solution); // this is stress

        // now calculate the internal force vector
        Bmat->vector_mult_transpose(tmp_vec2_n2, tmp_vec1_n1);
        f.add(JxW[qp], tmp_vec2_n2);

        // add the prestress

        if (request_jacobian) {
            
            Bmat->right_multiply_transpose(tmp_mat2_n2n2, tmp_mat1_n1n2);
            jac.add(JxW[qp], tmp_mat2_n2n2);
        }
    }
    
    return request_jacobian;
}


#endif // __MAST_structural_element_3d_h__
