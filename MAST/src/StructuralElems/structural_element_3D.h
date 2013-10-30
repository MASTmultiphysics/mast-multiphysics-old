//
//  structural_element_3d.h
//  MAST
//
//  Created by Manav Bhatia on 10/21/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_structural_element_3d_h__
#define __MAST_structural_element_3d_h__

// libMesh includes
#include "libmesh/quadrature.h"


// MAST includes
#include "StructuralElems/structural_elem_base.h"
#include "PropertyCards/element_property_card_base.h"
#include "Numerics/fem_operator_matrix.h"
#include "ThermalElems/temperature_function.h"


namespace MAST {
    
    class StructuralElement3D: public MAST::StructuralElementBase {
        
    public:
        StructuralElement3D(System& sys,
                            const Elem& elem,
                            const MAST::ElementPropertyCardBase& p):
        StructuralElementBase(sys, elem, p)
        {
            _init_fe_and_qrule(elem);
        }
        
        /*!
         *    Calculates the internal force vector and Jacobian due to
         *    strain energy
         */
        virtual bool internal_force(bool request_jacobian,
                                    DenseVector<Real>& f,
                                    DenseMatrix<Real>& jac);
        
        virtual bool prestress_force (bool request_jacobian,
                                      DenseVector<Real>& f,
                                      DenseMatrix<Real>& jac)
        { libmesh_error();}

    protected:
        
        /*!
         *    Calculates the force vector and Jacobian due to thermal stresses
         */
        virtual bool thermal_force(bool request_jacobian,
                                   DenseVector<Real>& f,
                                   DenseMatrix<Real>& jac);
        
        /*!
         *   initialize strain operator matrix
         */
        void initialize_strain_operator (const unsigned int qp,
                                         FEMOperatorMatrix& Bmat);
        
        /*!
         *   matrix that transforms the global dofs to the local element coordinate
         *   system
         */
        virtual const DenseMatrix<Real>& _transformation_matrix() const {
            libmesh_error(); // should not be called for a 3D elem
        }
    };
    
    
    
    inline
    bool
    StructuralElement3D::internal_force (bool request_jacobian,
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
        
        for (unsigned int qp=0; qp<JxW.size(); qp++) {
            this->initialize_strain_operator(qp, Bmat);
            
            // if temperature is specified, the initialize it to the current location
            if (_temperature)
                _temperature->initialize(xyz[qp]);
            
            // get the material matrix
            _property.calculate_matrix(_elem,
                                       MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX,
                                       material_mat);
            
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
    
    
    
    inline
    bool
    StructuralElement3D::thermal_force (bool request_jacobian,
                                        DenseVector<Real>& f,
                                        DenseMatrix<Real>& jac)
    {
        if (!_temperature) // only if a temperature load is specified
            return false;
        
        FEMOperatorMatrix Bmat;
        
        const std::vector<Real>& JxW = _fe->get_JxW();
        const std::vector<Point>& xyz = _fe->get_xyz();
        const unsigned int n_phi = (unsigned int)JxW.size();
        const unsigned int n1=6, n2=6*n_phi;
        DenseMatrix<Real> material_mat, expansion_mat;
        DenseVector<Real>  phi, temperature, tmp_vec1_n1, tmp_vec2_n1, tmp_vec3_n2;
        
        phi.resize(n_phi); tmp_vec1_n1.resize(n1); tmp_vec2_n1.resize(n1);
        tmp_vec3_n2.resize(n2); temperature.resize(6);
        
        Bmat.reinit(6, _system.n_vars(), _elem.n_nodes()); // six stress-strain components
        
        Real temperature_value, ref_temperature;
        
        for (unsigned int qp=0; qp<JxW.size(); qp++) {
            this->initialize_strain_operator(qp, Bmat);
            
            // set the temperature vector to the value at this point
            _temperature->initialize(xyz[qp]);
            temperature_value = _temperature->value();
            ref_temperature   = _temperature->reference();
            
            // this is moved inside the domain since
            _property.calculate_matrix(_elem,
                                       MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX,
                                       material_mat);
            _property.calculate_matrix(_elem,
                                       MAST::SECTION_INTEGRATED_MATERIAL_THERMAL_EXPANSION_MATRIX,
                                       expansion_mat);
            
            // calculate the strain
            expansion_mat.vector_mult(tmp_vec2_n1, tmp_vec1_n1);
            
            // calculate the stress
            material_mat.vector_mult(tmp_vec1_n1, tmp_vec2_n1);
            
            // now calculate the internal force vector
            Bmat.vector_mult_transpose(tmp_vec3_n2, tmp_vec1_n1);
            f.add(JxW[qp], tmp_vec3_n2);
        }
        
        return false;
    }
    
    
    
    inline
    void
    StructuralElement3D::initialize_strain_operator(const unsigned int qp,
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
    
    
}


#endif // __MAST_structural_element_3d_h__
