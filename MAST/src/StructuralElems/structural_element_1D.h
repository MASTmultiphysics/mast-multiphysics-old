//
//  structural_element_1D.h
//  MAST
//
//  Created by Manav Bhatia on 10/21/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_structural_element_1D_h__
#define __MAST_structural_element_1D_h__

// C++ includes
#include <memory>


// libMesh includes
#include "libmesh/quadrature.h"


// MAST includes
#include "StructuralElems/structural_elem_base.h"
#include "PropertyCards/element_property_card_1D.h"
#include "Numerics/fem_operator_matrix.h"
#include "ThermalElems/temperature_function.h"
#include "StructuralElems/bernoulli_beam_elem.h"


namespace MAST {
    
    /*!
     *   class provides a simple mechanism to create a geometric element
     *   in the local coordinate system.
     */
    class Local1DELem {
    public:
        Local1DELem(const Elem& elem):
        _elem(elem),
        _local_elem(NULL)
        {
            this->create_local_elem();
        }
        
        
        /*!
         *   returns a constant reference to the global element.
         */
        const Elem& global_elem() const {
            return _elem;
        }
        
        
        /*!
         *   returns a constant reference to the local element.
         */
        const Elem& local_elem() const {
            return *_local_elem;
        }
        
    protected:
        
        void create_local_elem() {
            libmesh_error(); // to be implemented
        }
        
        /*!
         *   given element in global coordinate system
         */
        const Elem& _elem;
        
        /*!
         *   element created in local coordinate system
         */
        Elem* _local_elem;
        
    };
    
    
    class StructuralElement1D: public MAST::StructuralElementBase {
        
    public:
        StructuralElement1D(System& sys,
                            const Elem& elem,
                            const MAST::ElementPropertyCardBase& p);
        
        /*!
         *    Calculates the internal force vector and Jacobian due to
         *    strain energy
         */
        virtual bool internal_force(bool request_jacobian,
                                    DenseVector<Real>& f,
                                    DenseMatrix<Real>& jac);
        
        
    protected:
        
        /*!
         *    Calculates the force vector and Jacobian due to thermal stresses
         */
        virtual bool thermal_force(bool request_jacobian,
                                   DenseVector<Real>& f,
                                   DenseMatrix<Real>& jac);
        
        /*!
         *   initialize membrane strain operator matrix
         */
        void initialize_extension_strain_operator(const unsigned int qp,
                                                  FEMOperatorMatrix& Bmat);
        
        /*!
         *   initialze the bending strain operator for Mindling plate model.
         *   The bending part is included in \par Bend_par and the transverse
         *   shear part is included in \par Bmat_trans;
         */
        void initialize_timoshenko_bending_strain_operator(const unsigned int qp,
                                                           FEMOperatorMatrix& Bmat_bend,
                                                           FEMOperatorMatrix& Bmat_trans);
        
        
        /*!
         *   initialze the von Karman strain in \par vK_strain, the operator
         *   matrices needed for Jacobian calculation.
         *   vk_strain = [dw/dx 0; 0 dw/dy; dw/dy dw/dx]
         *   Bmat_vk   = [dw/dx; dw/dy]
         */
        void initialize_von_karman_strain_operator(const unsigned int qp,
                                                   DenseVector<Real>& vk_strain,
                                                   DenseMatrix<Real>& vk_dwdxi_mat,
                                                   FEMOperatorMatrix& Bmat_vk);
        
        /*!
         *    DKT element, initialized if the element needs it
         */
        std::auto_ptr<MAST::DKTElem> _dkt_elem;
        
        
        /*!
         *   element in local coordinate system
         */
        MAST::Local1DELem _local_elem;
    };
    
    
    
    MAST::StructuralElement1D::StructuralElement1D(System& sys,
                                                   const Elem& elem,
                                                   const MAST::ElementPropertyCardBase& p):
    _local_elem(elem),
    MAST::StructuralElementBase(sys, _local_elem.local_elem(), p)
    {
        // if the DKT element is being used, set it up now
        MAST::BendingModel1D bending_model =
        dynamic_cast<const MAST::ElementPropertyCard1D&>
        (_property).bending_model(elem, _fe->get_fe_type());
        
        if (bending_model == MAST::BERNOULLI)
            _dkt_elem.reset(new MAST::BernoulliElem(_local_elem.local_elem(),
                                                    *_qrule));
    }
    
    
    inline
    bool
    StructuralElement1D::internal_force (bool request_jacobian,
                                         DenseVector<Real>& f,
                                         DenseMatrix<Real>& jac)
    {
        FEMOperatorMatrix Bmat_mem, Bmat_bend, Bmat_trans, Bmat_vk;
        
        const std::vector<Real>& JxW = _fe->get_JxW();
        const std::vector<Point>& xyz = _fe->get_xyz();
        const unsigned int n_phi = (unsigned int)JxW.size();
        const unsigned int n1=3, n2=6*n_phi;
        DenseMatrix<Real> material_A_mat, material_B_mat, material_D_mat,
        material_trans_shear_mat, tmp_mat1_n1n2, tmp_mat2_n2n2, tmp_mat3,
        vk_dwdxi_mat, stress;
        DenseVector<Real>  phi, tmp_vec1_n1, tmp_vec2_n1, tmp_vec3_n2,
        tmp_vec4_2, tmp_vec5_2;
        
        tmp_mat1_n1n2.resize(n1, n2); tmp_mat2_n2n2.resize(n2, n2);
        vk_dwdxi_mat.resize(n1,2); stress.resize(2,2);
        phi.resize(n_phi); tmp_vec1_n1.resize(n1); tmp_vec2_n1.resize(n1);
        tmp_vec3_n2.resize(n2); tmp_vec4_2.resize(2); tmp_vec5_2.resize(2);
        
        
        Bmat_mem.reinit(3, _system.n_vars(), _elem.n_nodes()); // three stress-strain components
        
        const MAST::ElementPropertyCard1D& property_1D =
        dynamic_cast<const MAST::ElementPropertyCard1D&>(_property);
        bool if_vk = (property_1D.strain_type() == MAST::VON_KARMAN_STRAIN),
        if_bending = (property_1D.bending_model(_elem, _fe->get_fe_type()) != MAST::NO_BENDING_1D),
        if_dkt = (property_1D.bending_model(_elem, _fe->get_fe_type()) == MAST::DKT);
        
        for (unsigned int qp=0; qp<JxW.size(); qp++) {
            this->initialize_membrane_strain_operator(qp, Bmat_mem);
            
            // get the bending strain operator if needed
            tmp_vec2_n1.zero(); // used to store vk strain, if applicable
            if (if_bending) {
                if (if_dkt)
                    _dkt_elem->initialize_bending_strain_operator(qp, Bmat_bend);
                else
                    this->initialize_timoshenko_bending_strain_operator(qp,
                                                                        Bmat_bend,
                                                                        Bmat_trans);
                if (if_vk)  // get the vonKarman strain operator if needed
                    this->initialize_von_karman_strain_operator(qp,
                                                                tmp_vec2_n1,
                                                                vk_dwdxi_mat,
                                                                Bmat_vk);
            }
            
            // if temperature is specified, the initialize it to the current location
            if (_temperature)
                _temperature->initialize(xyz[qp]);
            
            // get the material matrix
            _property.calculate_matrix(_elem,
                                       MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX,
                                       material_A_mat);
            if (if_bending) {
                _property.calculate_matrix(_elem,
                                           MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_B_MATRIX,
                                           material_B_mat);
                _property.calculate_matrix(_elem,
                                           MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_D_MATRIX,
                                           material_D_mat);
                if (!if_dkt)
                    _property.calculate_matrix(_elem,
                                               MAST::SECTION_INTEGRATED_MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX,
                                               material_trans_shear_mat);
                
            }
            
            // first handle constant throught the thickness stresses: membrane and vonKarman
            Bmat_mem.vector_mult(tmp_vec1_n1, local_solution);
            tmp_vec2_n1.add(1., tmp_vec1_n1);  // epsilon_mem + epsilon_vk
            
            
            // multiply this with the constant through the thickness strain
            // membrane strain
            material_A_mat.vector_mult(tmp_vec1_n1, tmp_vec2_n1); // stress
            Bmat_mem.vector_mult_transpose(tmp_vec3_n2, tmp_vec1_n1);
            f.add(-JxW[qp], tmp_vec3_n2);
            // copy the stress values to a matrix if needed
            stress(0,0) = tmp_vec1_n1(0); // sigma_xx
            stress(0,1) = tmp_vec1_n1(2); // sigma_xy
            stress(1,0) = tmp_vec1_n1(2); // sigma_yx
            stress(1,1) = tmp_vec1_n1(1); // sigma_yy
            
            if (if_bending) {
                if (if_vk) {
                    // von Karman strain
                    vk_dwdxi_mat.vector_mult_transpose(tmp_vec4_2, tmp_vec1_n1);
                    Bmat_vk.vector_mult_transpose(tmp_vec3_n2, tmp_vec4_2);
                    f.add(-JxW[qp], tmp_vec3_n2);
                }
                
                // now coupling with the bending strain
                material_B_mat.vector_mult(tmp_vec1_n1, tmp_vec2_n1);
                Bmat_bend.vector_mult_transpose(tmp_vec3_n2, tmp_vec1_n1);
                f.add(-JxW[qp], tmp_vec3_n2);
                
                // now bending stress and its coupling
                Bmat_bend.vector_mult(tmp_vec2_n1, local_solution);
                
                // now get its projection onto the constant through thickness
                // and membrane operators
                material_B_mat.vector_mult(tmp_vec1_n1, tmp_vec2_n1);
                Bmat_mem.vector_mult_transpose(tmp_vec3_n2, tmp_vec1_n1);
                f.add(-JxW[qp], tmp_vec3_n2);
                
                if (if_vk) {
                    // von Karman strain
                    vk_dwdxi_mat.vector_mult_transpose(tmp_vec4_2, tmp_vec1_n1);
                    Bmat_vk.vector_mult_transpose(tmp_vec3_n2, tmp_vec4_2);
                    f.add(-JxW[qp], tmp_vec3_n2);
                }
                
                // and the bending operator
                material_D_mat.vector_mult(tmp_vec1_n1, tmp_vec2_n1);
                Bmat_bend.vector_mult_transpose(tmp_vec3_n2, tmp_vec1_n1);
                f.add(-JxW[qp], tmp_vec3_n2);
                
                if (!if_dkt) {
                    // now add the transverse shear component
                    Bmat_trans.vector_mult(tmp_vec4_2, local_solution);
                    material_trans_shear_mat.vector_mult(tmp_vec5_2, tmp_vec4_2);
                    Bmat_trans.vector_mult_transpose(tmp_vec3_n2, tmp_vec5_2);
                    f.add(-JxW[qp], tmp_vec3_n2);
                }
            }
            
            // add the prestress
            
            if (request_jacobian) {
                // membrane - membrane
                Bmat_mem.left_multiply(tmp_mat1_n1n2, material_A_mat);
                Bmat_mem.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat1_n1n2);
                jac.add(-JxW[qp], tmp_mat2_n2n2);
                
                if (if_bending) {
                    if (if_vk) {
                        // membrane - vk
                        tmp_mat3.resize(2, n2);
                        Bmat_vk.left_multiply(tmp_mat3, vk_dwdxi_mat);
                        tmp_mat3.left_multiply(material_A_mat);
                        Bmat_mem.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat3);
                        jac.add(-JxW[qp], tmp_mat2_n2n2);
                        
                        // vk - membrane
                        Bmat_mem.left_multiply(tmp_mat1_n1n2, material_A_mat);
                        vk_dwdxi_mat.get_transpose(tmp_mat3);
                        tmp_mat3.right_multiply(tmp_mat1_n1n2);
                        Bmat_vk.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat3);
                        jac.add(-JxW[qp], tmp_mat2_n2n2);
                        
                        // vk - vk
                        tmp_mat3.resize(2, n2);
                        Bmat_vk.left_multiply(tmp_mat3, stress);
                        Bmat_vk.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat3);
                        jac.add(-JxW[qp], tmp_mat2_n2n2);
                        
                        Bmat_vk.left_multiply(tmp_mat3, vk_dwdxi_mat);
                        tmp_mat3.left_multiply(material_A_mat);
                        tmp_mat3.left_multiply_transpose(vk_dwdxi_mat);
                        Bmat_vk.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat3);
                        jac.add(-JxW[qp], tmp_mat2_n2n2);
                        
                        // bending - vk
                        tmp_mat3.resize(2, n2);
                        Bmat_vk.left_multiply(tmp_mat3, vk_dwdxi_mat);
                        tmp_mat3.left_multiply(material_B_mat);
                        Bmat_bend.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat3);
                        jac.add(-JxW[qp], tmp_mat2_n2n2);
                        
                        // vk - bending
                        Bmat_bend.left_multiply(tmp_mat1_n1n2, material_B_mat);
                        vk_dwdxi_mat.get_transpose(tmp_mat3);
                        tmp_mat3.right_multiply(tmp_mat1_n1n2);
                        Bmat_vk.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat3);
                        jac.add(-JxW[qp], tmp_mat2_n2n2);
                    }
                    
                    // bending - membrane
                    Bmat_mem.left_multiply(tmp_mat1_n1n2, material_B_mat);
                    Bmat_bend.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat1_n1n2);
                    jac.add(-JxW[qp], tmp_mat2_n2n2);
                    
                    // membrane - bending
                    Bmat_bend.left_multiply(tmp_mat1_n1n2, material_B_mat);
                    Bmat_mem.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat1_n1n2);
                    jac.add(-JxW[qp], tmp_mat2_n2n2);
                    
                    // bending - bending
                    Bmat_bend.left_multiply(tmp_mat1_n1n2, material_D_mat);
                    Bmat_bend.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat1_n1n2);
                    jac.add(-JxW[qp], tmp_mat2_n2n2);
                    
                    if (!if_dkt) {
                        // now add the transverse shear component
                        Bmat_trans.left_multiply(tmp_mat1_n1n2, material_trans_shear_mat);
                        Bmat_trans.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat1_n1n2);
                        jac.add(-JxW[qp], tmp_mat2_n2n2);
                    }
                }
            }
        }
        
        return request_jacobian;
    }
    
    
    
    inline
    bool
    StructuralElement1D::thermal_force (bool request_jacobian,
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
        const std::vector<std::vector<RealVectorValue> >& dphi = _fe->get_dphi();
        
        Real temperature_value, ref_temperature;
        
        for (unsigned int qp=0; qp<JxW.size(); qp++) {
            this->initialize_membrane_strain_operator(qp, Bmat);
            
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
    StructuralElement1D::initialize_membrane_strain_operator(const unsigned int qp,
                                                             FEMOperatorMatrix& Bmat) {
        
        const std::vector<std::vector<RealVectorValue> >& dphi = _fe->get_dphi();
        
        unsigned int n_phi = (unsigned int)dphi.size();
        DenseVector<Real> phi; phi.resize(n_phi);
        
        // now set the shape function values
        // dN/dx
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi(i_nd) = dphi[i_nd][qp](0);
        Bmat.set_shape_function(0, 0, phi); //  epsilon_xx = du/dx
        Bmat.set_shape_function(2, 1, phi); //  gamma_xy = dv/dx + ...
        
        // dN/dy
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi(i_nd) = dphi[i_nd][qp](1);
        Bmat.set_shape_function(1, 1, phi); //  epsilon_yy = dv/dy
        Bmat.set_shape_function(2, 0, phi); //  gamma_xy = du/dy + ...
    }
    
    
    inline
    void
    StructuralElement1D::initialize_mindlin_bending_strain_operator(const unsigned int qp,
                                                                    FEMOperatorMatrix& Bmat_bend,
                                                                    FEMOperatorMatrix& Bmat_trans) {
        
        const std::vector<std::vector<RealVectorValue> >& dphi = _fe->get_dphi();
        const std::vector<std::vector<Real> >& phi = _fe->get_phi();
        
        const unsigned int n_phi = (unsigned int)phi.size();
        
        DenseVector<Real> phi_vec; phi_vec.resize(n_phi);
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = dphi[i_nd][qp](0);  // dphi/dx
        
        Bmat_bend.set_shape_function(0, 4, phi_vec); // epsilon-x: thetay
        Bmat_trans.set_shape_function(0, 2, phi_vec); // gamma-xz:  w
        phi_vec.scale(-1.0);
        Bmat_bend.set_shape_function(2, 3, phi_vec); // gamma-xy : thetax
        
        
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = dphi[i_nd][qp](1);  // dphi/dy
        
        Bmat_bend.set_shape_function(2, 4, phi_vec); // gamma-xy : thetay
        Bmat_trans.set_shape_function(1, 2, phi_vec); // gamma-yz : w
        phi_vec.scale(-1.0);
        Bmat_bend.set_shape_function(1, 3, phi_vec); // epsilon-y: thetax
        
        
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];  // phi
        
        Bmat_trans.set_shape_function(0, 4, phi_vec); // gamma-xz:  thetay
        phi_vec.scale(-1.0);
        Bmat_trans.set_shape_function(1, 3, phi_vec); // gamma-yz : thetax
    }
    
    
    inline
    void
    StructuralElement1D::initialize_von_karman_strain_operator(const unsigned int qp,
                                                               DenseVector<Real>& vk_strain,
                                                               DenseMatrix<Real>& vk_dwdxi_mat,
                                                               FEMOperatorMatrix& Bmat_vk) {
        
        const std::vector<std::vector<RealVectorValue> >& dphi = _fe->get_dphi();
        
        const unsigned int n_phi = (unsigned int)dphi.size();
        Real dw=0.;
        vk_dwdxi_mat.zero();
        
        DenseVector<Real> phi_vec; phi_vec.resize(n_phi);
        
        dw = 0.;
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ ) {
            phi_vec(i_nd) = dphi[i_nd][qp](0);  // dphi/dx
            dw += phi_vec(i_nd)*local_solution(2*n_phi+i_nd); // dw/dx
        }
        Bmat_vk.set_shape_function(0, 2, phi_vec); // dw/dx
        vk_dwdxi_mat(0, 0) = dw;  // epsilon-xx : dw/dx
        vk_dwdxi_mat(2, 1) = dw;  // gamma-xy : dw/dx
        
        dw = 0.;
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ ) {
            phi_vec(i_nd) = dphi[i_nd][qp](1);  // dphi/dy
            dw += phi_vec(i_nd)*local_solution(2*n_phi+i_nd); // dw/dy
        }
        Bmat_vk.set_shape_function(1, 2, phi_vec); // dw/dy
        vk_dwdxi_mat(1, 1) = dw;  // epsilon-yy : dw/dy
        vk_dwdxi_mat(2, 0) = dw;  // gamma-xy : dw/dy
    }
    
}


#endif // __MAST_structural_element_1d_h__
