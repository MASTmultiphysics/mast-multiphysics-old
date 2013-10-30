//
//  structural_element_2D.h
//  MAST
//
//  Created by Manav Bhatia on 10/21/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_structural_element_2D_h__
#define __MAST_structural_element_2D_h__

// C++ includes
#include <memory>


// libMesh includes
#include "libmesh/quadrature.h"


// MAST includes
#include "StructuralElems/structural_elem_base.h"
#include "PropertyCards/element_property_card_2D.h"
#include "Numerics/fem_operator_matrix.h"
#include "ThermalElems/temperature_function.h"
#include "StructuralElems/dkt_elem.h"


namespace MAST {
    
    /*!
     *   class provides a simple mechanism to create a geometric element
     *   in the local coordinate system.
     */
    class Local2DElem {
    public:
        Local2DElem(const Elem& elem):
        _elem(elem),
        _local_elem(NULL)
        {
            _create_local_elem();
        }
        
        
        ~Local2DElem() {
            delete _local_elem;
            for (unsigned int i=0; i<_local_nodes.size(); i++)
                delete _local_nodes[i];
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

        /*!
         *    returns the transformation matrix for this element. This is used 
         *    to map the coordinates from local to global coordinate system
         */
        const DenseMatrix<Real>& T_matrix() const {
            return _T_mat;
        }
        
    protected:
        
        /*!
         *    creation of an element in the local coordinate system
         */
        void _create_local_elem();
        
        /*!
         *   given element in global coordinate system
         */
        const Elem& _elem;
        
        /*!
         *   element created in local coordinate system
         */
        Elem* _local_elem;
        
        /*!
         *   nodes for local element
         */
        std::vector<Node*> _local_nodes;
        
        /*!
         *    Transformation matrix defines T_ij = V_i^t . Vn_j, where
         *    V_i are the unit vectors of the global cs, and Vn_j are the 
         *    unit vectors of the local cs. To transform a vector from global to 
         *    local cs,    an_j = T^t a_i, and the reverse transformation is 
         *    obtained as  a_j  = T  an_i
         */
        DenseMatrix<Real> _T_mat;
        
    };
    
    
    class StructuralElement2D: public MAST::StructuralElementBase {
        
    public:
        StructuralElement2D(System& sys,
                            const Elem& elem,
                            const MAST::ElementPropertyCardBase& p);
        
        /*!
         *    Calculates the internal force vector and Jacobian due to
         *    strain energy
         */
        virtual bool internal_force(bool request_jacobian,
                                    DenseVector<Real>& f,
                                    DenseMatrix<Real>& jac);
        
        /*!
         *    Calculates the internal force vector and Jacobian due to
         *    strain energy coming from a prestress
         */
        virtual bool prestress_force (bool request_jacobian,
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
        void initialize_membrane_strain_operator(const unsigned int qp,
                                                 FEMOperatorMatrix& Bmat);
        
        /*!
         *   initialze the bending strain operator for Mindling plate model.
         *   The bending part is included in \par Bend_par and the transverse
         *   shear part is included in \par Bmat_trans;
         */
        void initialize_mindlin_bending_strain_operator(const unsigned int qp,
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
         *   matrix that transforms the global dofs to the local element coordinate
         *   system
         */
        virtual const DenseMatrix<Real>& _transformation_matrix() const {
            return _local_elem.T_matrix();
        }
            
        /*!
         *    DKT element, initialized if the element needs it
         */
        std::auto_ptr<MAST::DKTElem> _dkt_elem;
        
        
        /*!
         *   element in local coordinate system
         */
        MAST::Local2DElem _local_elem;
        
    };
    
    
    
    inline
    void
    MAST::Local2DElem::_create_local_elem() {
        
        libmesh_assert(_elem.dim() == 2);
        _local_elem = Elem::build(_elem.type()).release();
        _local_nodes.resize(_elem.n_nodes());
        for (unsigned int i=0; i<_elem.n_nodes(); i++) {
            _local_nodes[i] = new Node;
            _local_nodes[i]->set_id() = _elem.get_node(i)->id();
            _local_elem->set_node(i) = _local_nodes[i];
        }
        
        // first node is the origin of the new cs
        // calculate the coordinate system for the plane of the element
        Point v1, v2, v3, p;
        v1 = *_elem.get_node(1); v1 -= *_elem.get_node(0); v1 /= v1.size(); // local x
        v2 = *_elem.get_node(2); v2 -= *_elem.get_node(0); v2 /= v2.size();
        v3 = v1.cross(v2); v3 /= v3.size();      // local z
        v2 = v3.cross(v1); v2 /= v2.size();      // local y
        
        // now the transformation matrix from old to new cs
        //        an_i vn_i = a_i v_i
        //        an_j = a_i v_i.vn_j  = a_i t_ij = T^t a_i
        //        t_ij = v_i.vn_j

        _T_mat.resize(3,3);
        for (unsigned int i=0; i<3; i++) {
            _T_mat(i,0) = v1(i);
            _T_mat(i,1) = v2(i);
            _T_mat(i,2) = v3(i);
        }
        
        // now calculate the new coordinates with respect to the origin
        for (unsigned int i=0; i<_local_nodes.size(); i++) {
            p = *_elem.get_node(i);
            p -= *_elem.get_node(0); // local wrt origin
            for (unsigned int j=0; j<3; j++)
                for (unsigned int k=0; k<3; k++)
                    (*_local_nodes[i])(j) += _T_mat(k,j)*p(k);
        }
    }

    
    
    MAST::StructuralElement2D::StructuralElement2D(System& sys,
                                                   const Elem& elem,
                                                   const MAST::ElementPropertyCardBase& p):
    MAST::StructuralElementBase(sys, elem, p),
    _local_elem(elem)
    {
        _init_fe_and_qrule(_local_elem.local_elem());
        
        // if the DKT element is being used, set it up now
        MAST::BendingModel2D bending_model =
        dynamic_cast<const MAST::ElementPropertyCard2D&>
        (_property).bending_model(elem, _fe->get_fe_type());
        
        if (bending_model == MAST::DKT)
            _dkt_elem.reset(new MAST::DKTElem(_local_elem.local_elem(),
                                              *_qrule));
    }
    
    
    inline
    bool
    StructuralElement2D::internal_force (bool request_jacobian,
                                         DenseVector<Real>& f,
                                         DenseMatrix<Real>& jac)
    {
        FEMOperatorMatrix Bmat_mem, Bmat_bend, Bmat_trans, Bmat_vk;
        
        const std::vector<Real>& JxW = _fe->get_JxW();
        const std::vector<Point>& xyz = _fe->get_xyz();
        const unsigned int n_phi = (unsigned int)_fe->get_phi().size();
        const unsigned int n1=3, n2=6*n_phi;
        DenseMatrix<Real> material_A_mat, material_B_mat, material_D_mat,
        material_trans_shear_mat, tmp_mat1_n1n2, tmp_mat2_n2n2, tmp_mat3,
        tmp_mat4_2n2, vk_dwdxi_mat, stress, local_jac;
        DenseVector<Real>  tmp_vec1_n1, tmp_vec2_n1, tmp_vec3_n2,
        tmp_vec4_2, tmp_vec5_2, local_f;
        
        tmp_mat1_n1n2.resize(n1, n2); tmp_mat2_n2n2.resize(n2, n2);
        tmp_mat4_2n2.resize(2, n2); local_jac.resize(n2, n2);
        vk_dwdxi_mat.resize(n1,2); stress.resize(2,2); local_f.resize(n2);
        tmp_vec1_n1.resize(n1); tmp_vec2_n1.resize(n1);
        tmp_vec3_n2.resize(n2); tmp_vec4_2.resize(2); tmp_vec5_2.resize(2);
        
        
        Bmat_mem.reinit(3, _system.n_vars(), n_phi); // three stress-strain components
        Bmat_bend.reinit(3, _system.n_vars(), n_phi);
        Bmat_trans.reinit(2, _system.n_vars(), n_phi); // only two shear stresses
        Bmat_vk.reinit(2, _system.n_vars(), n_phi); // only dw/dx and dw/dy
        
        const MAST::ElementPropertyCard2D& property_2D =
        dynamic_cast<const MAST::ElementPropertyCard2D&>(_property);
        bool if_vk = (property_2D.strain_type() == MAST::VON_KARMAN_STRAIN),
        if_bending = (property_2D.bending_model(_elem, _fe->get_fe_type()) != MAST::NO_BENDING_2D),
        if_dkt = (property_2D.bending_model(_elem, _fe->get_fe_type()) == MAST::DKT);
        
        for (unsigned int qp=0; qp<JxW.size(); qp++) {
            this->initialize_membrane_strain_operator(qp, Bmat_mem);
            
            // get the bending strain operator if needed
            tmp_vec2_n1.zero(); // used to store vk strain, if applicable
            if (if_bending) {
                if (if_dkt)
                    _dkt_elem->initialize_dkt_bending_strain_operator(qp, Bmat_bend);
                else
                    this->initialize_mindlin_bending_strain_operator(qp,
                                                                     Bmat_bend,
                                                                     Bmat_trans);
                if (if_vk)  // get the vonKarman strain operator if needed
                    this->initialize_von_karman_strain_operator(qp,
                                                                tmp_vec2_n1, // epsilon_vk
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
            local_f.add(-JxW[qp], tmp_vec3_n2);
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
                    local_f.add(-JxW[qp], tmp_vec3_n2);
                }
                
                // now coupling with the bending strain
                material_B_mat.vector_mult(tmp_vec1_n1, tmp_vec2_n1);
                Bmat_bend.vector_mult_transpose(tmp_vec3_n2, tmp_vec1_n1);
                local_f.add(-JxW[qp], tmp_vec3_n2);
                
                // now bending stress and its coupling
                Bmat_bend.vector_mult(tmp_vec2_n1, local_solution);
                
                // now get its projection onto the constant through thickness
                // and membrane operators
                material_B_mat.vector_mult(tmp_vec1_n1, tmp_vec2_n1);
                Bmat_mem.vector_mult_transpose(tmp_vec3_n2, tmp_vec1_n1);
                local_f.add(-JxW[qp], tmp_vec3_n2);
                
                if (if_vk) {
                    // von Karman strain
                    vk_dwdxi_mat.vector_mult_transpose(tmp_vec4_2, tmp_vec1_n1);
                    Bmat_vk.vector_mult_transpose(tmp_vec3_n2, tmp_vec4_2);
                    local_f.add(-JxW[qp], tmp_vec3_n2);
                }
                
                // and the bending operator
                material_D_mat.vector_mult(tmp_vec1_n1, tmp_vec2_n1);
                Bmat_bend.vector_mult_transpose(tmp_vec3_n2, tmp_vec1_n1);
                local_f.add(-JxW[qp], tmp_vec3_n2);
                
                if (!if_dkt) {
                    // now add the transverse shear component
                    Bmat_trans.vector_mult(tmp_vec4_2, local_solution);
                    material_trans_shear_mat.vector_mult(tmp_vec5_2, tmp_vec4_2);
                    Bmat_trans.vector_mult_transpose(tmp_vec3_n2, tmp_vec5_2);
                    local_f.add(-JxW[qp], tmp_vec3_n2);
                }
            }
            
            if (request_jacobian) {
                // membrane - membrane
                Bmat_mem.left_multiply(tmp_mat1_n1n2, material_A_mat);
                Bmat_mem.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat1_n1n2);
                local_jac.add(-JxW[qp], tmp_mat2_n2n2);
                
                if (if_bending) {
                    if (if_vk) {
                        // membrane - vk
                        tmp_mat3.resize(2, n2);
                        Bmat_vk.left_multiply(tmp_mat3, vk_dwdxi_mat);
                        tmp_mat3.left_multiply(material_A_mat);
                        Bmat_mem.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat3);
                        local_jac.add(-JxW[qp], tmp_mat2_n2n2);
                        
                        // vk - membrane
                        Bmat_mem.left_multiply(tmp_mat1_n1n2, material_A_mat);
                        vk_dwdxi_mat.get_transpose(tmp_mat3);
                        tmp_mat3.right_multiply(tmp_mat1_n1n2);
                        Bmat_vk.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat3);
                        local_jac.add(-JxW[qp], tmp_mat2_n2n2);
                        
                        // vk - vk
                        tmp_mat3.resize(2, n2);
                        Bmat_vk.left_multiply(tmp_mat3, stress);
                        Bmat_vk.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat3);
                        local_jac.add(-JxW[qp], tmp_mat2_n2n2);
                        
                        Bmat_vk.left_multiply(tmp_mat3, vk_dwdxi_mat);
                        tmp_mat3.left_multiply(material_A_mat);
                        tmp_mat3.left_multiply_transpose(vk_dwdxi_mat);
                        Bmat_vk.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat3);
                        local_jac.add(-JxW[qp], tmp_mat2_n2n2);
                        
                        // bending - vk
                        tmp_mat3.resize(2, n2);
                        Bmat_vk.left_multiply(tmp_mat3, vk_dwdxi_mat);
                        tmp_mat3.left_multiply(material_B_mat);
                        Bmat_bend.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat3);
                        local_jac.add(-JxW[qp], tmp_mat2_n2n2);
                        
                        // vk - bending
                        Bmat_bend.left_multiply(tmp_mat1_n1n2, material_B_mat);
                        vk_dwdxi_mat.get_transpose(tmp_mat3);
                        tmp_mat3.right_multiply(tmp_mat1_n1n2);
                        Bmat_vk.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat3);
                        local_jac.add(-JxW[qp], tmp_mat2_n2n2);
                    }
                    
                    // bending - membrane
                    Bmat_mem.left_multiply(tmp_mat1_n1n2, material_B_mat);
                    Bmat_bend.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat1_n1n2);
                    local_jac.add(-JxW[qp], tmp_mat2_n2n2);
                    
                    // membrane - bending
                    Bmat_bend.left_multiply(tmp_mat1_n1n2, material_B_mat);
                    Bmat_mem.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat1_n1n2);
                    local_jac.add(-JxW[qp], tmp_mat2_n2n2);
                    
                    // bending - bending
                    Bmat_bend.left_multiply(tmp_mat1_n1n2, material_D_mat);
                    Bmat_bend.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat1_n1n2);
                    local_jac.add(-JxW[qp], tmp_mat2_n2n2);
                    
                    if (!if_dkt) {
                        // now add the transverse shear component
                        Bmat_trans.left_multiply(tmp_mat4_2n2, material_trans_shear_mat);
                        Bmat_trans.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat4_2n2);
                        local_jac.add(-JxW[qp], tmp_mat2_n2n2);
                    }
                }
            }
        }
        
        // now transform to the global coorodinate system
        _transform_to_global_system(local_f, tmp_vec3_n2);
        f.add(1., tmp_vec3_n2);
        if (request_jacobian) {
            _transform_to_global_system(local_jac, tmp_mat2_n2n2);
            jac.add(1., tmp_mat2_n2n2);
        }
        
        return request_jacobian;
    }

    
    
    inline
    bool
    StructuralElement2D::prestress_force (bool request_jacobian,
                                          DenseVector<Real>& f,
                                          DenseMatrix<Real>& jac)
    {
        if (!_property.if_prestressed())
            return false;
        
        FEMOperatorMatrix Bmat_mem, Bmat_bend, Bmat_trans, Bmat_vk;
        
        const std::vector<Real>& JxW = _fe->get_JxW();
        const unsigned int n_phi = (unsigned int)_fe->get_phi().size();
        const unsigned int n1=3, n2=6*n_phi;
        DenseMatrix<Real> tmp_mat2_n2n2, tmp_mat3, vk_dwdxi_mat, local_jac,
        prestress_mat_A;
        DenseVector<Real> tmp_vec2_n1, tmp_vec3_n2, tmp_vec4_2, tmp_vec5_2,
        local_f, prestress_vec_A, prestress_vec_B;
        
        tmp_mat2_n2n2.resize(n2, n2); local_jac.resize(n2, n2);
        vk_dwdxi_mat.resize(n1,2); local_f.resize(n2); tmp_vec2_n1.resize(n1);
        tmp_vec3_n2.resize(n2); tmp_vec4_2.resize(2); tmp_vec5_2.resize(2);
        
        
        Bmat_mem.reinit(3, _system.n_vars(), n_phi); // three stress-strain components
        Bmat_bend.reinit(3, _system.n_vars(), n_phi);
        Bmat_trans.reinit(2, _system.n_vars(), n_phi); // only two shear stresses
        Bmat_vk.reinit(2, _system.n_vars(), n_phi); // only dw/dx and dw/dy
        
        const MAST::ElementPropertyCard2D& property_2D =
        dynamic_cast<const MAST::ElementPropertyCard2D&>(_property);
        bool if_vk = (property_2D.strain_type() == MAST::VON_KARMAN_STRAIN),
        if_bending = (property_2D.bending_model(_elem, _fe->get_fe_type()) != MAST::NO_BENDING_2D),
        if_dkt = (property_2D.bending_model(_elem, _fe->get_fe_type()) == MAST::DKT);
        
        
        // get the element prestress
        _property.prestress_matrix(SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX,
                                   prestress_mat_A);
        _property.prestress_vector(SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX,
                                   prestress_vec_A);

        _property.prestress_vector(SECTION_INTEGRATED_MATERIAL_STIFFNESS_B_MATRIX,
                                   prestress_vec_B);
        // transform to the local coordinate system
        

        for (unsigned int qp=0; qp<JxW.size(); qp++) {
            this->initialize_membrane_strain_operator(qp, Bmat_mem);
            
            // get the bending strain operator if needed
            tmp_vec2_n1.zero(); // used to store vk strain, if applicable
            if (if_bending) {
                if (if_dkt)
                    _dkt_elem->initialize_dkt_bending_strain_operator(qp, Bmat_bend);
                else
                    this->initialize_mindlin_bending_strain_operator(qp,
                                                                     Bmat_bend,
                                                                     Bmat_trans);
                if (if_vk)  // get the vonKarman strain operator if needed
                    this->initialize_von_karman_strain_operator(qp,
                                                                tmp_vec2_n1,
                                                                vk_dwdxi_mat,
                                                                Bmat_vk);
            }
            
            // first handle constant throught the thickness stresses: membrane and vonKarman
            // multiply this with the constant through the thickness strain
            // membrane strain
            Bmat_mem.vector_mult_transpose(tmp_vec3_n2, prestress_vec_A);
            local_f.add(-JxW[qp], tmp_vec3_n2); // epsilon_mem * sigma_0
            
            if (if_bending) {
                if (if_vk) {
                    // von Karman strain
                    vk_dwdxi_mat.vector_mult_transpose(tmp_vec4_2, prestress_vec_A);
                    Bmat_vk.vector_mult_transpose(tmp_vec3_n2, tmp_vec4_2);
                    local_f.add(-JxW[qp], tmp_vec3_n2); // epsilon_vk * sigma_0
                }
                
                // now coupling with the bending strain
                Bmat_bend.vector_mult_transpose(tmp_vec3_n2, prestress_vec_B);
                local_f.add(-JxW[qp], tmp_vec3_n2); // epsilon_bend * sigma_0
            }
            
            if (request_jacobian) {
                if (if_bending && if_vk) {
                    tmp_mat3.resize(2, n2);
                    Bmat_vk.left_multiply(tmp_mat3, prestress_mat_A);
                    Bmat_vk.right_multiply_transpose(tmp_mat2_n2n2, tmp_mat3);
                    local_jac.add(-JxW[qp], tmp_mat2_n2n2);
                }
            }
        }
        
        // now transform to the global coorodinate system
        _transform_to_global_system(local_f, tmp_vec3_n2);
        f.add(1., tmp_vec3_n2);
        if (request_jacobian && if_vk) {
            _transform_to_global_system(local_jac, tmp_mat2_n2n2);
            jac.add(1., tmp_mat2_n2n2);
        }
        
        // only the nonlinear strain returns a Jacobian for prestressing
        return (request_jacobian && if_vk);
    }

    
    
    
    inline
    bool
    StructuralElement2D::thermal_force (bool request_jacobian,
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
        DenseVector<Real>  phi, temperature, tmp_vec1_n1, tmp_vec2_n1,
        tmp_vec3_n2, local_f;
        
        phi.resize(n_phi); tmp_vec1_n1.resize(n1); tmp_vec2_n1.resize(n1);
        tmp_vec3_n2.resize(n2); temperature.resize(6); local_f.resize(n2);
        
        Bmat.reinit(3, _system.n_vars(), n_phi); // three stress-strain components
        
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
            local_f.add(JxW[qp], tmp_vec3_n2);
        }

        // now transform to the global system and add
        _transform_to_global_system(local_f, tmp_vec3_n2);
        f.add(1., tmp_vec3_n2);
        
        return false;
    }
    
    
    
    inline
    void
    StructuralElement2D::initialize_membrane_strain_operator(const unsigned int qp,
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
    StructuralElement2D::initialize_mindlin_bending_strain_operator(const unsigned int qp,
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
    StructuralElement2D::initialize_von_karman_strain_operator(const unsigned int qp,
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



#endif  // __MAST_structural_element_2D_h__
