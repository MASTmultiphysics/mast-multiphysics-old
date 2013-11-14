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


// MAST includes
#include "StructuralElems/structural_elem_base.h"


// Forward declerations
class FEMOperatorMatrix;


namespace MAST {

    // Forward declerations
    class DKTElem;
    class BoundaryCondition;
    
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
        
        
        /*!
         *   @returns the unit normal vector out of the plane of this element.
         *   This is valid only for flat elements
         */
        const Point& normal() const {
            return _normal;
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
        
        /*!
         *   surface normal
         */
        Point _normal;
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

        
        /*!
         *    Calculates the sensitivity internal force vector and Jacobian due to
         *    strain energy
         */
        virtual bool internal_force_sensitivity(bool request_jacobian,
                                                DenseVector<Real>& f,
                                                DenseMatrix<Real>& jac);
        
        /*!
         *    Calculates the internal force vector and Jacobian due to
         *    strain energy coming from a prestress
         */
        virtual bool prestress_force_sensitivity (bool request_jacobian,
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
         *    Calculates the sensitivity fo force vector and Jacobian due
         *    to thermal stresses
         */
        virtual bool thermal_force_sensitivity(bool request_jacobian,
                                               DenseVector<Real>& f,
                                               DenseMatrix<Real>& jac);

        /*!
         *    Calculates the force vector and Jacobian due to surface pressure.
         */
        virtual bool surface_pressure_force(bool request_jacobian,
                                            DenseVector<Real>& f,
                                            DenseMatrix<Real>& jac,
                                            const unsigned int side,
                                            MAST::BoundaryCondition& p);

        /*!
         *    Calculates the force vector and Jacobian due to surface pressure.
         */
        virtual bool surface_pressure_force(bool request_jacobian,
                                            DenseVector<Real>& f,
                                            DenseMatrix<Real>& jac,
                                            MAST::BoundaryCondition& p);

        /*!
         *    Calculates the force vector and Jacobian due to surface pressure.
         */
        virtual bool surface_pressure_force_sensitivity(bool request_jacobian,
                                                        DenseVector<Real>& f,
                                                        DenseMatrix<Real>& jac,
                                                        const unsigned int side,
                                                        MAST::BoundaryCondition& p);
        
        /*!
         *    Calculates the force vector and Jacobian due to surface pressure.
         */
        virtual bool surface_pressure_force_sensitivity(bool request_jacobian,
                                                        DenseVector<Real>& f,
                                                        DenseMatrix<Real>& jac,
                                                        MAST::BoundaryCondition& p);


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
         *   performs integration at the quadrature point for the provided 
         *   matrices. The temperature vector and matrix entities are provided for
         *   integration
         */
        void _internal_force_operation(bool if_bending,
                                       bool if_dkt,
                                       bool if_vk,
                                       const unsigned int n2,
                                       const unsigned int qp,
                                       const std::vector<Real>& JxW,
                                       bool request_jacobian,
                                       DenseVector<Real>& local_f,
                                       DenseMatrix<Real>& local_jac,
                                       FEMOperatorMatrix& Bmat_mem,
                                       FEMOperatorMatrix& Bmat_bend,
                                       FEMOperatorMatrix& Bmat_trans,
                                       FEMOperatorMatrix& Bmat_vk,
                                       DenseMatrix<Real>& stress,
                                       DenseMatrix<Real>& vk_dwdxi_mat,
                                       DenseMatrix<Real>& material_A_mat,
                                       DenseMatrix<Real>& material_B_mat,
                                       DenseMatrix<Real>& material_D_mat,
                                       DenseMatrix<Real>& material_trans_shear_mat,
                                       DenseVector<Real>& tmp_vec1_n1,
                                       DenseVector<Real>& tmp_vec2_n1,
                                       DenseVector<Real>& tmp_vec3_n2,
                                       DenseVector<Real>& tmp_vec4_2,
                                       DenseVector<Real>& tmp_vec5_2,
                                       DenseMatrix<Real>& tmp_mat1_n1n2,
                                       DenseMatrix<Real>& tmp_mat2_n2n2,
                                       DenseMatrix<Real>& tmp_mat3,
                                       DenseMatrix<Real>& tmp_mat4_2n2);
        
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
}



#endif  // __MAST_structural_element_2D_h__
