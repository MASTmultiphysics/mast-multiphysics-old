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


// MAST includes
#include "StructuralElems/bending_structural_elem.h"


// Forward declerations
class FEMOperatorMatrix;


namespace MAST {
    
    // Forward declerations
    class BoundaryCondition;
    
    /*!
     *   class provides a simple mechanism to create a geometric element
     *   in the local coordinate system.
     */
    class Local1DElem : public MAST::LocalElemBase {
    public:
        /*!
         *    constructor takes a reference to the element for this which this 
         *    local element is to be created, and a vector in the element 
         *    x-y plane. The x-axis, by default, is along the element length.
         */
        Local1DElem(const Elem& elem,
                    const Point& y):
        MAST::LocalElemBase(elem),
        _local_y(y)
        {
            _create_local_elem();
            libmesh_assert_greater(_local_y.size(), 0.);
        }
        
        
        virtual ~Local1DElem() {
            // the local element may not have been created
            // for cases where the original element lies in the xy-plane
            if (_local_elem) {
                delete _local_elem;
                for (unsigned int i=0; i<_local_nodes.size(); i++)
                    delete _local_nodes[i];
            }
        }
        
        
    protected:
        
        /*!
         *    creation of an element in the local coordinate system
         */
        void _create_local_elem();
        
        /*!
         *    orientation of element local y-axis
         */
        Point _local_y;
    };
    
    
    
    class StructuralElement1D: public MAST::BendingStructuralElem {
        
    public:
        StructuralElement1D(System& sys,
                            const Elem& elem,
                            const MAST::ElementPropertyCardBase& p);
        
        /*!
         *    row dimension of the direct strain matrix, also used for the
         *    bending operator row dimension
         */
        virtual unsigned int n_direct_strain_components() {
            return 2;
        }
        
        /*!
         *    row dimension of the von Karman strain matrix
         */
        virtual unsigned int n_von_karman_strain_components() {
            return 2;
        }

        /*!
         *    Calculates the internal force vector and Jacobian due to
         *    strain energy
         */
        virtual bool internal_force(bool request_jacobian,
                                    DenseVector<Real>& f,
                                    DenseMatrix<Real>& jac,
                                    bool if_ignore_ho_jac);
        
        /*!
         *    Calculates the sensitivity internal force vector and Jacobian due to
         *    strain energy
         */
        virtual bool internal_force_sensitivity(bool request_jacobian,
                                                DenseVector<Real>& f,
                                                DenseMatrix<Real>& jac,
                                                bool if_ignore_ho_jac);
        
        /*!
         *    Calculates the internal force vector and Jacobian due to
         *    strain energy coming from a prestress
         */
        virtual bool prestress_force (bool request_jacobian,
                                      DenseVector<Real>& f,
                                      DenseMatrix<Real>& jac);
        
        /*!
         *    Calculates the internal force vector and Jacobian due to
         *    strain energy coming from a prestress
         */
        virtual bool prestress_force_sensitivity (bool request_jacobian,
                                                  DenseVector<Real>& f,
                                                  DenseMatrix<Real>& jac);

        /*!
         *   returns the value of maximum von Mises stress over the element
         */
        virtual Real max_von_mises_stress(){
            libmesh_error();
        }
        
        
        /*!
         *   returns the sensitivity of maximum von Mises stress over the element
         */
        virtual Real max_von_mises_stress_sensitivity(){
            libmesh_error();
        }

    protected:
        
        
        /*!
         *    Calculates the force vector and Jacobian due to thermal stresses
         */
        virtual bool thermal_force(bool request_jacobian,
                                   DenseVector<Real>& f,
                                   DenseMatrix<Real>& jac,
                                   MAST::BoundaryCondition& p);
        
        /*!
         *    Calculates the sensitivity of force vector and the Jacobian due to
         *    thermal stresses
         */
        virtual bool thermal_force_sensitivity(bool request_jacobian,
                                               DenseVector<Real>& f,
                                               DenseMatrix<Real>& jac,
                                               MAST::BoundaryCondition& p);

        /*!
         *   initialize membrane strain operator matrix
         */
        virtual void initialize_direct_strain_operator(const unsigned int qp,
                                                       FEMOperatorMatrix& Bmat);
        
        /*!
         *   initialze the von Karman strain in \par vK_strain, the operator
         *   matrices needed for Jacobian calculation.
         *   vk_strain = [dw/dx 0; 0 dw/dy; dw/dy dw/dx]
         *   Bmat_vk   = [dw/dx; dw/dy]
         */
        virtual void initialize_von_karman_strain_operator(const unsigned int qp,
                                                           DenseVector<Real>& vk_strain,
                                                           DenseMatrix<Real>& vk_dvdxi_mat,
                                                           DenseMatrix<Real>& vk_dwdxi_mat,
                                                           FEMOperatorMatrix& Bmat_v_vk,
                                                           FEMOperatorMatrix& Bmat_w_vk);
        
        /*!
         *   initialze the sensitivity of von Karman operator
         *   matrices needed for Jacobian calculation.
         *   vk_dwdxi_mat_sens = [dw/dx 0; 0 dw/dy; dw/dy dw/dx]
         */
        void initialize_von_karman_strain_operator_sensitivity
        (const unsigned int qp,
         DenseMatrix<Real>& vk_dvdxi_mat_sens,
         DenseMatrix<Real>& vk_dwdxi_mat_sens);
        

        /*!
         *   performs integration at the quadrature point for the provided
         *   matrices. The temperature vector and matrix entities are provided for
         *   integration
         */
        virtual void _internal_force_operation(bool if_bending,
                                               bool if_vk,
                                               const unsigned int n2,
                                               const unsigned int qp,
                                               const std::vector<Real>& JxW,
                                               bool request_jacobian,
                                               bool if_ignore_ho_jac,
                                               DenseVector<Real>& local_f,
                                               DenseMatrix<Real>& local_jac,
                                               FEMOperatorMatrix& Bmat_mem,
                                               FEMOperatorMatrix& Bmat_bend,
                                               FEMOperatorMatrix& Bmat_v_vk,
                                               FEMOperatorMatrix& Bmat_w_vk,
                                               DenseMatrix<Real>& stress,
                                               DenseMatrix<Real>& stress_l,
                                               DenseMatrix<Real>& vk_dvdxi_mat,
                                               DenseMatrix<Real>& vk_dwdxi_mat,
                                               DenseMatrix<Real>& material_A_mat,
                                               DenseMatrix<Real>& material_B_mat,
                                               DenseMatrix<Real>& material_D_mat,
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
         *   sensitivity of linear part of the geometric stiffness matrix
         */
        void _linearized_geometric_stiffness_sensitivity_with_static_solution
        (const unsigned int n2,
         const unsigned int qp,
         const std::vector<Real>& JxW,
         DenseMatrix<Real>& local_jac,
         FEMOperatorMatrix& Bmat_mem,
         FEMOperatorMatrix& Bmat_bend,
         FEMOperatorMatrix& Bmat_v_vk,
         FEMOperatorMatrix& Bmat_w_vk,
         DenseMatrix<Real>& stress_l,
         DenseMatrix<Real>& vk_dvdxi_mat,
         DenseMatrix<Real>& vk_dwdxi_mat,
         DenseMatrix<Real>& material_A_mat,
         DenseMatrix<Real>& material_B_mat,
         DenseVector<Real>& tmp_vec1_n1,
         DenseVector<Real>& tmp_vec2_n1,
         DenseMatrix<Real>& tmp_mat1_n1n2,
         DenseMatrix<Real>& tmp_mat2_n2n2,
         DenseMatrix<Real>& tmp_mat3);
        
    };
}



#endif // __MAST_structural_element_1d_h__
