/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2014  Manav Bhatia
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

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
        Local1DElem(const libMesh::Elem& elem,
                    const libMesh::Point& y):
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
        libMesh::Point _local_y;
    };
    
    
    
    class StructuralElement1D: public MAST::BendingStructuralElem {
        
    public:
        StructuralElement1D(libMesh::System& sys,
                            const libMesh::Elem& elem,
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
                                    DenseRealVector& f,
                                    DenseRealMatrix& jac,
                                    bool if_ignore_ho_jac);
        
        /*!
         *    Calculates the sensitivity internal force vector and Jacobian due to
         *    strain energy
         */
        virtual bool internal_force_sensitivity(bool request_jacobian,
                                                DenseRealVector& f,
                                                DenseRealMatrix& jac,
                                                bool if_ignore_ho_jac);
        
        /*!
         *   calculates d[J]/d{x} . d{x}/dp
         */
        virtual bool
        internal_force_jac_dot_state_sensitivity (DenseRealMatrix& jac);

        /*!
         *    Calculates the internal force vector and Jacobian due to
         *    strain energy coming from a prestress
         */
        virtual bool prestress_force (bool request_jacobian,
                                      DenseRealVector& f,
                                      DenseRealMatrix& jac);
        
        /*!
         *    Calculates the internal force vector and Jacobian due to
         *    strain energy coming from a prestress
         */
        virtual bool prestress_force_sensitivity (bool request_jacobian,
                                                  DenseRealVector& f,
                                                  DenseRealMatrix& jac);

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
                                   DenseRealVector& f,
                                   DenseRealMatrix& jac,
                                   MAST::BoundaryCondition& p);
        
        /*!
         *    Calculates the sensitivity of force vector and the Jacobian due to
         *    thermal stresses
         */
        virtual bool thermal_force_sensitivity(bool request_jacobian,
                                               DenseRealVector& f,
                                               DenseRealMatrix& jac,
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
                                                           DenseRealVector& vk_strain,
                                                           DenseRealMatrix& vk_dvdxi_mat,
                                                           DenseRealMatrix& vk_dwdxi_mat,
                                                           FEMOperatorMatrix& Bmat_v_vk,
                                                           FEMOperatorMatrix& Bmat_w_vk);
        
        /*!
         *   initialze the sensitivity of von Karman operator
         *   matrices needed for Jacobian calculation.
         *   vk_dwdxi_mat_sens = [dw/dx 0; 0 dw/dy; dw/dy dw/dx]
         */
        void initialize_von_karman_strain_operator_sensitivity
        (const unsigned int qp,
         DenseRealMatrix& vk_dvdxi_mat_sens,
         DenseRealMatrix& vk_dwdxi_mat_sens);
        

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
                                               DenseRealVector& local_f,
                                               DenseRealMatrix& local_jac,
                                               FEMOperatorMatrix& Bmat_mem,
                                               FEMOperatorMatrix& Bmat_bend,
                                               FEMOperatorMatrix& Bmat_v_vk,
                                               FEMOperatorMatrix& Bmat_w_vk,
                                               DenseRealMatrix& stress,
                                               DenseRealMatrix& stress_l,
                                               DenseRealMatrix& vk_dvdxi_mat,
                                               DenseRealMatrix& vk_dwdxi_mat,
                                               DenseRealMatrix& material_A_mat,
                                               DenseRealMatrix& material_B_mat,
                                               DenseRealMatrix& material_D_mat,
                                               DenseRealVector& vec1_n1,
                                               DenseRealVector& vec2_n1,
                                               DenseRealVector& vec3_n2,
                                               DenseRealVector& vec4_2,
                                               DenseRealVector& vec5_2,
                                               DenseRealMatrix& mat1_n1n2,
                                               DenseRealMatrix& mat2_n2n2,
                                               DenseRealMatrix& mat3,
                                               DenseRealMatrix& mat4_2n2);
        
        
        /*!
         *   sensitivity of linear part of the geometric stiffness matrix
         */
        void _linearized_geometric_stiffness_sensitivity_with_static_solution
        (const unsigned int n2,
         const unsigned int qp,
         const std::vector<Real>& JxW,
         DenseRealMatrix& local_jac,
         FEMOperatorMatrix& Bmat_mem,
         FEMOperatorMatrix& Bmat_bend,
         FEMOperatorMatrix& Bmat_v_vk,
         FEMOperatorMatrix& Bmat_w_vk,
         DenseRealMatrix& stress_l,
         DenseRealMatrix& vk_dvdxi_mat,
         DenseRealMatrix& vk_dwdxi_mat,
         DenseRealMatrix& material_A_mat,
         DenseRealMatrix& material_B_mat,
         DenseRealVector& vec1_n1,
         DenseRealVector& vec2_n1,
         DenseRealMatrix& mat1_n1n2,
         DenseRealMatrix& mat2_n2n2,
         DenseRealMatrix& mat3);
        
    };
}



#endif // __MAST_structural_element_1d_h__
