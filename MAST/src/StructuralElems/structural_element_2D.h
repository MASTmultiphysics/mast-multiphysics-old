//
//  structural_element_2D.h
//  MAST
//
//  Created by Manav Bhatia on 10/21/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_structural_element_2D_h__
#define __MAST_structural_element_2D_h__


// libMesh includes
#include "libmesh/quadrature.h"


// MAST includes
#include "StructuralElems/structural_elem_base.h"
#include "PropertyCards/element_property_card_2D.h"
#include "Numerics/fem_operator_matrix.h"
#include "ThermalElems/temperature_function.h"


namespace MAST {

    /*!
     *   class provides a simple mechanism to create a geometric element 
     *   in the local coordinate system.
     */
    class Local2DELem {
    public:
        Local2DELem(const Elem& elem):
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
    
    
    class StructuralElement2D: public MAST::StructuralElementBase {
        
    public:
        StructuralElement2D(System& sys,
                            const Elem& elem,
                            const MAST::ElementPropertyCardBase& p):
        _local_elem(elem),
        StructuralElementBase(sys, _local_elem.local_elem(), p)
        { }
        
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
         *   initialize strain operator matrix
         */
        void initialize_membrane_strain_operator
        (const std::vector<std::vector<RealVectorValue> >& dphi,
         const unsigned int qp,
         FEMOperatorMatrix& Bmat);
        
        /*!
         *   element in local coordinate system
         */
        MAST::Local2DELem _local_elem;
    };
    
    
    
    inline
    bool
    StructuralElement2D::internal_force (bool request_jacobian,
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
        const std::vector<std::vector<RealVectorValue> >& dphi = _fe->get_dphi();
        
        for (unsigned int qp=0; qp<JxW.size(); qp++) {
            this->initialize_membrane_strain_operator(dphi, qp, Bmat);
            
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
        DenseVector<Real>  phi, temperature, tmp_vec1_n1, tmp_vec2_n1, tmp_vec3_n2;
        
        phi.resize(n_phi); tmp_vec1_n1.resize(n1); tmp_vec2_n1.resize(n1);
        tmp_vec3_n2.resize(n2); temperature.resize(6);
        
        Bmat.reinit(6, _system.n_vars(), _elem.n_nodes()); // six stress-strain components
        const std::vector<std::vector<RealVectorValue> >& dphi = _fe->get_dphi();
        
        Real temperature_value, ref_temperature;
        
        for (unsigned int qp=0; qp<JxW.size(); qp++) {
            this->initialize_membrane_strain_operator(dphi, qp, Bmat);
            
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
    StructuralElement2D::initialize_membrane_strain_operator
    (const std::vector<std::vector<RealVectorValue> >& dphi,
     const unsigned int qp,
     FEMOperatorMatrix& Bmat) {
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
    
    
}



#endif  // __MAST_structural_element_2D_h__
