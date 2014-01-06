//
//  element_property_card_3D.h
//  MAST
//
//  Created by Manav Bhatia on 10/24/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_element_property_card_3D_h__
#define __MAST_element_property_card_3D_h__

// MAST includes
#include "PropertyCards/element_property_card_base.h"
#include "PropertyCards/material_property_card_base.h"


namespace MAST
{
    class ElementPropertyCard3D: public MAST::ElementPropertyCardBase {
        
    public:
    public:
        ElementPropertyCard3D(unsigned int pid):
        MAST::ElementPropertyCardBase(pid),
        _material(NULL)
        { }
        
        /*!
         *   virtual destructor
         */
        virtual ~ElementPropertyCard3D() { }
        
        /*!
         *   dimension of the element for which this property is defined
         */
        virtual unsigned int dim() const {
            return 3;
        }

        /*!
         *   calculates the material matrix in \par m of type \par t.
         */
        virtual void calculate_matrix(const Elem& elem,
                                      MAST::ElemenetPropertyMatrixType t,
                                      DenseMatrix<Real>& m) const;
        
        /*!
         *   calculates the material matrix in \par m of type \par t.
         */
        virtual void calculate_matrix_sensitivity(const Elem& elem,
                                                  MAST::ElemenetPropertyMatrixType t,
                                                  DenseMatrix<Real>& m,
                                                  const MAST::SensitivityParameters& params) const;

        /*!
         *   return true if the property is isotropic
         */
        virtual bool if_isotropic() const {
            return true;
        }

        
        /*!
         *    sets the material card
         */
        virtual void set_material(MAST::MaterialPropertyCardBase& mat) {
            _material = &mat;
        }
        
        
        /*!
         *    returns a reference to the material
         */
        const MAST::MaterialPropertyCardBase& get_material() const {
            libmesh_assert(_material); // make sure it has already been set
            return *_material;
        }
        
        /*!
         *    initializes the vector to the prestress in the element
         */
        virtual void prestress_vector(MAST::ElemenetPropertyMatrixType t,
                                      const DenseMatrix<Real>& T,
                                      DenseVector<Real>& v) const;
        
        /*!
         *    initializes the vector to the sensitivity of prestress in the element
         */
        virtual void prestress_vector_sensitivity(MAST::ElemenetPropertyMatrixType t,
                                                  const DenseMatrix<Real>& T,
                                                  DenseVector<Real>& v,
                                                  const MAST::SensitivityParameters& p) const;

        /*!
         *    initializes the matrix to the prestress in the element
         */
        virtual void prestress_matrix(MAST::ElemenetPropertyMatrixType t,
                                      const DenseMatrix<Real>& T,
                                      DenseMatrix<Real>& m) const;

        /*!
         *    initializes the matrix to the sensitivity of prestress in the element
         */
        virtual void prestress_matrix_sensitivity(MAST::ElemenetPropertyMatrixType t,
                                                  const DenseMatrix<Real>& T,
                                                  DenseMatrix<Real>& m,
                                                  const MAST::SensitivityParameters& p) const;
        
        /*!
         *  returns true if the property card depends on the function \p f
         */
        virtual bool depends_on(const FunctionBase& f) const {
            return _material->depends_on(f) ||            // check if the material property depends on the function
            MAST::ElementPropertyCardBase::depends_on(f); // check with this property card
        }
        
        
        /*!
         *  returns true if the property card depends on the functions in \p p
         */
        virtual bool depends_on(const MAST::SensitivityParameters& p) const {
            return _material->depends_on(p) ||            // check if the material property depends on the function
            MAST::ElementPropertyCardBase::depends_on(p); // check with this property card
        }

        protected:
        
        /*!
         *    pointer to the material property card
         */
        MAST::MaterialPropertyCardBase* _material;
    };
}


inline void
MAST::ElementPropertyCard3D::calculate_matrix(const libMesh::Elem &elem,
                                              MAST::ElemenetPropertyMatrixType t,
                                              DenseMatrix<Real>& m) const
{
    libmesh_assert(_material); // should have been set
    
    switch (elem.dim()) {
        case 3:
            switch (t) {
                case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX:
                    _material->calculate_3d_matrix(MAST::MATERIAL_STIFFNESS_MATRIX, m);
                    break;
                    
                case MAST::SECTION_INTEGRATED_MATERIAL_INERTIA_MATRIX:
                    _material->calculate_3d_matrix(MAST::MATERIAL_INERTIA_MATRIX, m);
                    // now scale the rotation dofs with small factors
                    for (unsigned int i=0; i<3; i++) {
                        m(i+3, i+3) *= 1.0e-6;
                    }
                    break;

                    break;

                case MAST::SECTION_INTEGRATED_MATERIAL_THERMAL_EXPANSION_A_MATRIX: {
                    DenseMatrix<Real> a;
                    _material->calculate_3d_matrix(MAST::MATERIAL_STIFFNESS_MATRIX, m);
                    _material->calculate_3d_matrix(MAST::MATERIAL_THERMAL_EXPANSION_MATRIX, a);
                    m.right_multiply(a);
                }
                    break;
                    
                case MAST::SECTION_INTEGRATED_MATERIAL_THERMAL_EXPANSION_B_MATRIX:
                default:
                    libmesh_error(); // others need to be implemented
                    break;
            }
            break;
            
        case 1:
        case 2:
        default:
            libmesh_error();
            break;
    }
}



inline void
MAST::ElementPropertyCard3D::calculate_matrix_sensitivity(const libMesh::Elem &elem,
                                                          MAST::ElemenetPropertyMatrixType t,
                                                          DenseMatrix<Real>& m,
                                                          const MAST::SensitivityParameters& p) const
{
    libmesh_assert(_material); // should have been set
    
    // currently only first order sensitivity is provided
    libmesh_assert_equal_to(p.total_order(), 1);
    
    switch (elem.dim()) {
        case 3:
            switch (t) {
                case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX:
                    _material->calculate_3d_matrix_sensitivity(MAST::MATERIAL_STIFFNESS_MATRIX, m, p);
                    break;
                    
                case MAST::SECTION_INTEGRATED_MATERIAL_INERTIA_MATRIX:
                    _material->calculate_3d_matrix_sensitivity(MAST::MATERIAL_INERTIA_MATRIX, m, p);
                    // now scale the rotation dofs with small factors
                    for (unsigned int i=0; i<3; i++) {
                        m(i+3, i+3) *= 1.0e-6;
                    }
                    break;
                    
                    break;
                    
                case MAST::SECTION_INTEGRATED_MATERIAL_THERMAL_EXPANSION_A_MATRIX: {
                    DenseMatrix<Real> a, aprime;
                    _material->calculate_3d_matrix(MAST::MATERIAL_STIFFNESS_MATRIX, m);
                    _material->calculate_3d_matrix_sensitivity(MAST::MATERIAL_THERMAL_EXPANSION_MATRIX,
                                                               a, p);
                    m.right_multiply(a);
                }
                    break;

                case MAST::SECTION_INTEGRATED_MATERIAL_THERMAL_EXPANSION_B_MATRIX:
                default:
                    libmesh_error(); // others need to be implemented
                    break;
            }
            break;
            
        case 1:
        case 2:
        default:
            libmesh_error();
            break;
    }
}




inline void
MAST::ElementPropertyCard3D::prestress_vector(MAST::ElemenetPropertyMatrixType t,
                                              const DenseMatrix<Real>& T,
                                              DenseVector<Real>& v) const {
    v.resize(6); // zero, if the stress has not been defined
    if (_prestress.m() != 0) {
        for (unsigned int i=0; i<3; i++)
            v(i) = _prestress(i,i);
        v(3) = _prestress(0,1);  // tau_xy
        v(4) = _prestress(1,2);  // tau_yz
        v(5) = _prestress(0,2);  // tau_xz
    }
}



inline void
MAST::ElementPropertyCard3D::prestress_vector_sensitivity(MAST::ElemenetPropertyMatrixType t,
                                                          const DenseMatrix<Real>& T,
                                                          DenseVector<Real>& v,
                                                          const MAST::SensitivityParameters& p) const {
    v.resize(6); // sensitivity is always zero
}



inline void
MAST::ElementPropertyCard3D::prestress_matrix(MAST::ElemenetPropertyMatrixType t,
                                              const DenseMatrix<Real>& T,
                                              DenseMatrix<Real>& m) const {
    m.resize(3, 3);
    if (_prestress.m() != 0)
        m = _prestress;
}



inline void
MAST::ElementPropertyCard3D::prestress_matrix_sensitivity(MAST::ElemenetPropertyMatrixType t,
                                                          const DenseMatrix<Real>& T,
                                                          DenseMatrix<Real>& m,
                                                          const MAST::SensitivityParameters& p) const {
    m.resize(3, 3); // sensitivity is always zero
}



#endif // __MAST_element_property_card_3D_h__
