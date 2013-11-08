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
        ElementPropertyCard3D():
        MAST::ElementPropertyCardBase(),
        _material(NULL)
        { }
        
        /*!
         *   virtual destructor
         */
        virtual ~ElementPropertyCard3D() { }
        
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
         *    sets the material card
         */
        void set_material(MAST::MaterialPropertyCardBase& mat) {
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
                                      DenseVector<Real>& v) const;
        
        /*!
         *    initializes the matrix to the prestress in the element
         */
        virtual void prestress_matrix(MAST::ElemenetPropertyMatrixType t,
                                      DenseMatrix<Real>& m) const;

        
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

                case MAST::SECTION_INTEGRATED_MATERIAL_THERMAL_EXPANSION_MATRIX:
                    _material->calculate_3d_matrix(MAST::MATERIAL_THERMAL_EXPANSION_MATRIX, m);
                    break;
                    
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
                    
                case MAST::SECTION_INTEGRATED_MATERIAL_THERMAL_EXPANSION_MATRIX:
                    _material->calculate_3d_matrix_sensitivity(MAST::MATERIAL_THERMAL_EXPANSION_MATRIX, m, p);
                    break;
                    
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
                                              DenseVector<Real>& v) const {
    if (_prestress.size() == 0)
        v.resize(6); // zero, if the stress has not been defined
    else
        v = _prestress;
}


/*!
 *    initializes the matrix to the prestress in the element
 */
inline void
MAST::ElementPropertyCard3D::prestress_matrix(MAST::ElemenetPropertyMatrixType t,
                                              DenseMatrix<Real>& m) const {
    m.resize(3, 3);
    if (_prestress.size() == 6) {
        for (unsigned int i=0; i<3; i++)
            m(i,i) = _prestress(i);
        m(0,1) = _prestress(3);  // tau_xy
        m(1,0) = _prestress(3);  // tau_xy
        m(1,2) = _prestress(4);  // tau_yz
        m(2,1) = _prestress(4);  // tau_yz
        m(0,2) = _prestress(4);  // tau_zx
        m(2,0) = _prestress(4);  // tau_zx
    }
}



#endif // __MAST_element_property_card_3D_h__
