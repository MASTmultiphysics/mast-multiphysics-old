//
//  element_property_card_base.h
//  MAST
//
//  Created by Manav Bhatia on 10/15/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_element_property_card_base_h__
#define __MAST_element_property_card_base_h__

// MAST includes
#include "PropertyCards/property_card_base.h"

// libMesh includes
#include "libmesh/elem.h"


namespace MAST
{
    enum StrainType {
        LINEAR_STRAIN,
        VON_KARMAN_STRAIN
    };
    
    enum ElemenetPropertyMatrixType {
        SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX,
        SECTION_INTEGRATED_MATERIAL_STIFFNESS_B_MATRIX_1D_TRANSVERSE,
        SECTION_INTEGRATED_MATERIAL_STIFFNESS_D_MATRIX_1D_TRANSVERSE,
        SECTION_INTEGRATED_MATERIAL_STIFFNESS_B_MATRIX_1D_CHORDWISE,
        SECTION_INTEGRATED_MATERIAL_STIFFNESS_D_MATRIX_1D_CHORDWISE,
        SECTION_INTEGRATED_MATERIAL_STIFFNESS_B_MATRIX,
        SECTION_INTEGRATED_MATERIAL_STIFFNESS_D_MATRIX,
        SECTION_INTEGRATED_MATERIAL_DAMPING_MATRIX,
        SECTION_INTEGRATED_MATERIAL_INERTIA_MATRIX,
        SECTION_INTEGRATED_MATERIAL_THERMAL_EXPANSION_MATRIX,
        SECTION_INTEGRATED_MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX
    };
    
    
    class ElementPropertyCardBase: public MAST::PropertyCardBase {
        
    public:
        ElementPropertyCardBase():
        MAST::PropertyCardBase(),
        _strain_type(MAST::LINEAR_STRAIN),
        _diagonal_mass(false)
        { }
        
        /*!
         *   virtual destructor
         */
        virtual ~ElementPropertyCardBase() { }
        
        
        /*!
         *    returns the extra quadrature order (on top of the system) that 
         *    this element should use. By default this is zero, and can be
         *    changed by the derived classes
         */
        virtual unsigned int extra_quadrature_order(const Elem& elem,
                                                    const FEType& fe) const {
            return 0;
        }
        
        /*!
         *   calculates the material matrix in \par m of type \par t.
         */
        virtual void calculate_matrix(const Elem& elem,
                                      MAST::ElemenetPropertyMatrixType t,
                                      DenseMatrix<Real>& m) const = 0;

        /*!
         *   calculates the material matrix in \par m of type \par t.
         */
        virtual void calculate_matrix_sensitivity(const Elem& elem,
                                                  MAST::ElemenetPropertyMatrixType t,
                                                  DenseMatrix<Real>& m,
                                                  const MAST::SensitivityParameters& params) const = 0;

        /*!
         *    sets the type of strain to be used, which is LINEAR_STRAIN by
         *    default
         */
        void set_strain(MAST::StrainType strain) {
            _strain_type = strain;
        }
        
        
        /*!
         *    returns the type of strain to be used for this element
         */
        const MAST::StrainType strain_type() const {
            return _strain_type;
        }

        
        /*!
         *    sets the mass matrix to be diagonal or consistent
         */
        void set_diagonal_mass_matrix(bool m) {
            _diagonal_mass = m;
        }
        
        
        /*!
         *    returns the type of strain to be used for this element
         */
        bool if_diagonal_mass_matrix() const {
            return _diagonal_mass;
        }

        
        /*!
         *    returns true if the element prestress has been specified, false
         *    otherwise
         */
        virtual bool if_prestressed() const {
            if (_prestress.size())
                return true;
            else
                return false;
        }

        
        /*!
         *    initializes the prestress in the element. The vector should 
         *    contain six stress components: sigma_xx, sigma_yy, sigma_zz,
         *    tau_xy, tau_yz, tau_zx
         */
        virtual void prestress(const DenseVector<Real>& v) {
            libmesh_assert_equal_to(v.size(), 6);
            _prestress = v;
        }

        /*!
         *    initializes the vector to the prestress in the element
         */
        virtual void prestress_vector(MAST::ElemenetPropertyMatrixType t,
                                      DenseVector<Real>& v) const = 0;
        
        
        /*!
         *    initializes the matrix to the prestress in the element
         */
        virtual void prestress_matrix(MAST::ElemenetPropertyMatrixType t,
                                      DenseMatrix<Real>& m) const = 0;
        
    protected:
        
        /*!
         *    type of nonlinear strain to be used for analysis
         */
        MAST::StrainType _strain_type;

        /*!
         *    flag to use a diagonal mass matrix. By default, this is false
         */
        bool _diagonal_mass;
        
        /*!
         *   element prestress: six stress components: sigma_xx, sigma_yy, 
         *   sigma_zz, tau_xy, tau_yz, tau_zx
         */
        DenseVector<Real> _prestress;
    };
    
    
}



#endif // __MAST_element_property_card_base_h__
