//
//  element_property_card_base.h
//  MAST
//
//  Created by Manav Bhatia on 10/15/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_element_property_card_base_h__
#define __MAST_element_property_card_base_h__

// libMesh includes
#include "libmesh/elem.h"

// MAST includes
#include "PropertyCards/property_card_base.h"
#include "StructuralElems/bending_operator.h"



namespace MAST
{
    // forward decleration
    class MaterialPropertyCardBase;
    
    enum StrainType {
        LINEAR_STRAIN,
        VON_KARMAN_STRAIN
    };
 
    enum ElemenetPropertyMatrixType {
        SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX,
        SECTION_INTEGRATED_MATERIAL_STIFFNESS_B_MATRIX,
        SECTION_INTEGRATED_MATERIAL_STIFFNESS_D_MATRIX,
        SECTION_INTEGRATED_MATERIAL_DAMPING_MATRIX,
        SECTION_INTEGRATED_MATERIAL_INERTIA_MATRIX,
        SECTION_INTEGRATED_MATERIAL_THERMAL_EXPANSION_MATRIX,
        SECTION_INTEGRATED_MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX
    };
    
    
    class ElementPropertyCardBase: public MAST::PropertyCardBase {
        
    public:
        ElementPropertyCardBase(unsigned int pid):
        MAST::PropertyCardBase(),
        _pid(pid),
        _strain_type(MAST::LINEAR_STRAIN),
        _diagonal_mass(false)
        { }
        
        /*!
         *   virtual destructor
         */
        virtual ~ElementPropertyCardBase() { }

        /*!
         *   returns the bending model to be used for the element. Should be 
         *   reimplemented in the derived classes
         */
        virtual MAST::BendingOperatorType bending_model(const Elem& elem,
                                                 const FEType& fe) const
        { libmesh_error(); }
        
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
         *    returns the id for this card
         */
        unsigned int id() const {
            return _pid;
        }
        
        
        /*!
         *   return true if the property is isotropic
         */
        virtual bool if_isotropic() const = 0;

        
        /*!
         *   return the material property. This needs to be reimplemented 
         *   for individual card type, and should be used only for isotropic 
         *   cards.
         */
        virtual const MAST::MaterialPropertyCardBase& get_material() const {
            libmesh_error();
        }

        
        /*!
         *   dimension of the element for which this property is defined
         */
        virtual unsigned int dim() const = 0;
        
        
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
            if (_prestress.m() == 0)
                return false;
            else
                return true;
        }

        
        /*!
         *    initializes the prestress in the element. The vector should 
         *    contain six stress components: sigma_xx, sigma_yy, sigma_zz,
         *    tau_xy, tau_yz, tau_zx. It is assumed that the stress values 
         *    are defined in the global coordinate system.
         */
        virtual void prestress(const DenseVector<Real>& v) {
            libmesh_assert_equal_to(v.size(), 6);
            _prestress.resize(3, 3);
            _prestress(0,0) = v(0); // sigma_xx
            _prestress(1,1) = v(1); // sigma_yy
            _prestress(2,2) = v(2); // sigma_zz
            _prestress(0,1) = v(3); // sigma_xy
            _prestress(1,0) = v(3);
            _prestress(1,2) = v(4); // sigma_yz
            _prestress(2,1) = v(4);
            _prestress(0,2) = v(5); // sigma_zx
            _prestress(2,0) = v(5);
            
        }

        /*!
         *    initializes the vector to the prestress in the element. The 
         *    stress value is defined in the global coordinate system. 
         *    Hence, the given matrix \par T is used to transform the vector
         *    to the local coordinate defined for \par T. Note that 
         *    T_ij = V_i^t . Vn_j, where
         *    V_i are the unit vectors of the global cs, and Vn_j are the
         *    unit vectors of the local cs. To transform a vector from global to
         *    local cs,    an_j = T^t a_i, and the reverse transformation is
         *    obtained as  a_j  = T  an_i
         */
        virtual void prestress_vector(MAST::ElemenetPropertyMatrixType t,
                                      const DenseMatrix<Real>& T,
                                      DenseVector<Real>& v) const = 0;
        
        
        /*!
         *    initializes the matrix to the prestress in the element. The
         *    stress value is defined in the global coordinate system.
         *    Hence, the given matrix \par T is used to transform the vector
         *    to the local coordinate defined for \par T. Note that
         *    T_ij = V_i^t . Vn_j, where
         *    V_i are the unit vectors of the global cs, and Vn_j are the
         *    unit vectors of the local cs. To transform a vector from global to
         *    local cs,    an_j = T^t a_i, and the reverse transformation is
         *    obtained as  a_j  = T  an_i
         */
        virtual void prestress_matrix(MAST::ElemenetPropertyMatrixType t,
                                      const DenseMatrix<Real>& T,
                                      DenseMatrix<Real>& m) const = 0;
        
        
    protected:
        
        /*!
         *    property card id
         */
        unsigned int _pid;
        
        /*!
         *    type of nonlinear strain to be used for analysis
         */
        MAST::StrainType _strain_type;

        /*!
         *    flag to use a diagonal mass matrix. By default, this is false
         */
        bool _diagonal_mass;
        
        /*!
         *   element prestress tensor: six stress components: sigma_xx, sigma_yy,
         *   sigma_zz, tau_xy, tau_yz, tau_zx
         */
        DenseMatrix<Real> _prestress;
    };
    
    
}



#endif // __MAST_element_property_card_base_h__
