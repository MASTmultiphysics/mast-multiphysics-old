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
        SECTION_INTEGRATED_MATERIAL_THERMAL_EXPANSION_A_MATRIX,
        SECTION_INTEGRATED_MATERIAL_THERMAL_EXPANSION_B_MATRIX,
        SECTION_INTEGRATED_MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX,
        SECTION_INTEGRATED_PRESTRESS_A_MATRIX,
        SECTION_INTEGRATED_PRESTRESS_B_MATRIX
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
        virtual MAST::BendingOperatorType bending_model(const libMesh::Elem& elem,
                                                        const libMesh::FEType& fe) const
        { libmesh_error(); }
        
        /*!
         *    returns the extra quadrature order (on top of the system) that
         *    this element should use. By default this is zero, and can be
         *    changed by the derived classes
         */
        virtual int extra_quadrature_order(const libMesh::Elem& elem,
                                           const libMesh::FEType& fe) const {
            return 0;
        }
        
        /*!
         *   returns a function to evaluate the specified quantitys
         *   type \par t.
         */
        virtual std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > >
        get_property(MAST::ElemenetPropertyMatrixType t,
                     const MAST::StructuralElementBase& e) const = 0;
        
        
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
            return this->contains("prestress");
        }
        
        
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
    };
    
    
    /*!
     *    base class for prestress matrix that requires conversion of the matrix to a
     *    vector form
     *
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
    class SectionIntegratedPrestressMatrixBase: public MAST::FieldFunction<DenseRealMatrix > {
    public:
        SectionIntegratedPrestressMatrixBase(const std::string& nm):
        MAST::FieldFunction<DenseRealMatrix >(nm)
        { }

        SectionIntegratedPrestressMatrixBase(const MAST::SectionIntegratedPrestressMatrixBase& f):
        MAST::FieldFunction<DenseRealMatrix >(f)
        { }

        virtual ~SectionIntegratedPrestressMatrixBase() { }
        
        virtual void convert_to_vector(const DenseRealMatrix& m, DenseRealVector& v) const = 0;
    };
    
    
}



#endif // __MAST_element_property_card_base_h__
