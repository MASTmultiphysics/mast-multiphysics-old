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
        
        
        
        class SectionIntegratedStiffnessMatrix: public MAST::FieldFunction<DenseMatrix<Real> > {
        public:
            SectionIntegratedStiffnessMatrix(MAST::FieldFunction<DenseMatrix<Real> > &mat);
            
            virtual ~SectionIntegratedStiffnessMatrix() { }
            
            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& v) const;
            
            virtual void partial_derivative (const MAST::SensitivityParameters& par,
                                             const Point& p, const Real t, DenseMatrix<Real>& v) const;
            
            virtual void total_derivative (const MAST::SensitivityParameters& par,
                                           const Point& p, const Real t, DenseMatrix<Real>& v) const;
            
        protected:
            
            MAST::FieldFunction<DenseMatrix<Real> > &_material_stiffness;
        };
        
        
        
        class SectionIntegratedInertiaMatrix: public MAST::FieldFunction<DenseMatrix<Real> > {
        public:
            SectionIntegratedInertiaMatrix(MAST::FieldFunction<DenseMatrix<Real> > &mat);
            
            virtual ~SectionIntegratedInertiaMatrix() { }
            
            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& v) const;
            
            virtual void partial_derivative (const MAST::SensitivityParameters& par,
                                             const Point& p, const Real t, DenseMatrix<Real>& v) const;
            
            virtual void total_derivative (const MAST::SensitivityParameters& par,
                                           const Point& p, const Real t, DenseMatrix<Real>& v) const;
            
        protected:
            
            MAST::FieldFunction<DenseMatrix<Real> > &_material_inertia;
        };
        
        
        
        class SectionIntegratedThermalExpansionMatrix: public MAST::FieldFunction<DenseMatrix<Real> > {
        public:
            SectionIntegratedThermalExpansionMatrix(MAST::FieldFunction<DenseMatrix<Real> > &mat_stiff,
                                                    MAST::FieldFunction<DenseMatrix<Real> > &mat_expansion);
            
            virtual ~SectionIntegratedThermalExpansionMatrix() { }
            
            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& v) const;
            
            virtual void partial_derivative (const MAST::SensitivityParameters& par,
                                             const Point& p, const Real t, DenseMatrix<Real>& v) const;
            
            virtual void total_derivative (const MAST::SensitivityParameters& par,
                                           const Point& p, const Real t, DenseMatrix<Real>& v) const;
            
        protected:
            
            MAST::FieldFunction<DenseMatrix<Real> > &_material_stiffness;
            MAST::FieldFunction<DenseMatrix<Real> > &_material_expansion;
        };
        
        
        
        
        class SectionIntegratedPrestressMatrix: public MAST::FieldFunction<DenseMatrix<Real> > {
        public:
            SectionIntegratedPrestressMatrix(MAST::FieldFunction<DenseMatrix<Real> > &prestress);
            
            virtual ~SectionIntegratedPrestressMatrix() { }
            
            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& v) const;
            
            virtual void partial_derivative (const MAST::SensitivityParameters& par,
                                             const Point& p, const Real t, DenseMatrix<Real>& v) const;
            
            virtual void total_derivative (const MAST::SensitivityParameters& par,
                                           const Point& p, const Real t, DenseMatrix<Real>& v) const;
            
        protected:
            
            MAST::FieldFunction<DenseMatrix<Real> > &_prestress;
        };


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
                case MAST::SECTION_INTEGRATED_MATERIAL_INERTIA_MATRIX:
                case MAST::SECTION_INTEGRATED_MATERIAL_THERMAL_EXPANSION_A_MATRIX:
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



#endif // __MAST_element_property_card_3D_h__
