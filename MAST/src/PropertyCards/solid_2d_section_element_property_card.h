//
//  solid_2d_section_element_property_card.h
//  MAST
//
//  Created by Manav Bhatia on 1/30/14.
//  Copyright (c) 2014 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_solid_2d_section_element_property_card_h__
#define __MAST_solid_2d_section_element_property_card_h__

// MAST includes
#include "PropertyCards/element_property_card_2D.h"


namespace MAST {
    
    class Solid2DSectionElementPropertyCard : public MAST::ElementPropertyCard2D {
    public:
        Solid2DSectionElementPropertyCard(unsigned int pid):
        MAST::ElementPropertyCard2D(pid),
        _material(NULL)
        { }
        
        /*!
         *   virtual destructor
         */
        virtual ~Solid2DSectionElementPropertyCard() { }
        
        
        class SectionIntegratedExtensionStiffnessMatrix: public MAST::FieldFunction<DenseMatrix<Real> > {
        public:
            SectionIntegratedExtensionStiffnessMatrix(MAST::FieldFunction<DenseMatrix<Real> > &mat,
                                                      MAST::FieldFunction<Real>& h);
            
            virtual ~SectionIntegratedExtensionStiffnessMatrix() { }
            
            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& v) const;
            
            virtual void partial_derivative (const MAST::SensitivityParameters& par,
                                             const Point& p, const Real t, DenseMatrix<Real>& v) const;
            
            virtual void total_derivative (const MAST::SensitivityParameters& par,
                                           const Point& p, const Real t, DenseMatrix<Real>& v) const;
            
        protected:
            
            MAST::FieldFunction<DenseMatrix<Real> > &_material_stiffness;
            MAST::FieldFunction<Real> &_h;
        };
        
        
        
        class SectionIntegratedExtensionBendingStiffnessMatrix: public MAST::FieldFunction<DenseMatrix<Real> > {
        public:
            SectionIntegratedExtensionBendingStiffnessMatrix(MAST::FieldFunction<DenseMatrix<Real> > &mat,
                                                             MAST::FieldFunction<Real>& h);
            
            virtual ~SectionIntegratedExtensionBendingStiffnessMatrix() { }
            
            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& v) const;
            
            virtual void partial_derivative (const MAST::SensitivityParameters& par,
                                             const Point& p, const Real t, DenseMatrix<Real>& v) const;
            
            virtual void total_derivative (const MAST::SensitivityParameters& par,
                                           const Point& p, const Real t, DenseMatrix<Real>& v) const;
            
        protected:
            
            MAST::FieldFunction<DenseMatrix<Real> > &_material_stiffness;
            MAST::FieldFunction<Real> &_h;
        };
        
        
        class SectionIntegratedBendingStiffnessMatrix: public MAST::FieldFunction<DenseMatrix<Real> > {
        public:
            SectionIntegratedBendingStiffnessMatrix(MAST::FieldFunction<DenseMatrix<Real> > &mat,
                                                    MAST::FieldFunction<Real>& h);
            
            virtual ~SectionIntegratedBendingStiffnessMatrix() { }
            
            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& v) const;
            
            virtual void partial_derivative (const MAST::SensitivityParameters& par,
                                             const Point& p, const Real t, DenseMatrix<Real>& v) const;
            
            virtual void total_derivative (const MAST::SensitivityParameters& par,
                                           const Point& p, const Real t, DenseMatrix<Real>& v) const;
            
        protected:
            
            MAST::FieldFunction<DenseMatrix<Real> > &_material_stiffness;
            MAST::FieldFunction<Real> &_h;
        };
        
        
        
        class SectionIntegratedInertiaMatrix: public MAST::FieldFunction<DenseMatrix<Real> > {
        public:
            SectionIntegratedInertiaMatrix(MAST::FieldFunction<DenseMatrix<Real> > &mat,
                                           MAST::FieldFunction<Real>& h);
            
            virtual ~SectionIntegratedInertiaMatrix() { }
            
            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& v) const;
            
            virtual void partial_derivative (const MAST::SensitivityParameters& par,
                                             const Point& p, const Real t, DenseMatrix<Real>& v) const;
            
            virtual void total_derivative (const MAST::SensitivityParameters& par,
                                           const Point& p, const Real t, DenseMatrix<Real>& v) const;
            
        protected:
            
            MAST::FieldFunction<DenseMatrix<Real> > &_material_inertia;
            MAST::FieldFunction<Real> &_h;
        };
        
        
        
        class SectionIntegratedThermalExpansionMatrix: public MAST::FieldFunction<DenseMatrix<Real> > {
        public:
            SectionIntegratedThermalExpansionMatrix(MAST::FieldFunction<DenseMatrix<Real> > &mat_stiff,
                                                    MAST::FieldFunction<DenseMatrix<Real> > &mat_expansion,
                                                    MAST::FieldFunction<Real>& h);
            
            virtual ~SectionIntegratedThermalExpansionMatrix() { }
            
            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& v) const;
            
            virtual void partial_derivative (const MAST::SensitivityParameters& par,
                                             const Point& p, const Real t, DenseMatrix<Real>& v) const;
            
            virtual void total_derivative (const MAST::SensitivityParameters& par,
                                           const Point& p, const Real t, DenseMatrix<Real>& v) const;
            
        protected:
            
            MAST::FieldFunction<DenseMatrix<Real> > &_material_stiffness;
            MAST::FieldFunction<DenseMatrix<Real> > &_material_expansion;
            MAST::FieldFunction<Real> &_h;
        };
        
        
        
        
        class SectionIntegratedPrestressMatrix: public MAST::FieldFunction<DenseMatrix<Real> > {
        public:
            SectionIntegratedPrestressMatrix(MAST::FieldFunction<DenseMatrix<Real> > &prestress,
                                             MAST::FieldFunction<Real>& h);
            
            virtual ~SectionIntegratedPrestressMatrix() { }
            
            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& v) const;
            
            virtual void partial_derivative (const MAST::SensitivityParameters& par,
                                             const Point& p, const Real t, DenseMatrix<Real>& v) const;
            
            virtual void total_derivative (const MAST::SensitivityParameters& par,
                                           const Point& p, const Real t, DenseMatrix<Real>& v) const;
            
        protected:
            
            MAST::FieldFunction<DenseMatrix<Real> > &_prestress;
            MAST::FieldFunction<Real> &_h;
        };
        
        
        /*!
         *   dimension of the element for which this property is defined
         */
        virtual unsigned int dim() const {
            return 2;
        }
        
        /*!
         *   return true if the property is isotropic
         */
        virtual bool if_isotropic() const {
            return true;
        }
        
        /*!
         *    sets the material card
         */
        void set_material(MAST::MaterialPropertyCardBase& mat) {
            _material = &mat;
        }
        
        
        /*!
         *    returns a reference to the material
         */
        virtual const MAST::MaterialPropertyCardBase& get_material() const {
            libmesh_assert(_material); // make sure it has already been set
            return *_material;
        }
        
        
        /*!
         *   calculates the matrix in \par m of type \par t.
         */
        virtual void calculate_matrix(const Elem& elem,
                                      MAST::ElemenetPropertyMatrixType t,
                                      DenseMatrix<Real>& m) const;
        
        /*!
         *   calculates the sensitivity of matrix in \par m of type \par t.
         */
        virtual void calculate_matrix_sensitivity(const Elem& elem,
                                                  MAST::ElemenetPropertyMatrixType t,
                                                  DenseMatrix<Real>& m,
                                                  const MAST::SensitivityParameters& p) const;
        
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
         *   material property card
         */
        MAST::MaterialPropertyCardBase *_material;
    };
    
}



#endif // __MAST_solid_2d_section_element_property_card_h__
