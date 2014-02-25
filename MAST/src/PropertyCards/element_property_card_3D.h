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
            SectionIntegratedStiffnessMatrix(MAST::FieldFunction<DenseMatrix<Real> > *mat);
            
            SectionIntegratedStiffnessMatrix(const MAST::ElementPropertyCard3D::SectionIntegratedStiffnessMatrix& f):
            MAST::FieldFunction<DenseMatrix<Real> >(f),
            _material_stiffness(f._material_stiffness->clone().release()) {
                _functions.insert(_material_stiffness);
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real> > > clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real> > >
                (new MAST::ElementPropertyCard3D::SectionIntegratedStiffnessMatrix(*this));
            }

            virtual ~SectionIntegratedStiffnessMatrix() { delete _material_stiffness;}
            
            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                             const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                           const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
        protected:
            
            MAST::FieldFunction<DenseMatrix<Real> > *_material_stiffness;
        };
        
        
        
        class SectionIntegratedInertiaMatrix: public MAST::FieldFunction<DenseMatrix<Real> > {
        public:
            SectionIntegratedInertiaMatrix(MAST::FieldFunction<DenseMatrix<Real> > *mat);
            
            SectionIntegratedInertiaMatrix(const MAST::ElementPropertyCard3D::SectionIntegratedInertiaMatrix& f):
            MAST::FieldFunction<DenseMatrix<Real> >(f),
            _material_inertia(f._material_inertia->clone().release()) {
                _functions.insert(_material_inertia);
            }

            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real> > > clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real> > >
                (new MAST::ElementPropertyCard3D::SectionIntegratedInertiaMatrix(*this));
            }

            virtual ~SectionIntegratedInertiaMatrix() { delete _material_inertia;}
            
            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                             const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                           const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
        protected:
            
            MAST::FieldFunction<DenseMatrix<Real> > *_material_inertia;
        };
        
        
        
        class SectionIntegratedThermalExpansionMatrix: public MAST::FieldFunction<DenseMatrix<Real> > {
        public:
            SectionIntegratedThermalExpansionMatrix(MAST::FieldFunction<DenseMatrix<Real> > *mat_stiff,
                                                    MAST::FieldFunction<DenseMatrix<Real> > *mat_expansion);
            
            SectionIntegratedThermalExpansionMatrix(const MAST::ElementPropertyCard3D::
                                                    SectionIntegratedThermalExpansionMatrix& f):
            MAST::FieldFunction<DenseMatrix<Real> >(f),
            _material_stiffness(f._material_stiffness->clone().release()),
            _material_expansion(f._material_expansion->clone().release()) {
                _functions.insert(_material_stiffness);
                _functions.insert(_material_expansion);
            }

            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real> > > clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real> > >
                (new MAST::ElementPropertyCard3D::SectionIntegratedThermalExpansionMatrix(*this));
            }

            virtual ~SectionIntegratedThermalExpansionMatrix() {
                delete _material_stiffness;
                delete _material_expansion;
            }
            
            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                             const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                           const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
        protected:
            
            MAST::FieldFunction<DenseMatrix<Real> > *_material_stiffness;
            MAST::FieldFunction<DenseMatrix<Real> > *_material_expansion;
        };
        
        
        
        
        class SectionIntegratedPrestressAMatrix: public MAST::SectionIntegratedPrestressMatrixBase {
        public:
            SectionIntegratedPrestressAMatrix(MAST::FieldFunction<DenseMatrix<Real> > *prestress);
            
            SectionIntegratedPrestressAMatrix(const MAST::ElementPropertyCard3D::SectionIntegratedPrestressAMatrix& f):
            MAST::SectionIntegratedPrestressMatrixBase(f),
            _prestress(f._prestress->clone().release()) {
                _functions.insert(_prestress);
            }

            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real> > > clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real> > >
                (new MAST::ElementPropertyCard3D::SectionIntegratedPrestressAMatrix(*this));
            }

            virtual ~SectionIntegratedPrestressAMatrix() { delete _prestress;}
            
            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                             const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                           const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void convert_to_vector(const DenseMatrix<Real>& m, DenseVector<Real>& v) const;
            
        protected:
            
            MAST::FieldFunction<DenseMatrix<Real> > *_prestress;
        };

        
        /*!
         *   dimension of the element for which this property is defined
         */
        virtual unsigned int dim() const {
            return 3;
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
         *   returns a function to evaluate the specified quantitys
         *   type \par t.
         */
        virtual std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real> > >
        get_property(MAST::ElemenetPropertyMatrixType t,
                     const MAST::StructuralElementBase& e) const;

        /*!
         *  returns true if the property card depends on the function \p f
         */
        virtual bool depends_on(const MAST::FieldFunctionBase& f) const {
            return _material->depends_on(f) ||            // check if the material property depends on the function
            MAST::ElementPropertyCardBase::depends_on(f); // check with this property card
        }
        
        
        protected:
        
        /*!
         *    pointer to the material property card
         */
        MAST::MaterialPropertyCardBase* _material;
    };
}


#endif // __MAST_element_property_card_3D_h__
