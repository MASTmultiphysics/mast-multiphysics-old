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
            SectionIntegratedExtensionStiffnessMatrix(MAST::FieldFunction<DenseMatrix<Real> > *mat,
                                                      MAST::FieldFunction<Real> *h);
            
            SectionIntegratedExtensionStiffnessMatrix(const MAST::Solid2DSectionElementPropertyCard::SectionIntegratedExtensionStiffnessMatrix& f):
            MAST::FieldFunction<DenseMatrix<Real> >(f),
            _material_stiffness(f._material_stiffness->clone().release()),
            _h(f._h->clone().release())
            { }

            
            virtual ~SectionIntegratedExtensionStiffnessMatrix() {
                delete _material_stiffness;
                delete _h;
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>> clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>>
                (new MAST::Solid2DSectionElementPropertyCard::SectionIntegratedExtensionStiffnessMatrix(*this));
            }

            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
        protected:
            
            MAST::FieldFunction<DenseMatrix<Real> > *_material_stiffness;
            MAST::FieldFunction<Real> *_h;
        };
        
        
        
        class SectionIntegratedExtensionBendingStiffnessMatrix: public MAST::FieldFunction<DenseMatrix<Real> > {
        public:
            SectionIntegratedExtensionBendingStiffnessMatrix(MAST::FieldFunction<DenseMatrix<Real> > *mat,
                                                             MAST::FieldFunction<Real> *h);
            
            SectionIntegratedExtensionBendingStiffnessMatrix(const MAST::Solid2DSectionElementPropertyCard::SectionIntegratedExtensionBendingStiffnessMatrix& f):
            MAST::FieldFunction<DenseMatrix<Real> >(f),
            _material_stiffness(f._material_stiffness->clone().release()),
            _h(f._h->clone().release())
            { }
            
            virtual ~SectionIntegratedExtensionBendingStiffnessMatrix() {
                delete _material_stiffness;
                delete _h;
            }

            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>> clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>>
                (new MAST::Solid2DSectionElementPropertyCard::SectionIntegratedExtensionBendingStiffnessMatrix(*this));
            }
            
            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
        protected:
            
            MAST::FieldFunction<DenseMatrix<Real> > *_material_stiffness;
            MAST::FieldFunction<Real> *_h;
        };
        
        
        class SectionIntegratedBendingStiffnessMatrix: public MAST::FieldFunction<DenseMatrix<Real> > {
        public:
            SectionIntegratedBendingStiffnessMatrix(MAST::FieldFunction<DenseMatrix<Real> > *mat,
                                                    MAST::FieldFunction<Real> *h);
            
            SectionIntegratedBendingStiffnessMatrix(const MAST::Solid2DSectionElementPropertyCard::SectionIntegratedBendingStiffnessMatrix& f):
            MAST::FieldFunction<DenseMatrix<Real> >(f),
            _material_stiffness(f._material_stiffness->clone().release()),
            _h(f._h->clone().release())
            { }
            

            virtual ~SectionIntegratedBendingStiffnessMatrix() {
                delete _material_stiffness;
                delete _h;
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>> clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>>
                (new MAST::Solid2DSectionElementPropertyCard::SectionIntegratedBendingStiffnessMatrix(*this));
            }

            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
        protected:
            
            MAST::FieldFunction<DenseMatrix<Real> > *_material_stiffness;
            MAST::FieldFunction<Real> *_h;
        };
        
        
        
        class SectionIntegratedInertiaMatrix: public MAST::FieldFunction<DenseMatrix<Real> > {
        public:
            SectionIntegratedInertiaMatrix(MAST::FieldFunction<Real> *rho,
                                           MAST::FieldFunction<Real> *h);
            
            SectionIntegratedInertiaMatrix(const MAST::Solid2DSectionElementPropertyCard::SectionIntegratedInertiaMatrix& f):
            MAST::FieldFunction<DenseMatrix<Real> >(f),
            _rho(f._rho->clone().release()),
            _h(f._h->clone().release())
            { }
            

            virtual ~SectionIntegratedInertiaMatrix() {
                delete _rho;
                delete _h;
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>> clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>>
                (new MAST::Solid2DSectionElementPropertyCard::SectionIntegratedInertiaMatrix(*this));
            }

            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
        protected:
            
            MAST::FieldFunction<Real> *_rho, *_h;
        };
        
        
        
        class SectionIntegratedThermalExpansionAMatrix: public MAST::FieldFunction<DenseMatrix<Real> > {
        public:
            SectionIntegratedThermalExpansionAMatrix(MAST::FieldFunction<DenseMatrix<Real> > *mat_stiff,
                                                     MAST::FieldFunction<DenseMatrix<Real> > *mat_expansion,
                                                     MAST::FieldFunction<Real> *h);
            
            SectionIntegratedThermalExpansionAMatrix(const MAST::Solid2DSectionElementPropertyCard::SectionIntegratedThermalExpansionAMatrix& f):
            MAST::FieldFunction<DenseMatrix<Real> >(f),
            _material_stiffness(f._material_stiffness->clone().release()),
            _material_expansion(f._material_expansion->clone().release()),
            _h(f._h->clone().release())
            { }
            
            
            
            virtual ~SectionIntegratedThermalExpansionAMatrix() {
                delete _material_stiffness;
                delete _material_expansion;
                delete _h;
            }
            
            
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>> clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>>
                (new MAST::Solid2DSectionElementPropertyCard::SectionIntegratedThermalExpansionAMatrix(*this));
            }
            
            
            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
        protected:
            
            MAST::FieldFunction<DenseMatrix<Real> > *_material_stiffness;
            MAST::FieldFunction<DenseMatrix<Real> > *_material_expansion;
            MAST::FieldFunction<Real> *_h;
        };

        
        
        class SectionIntegratedThermalExpansionBMatrix: public MAST::FieldFunction<DenseMatrix<Real> > {
        public:
            SectionIntegratedThermalExpansionBMatrix(MAST::FieldFunction<DenseMatrix<Real> > *mat_stiff,
                                                    MAST::FieldFunction<DenseMatrix<Real> > *mat_expansion,
                                                    MAST::FieldFunction<Real> *h);
            
            SectionIntegratedThermalExpansionBMatrix(const MAST::Solid2DSectionElementPropertyCard::SectionIntegratedThermalExpansionBMatrix& f):
            MAST::FieldFunction<DenseMatrix<Real> >(f),
            _material_stiffness(f._material_stiffness->clone().release()),
            _material_expansion(f._material_expansion->clone().release()),
            _h(f._h->clone().release())
            { }
            
            
            virtual ~SectionIntegratedThermalExpansionBMatrix() {
                delete _material_stiffness;
                delete _material_expansion;
                delete _h;
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>> clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>>
                (new MAST::Solid2DSectionElementPropertyCard::SectionIntegratedThermalExpansionBMatrix(*this));
            }
            
            
            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            
        protected:
            
            MAST::FieldFunction<DenseMatrix<Real> > *_material_stiffness;
            MAST::FieldFunction<DenseMatrix<Real> > *_material_expansion;
            MAST::FieldFunction<Real> *_h;
        };

        
        
        
        class SectionIntegratedPrestressAMatrix: public MAST::SectionIntegratedPrestressMatrixBase {
        public:
            SectionIntegratedPrestressAMatrix(MAST::FieldFunction<DenseMatrix<Real> > *prestress,
                                              MAST::FieldFunction<DenseMatrix<Real> > *T,
                                              MAST::FieldFunction<Real> *h);

            SectionIntegratedPrestressAMatrix(const MAST::Solid2DSectionElementPropertyCard::SectionIntegratedPrestressAMatrix& f):
            MAST::SectionIntegratedPrestressMatrixBase(f),
            _prestress(f._prestress->clone().release()),
            _T(f._T->clone().release()),
            _h(f._h->clone().release())
            { }
            
            
            virtual ~SectionIntegratedPrestressAMatrix() {
                delete _prestress;
                delete _T;
                delete _h;
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>> clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>>
                (new MAST::Solid2DSectionElementPropertyCard::SectionIntegratedPrestressAMatrix(*this));
            }

            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const Point& p, const Real t, DenseMatrix<Real>& m) const;

            virtual void convert_to_vector(const DenseMatrix<Real>& m, DenseVector<Real>& v) const;

        protected:
            
            MAST::FieldFunction<DenseMatrix<Real> > *_prestress, *_T;
            MAST::FieldFunction<Real> *_h;
        };

        
        
        class SectionIntegratedPrestressBMatrix: public MAST::SectionIntegratedPrestressMatrixBase {
        public:
            SectionIntegratedPrestressBMatrix(MAST::FieldFunction<DenseMatrix<Real> > *prestress,
                                              MAST::FieldFunction<DenseMatrix<Real> > *T,
                                              MAST::FieldFunction<Real> *h);

            SectionIntegratedPrestressBMatrix(const MAST::Solid2DSectionElementPropertyCard::SectionIntegratedPrestressBMatrix& f):
            MAST::SectionIntegratedPrestressMatrixBase(f),
            _prestress(f._prestress->clone().release()),
            _T(f._T->clone().release()),
            _h(f._h->clone().release())
            { }
            

            virtual ~SectionIntegratedPrestressBMatrix() {
                delete _prestress;
                delete _T;
                delete _h;
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>> clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>>
                (new MAST::Solid2DSectionElementPropertyCard::SectionIntegratedPrestressBMatrix(*this));
            }

            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void convert_to_vector(const DenseMatrix<Real>& m, DenseVector<Real>& v) const;
            
        protected:
            
            MAST::FieldFunction<DenseMatrix<Real> > *_prestress, *_T;
            MAST::FieldFunction<Real> *_h;
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
         *   returns a function to evaluate the specified quantitys
         *   type \par t.
         */
        virtual std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>>
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
         *   material property card
         */
        MAST::MaterialPropertyCardBase *_material;
    };
    
}



#endif // __MAST_solid_2d_section_element_property_card_h__
