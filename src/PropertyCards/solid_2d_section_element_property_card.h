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
        
        
        class SectionIntegratedExtensionStiffnessMatrix: public MAST::FieldFunction<DenseRealMatrix > {
        public:
            SectionIntegratedExtensionStiffnessMatrix(MAST::FieldFunction<DenseRealMatrix > *mat,
                                                      MAST::FieldFunction<Real> *h);
            
            SectionIntegratedExtensionStiffnessMatrix(const MAST::Solid2DSectionElementPropertyCard::SectionIntegratedExtensionStiffnessMatrix& f):
            MAST::FieldFunction<DenseRealMatrix >(f),
            _material_stiffness(f._material_stiffness->clone().release()),
            _h(f._h->clone().release()) {
                _functions.insert(_material_stiffness);
                _functions.insert(_h);
            }

            
            virtual ~SectionIntegratedExtensionStiffnessMatrix() {
                delete _material_stiffness;
                delete _h;
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > > clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > >
                (new MAST::Solid2DSectionElementPropertyCard::SectionIntegratedExtensionStiffnessMatrix(*this));
            }

            virtual void operator() (const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
        protected:
            
            MAST::FieldFunction<DenseRealMatrix > *_material_stiffness;
            MAST::FieldFunction<Real> *_h;
        };
        
        
        
        class SectionIntegratedExtensionBendingStiffnessMatrix: public MAST::FieldFunction<DenseRealMatrix > {
        public:
            SectionIntegratedExtensionBendingStiffnessMatrix(MAST::FieldFunction<DenseRealMatrix > *mat,
                                                             MAST::FieldFunction<Real> *h,
                                                             MAST::FieldFunction<Real> *off);
            
            SectionIntegratedExtensionBendingStiffnessMatrix(const MAST::Solid2DSectionElementPropertyCard::SectionIntegratedExtensionBendingStiffnessMatrix& f):
            MAST::FieldFunction<DenseRealMatrix >(f),
            _material_stiffness(f._material_stiffness->clone().release()),
            _h(f._h->clone().release()),
            _off(f._off->clone().release()) {
                _functions.insert(_material_stiffness);
                _functions.insert(_h);
                _functions.insert(_off);
            }
            
            virtual ~SectionIntegratedExtensionBendingStiffnessMatrix() {
                delete _material_stiffness;
                delete _h;
                delete _off;
            }

            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > > clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > >
                (new MAST::Solid2DSectionElementPropertyCard::SectionIntegratedExtensionBendingStiffnessMatrix(*this));
            }
            
            virtual void operator() (const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
        protected:
            
            MAST::FieldFunction<DenseRealMatrix > *_material_stiffness;
            MAST::FieldFunction<Real> *_h, *_off;
        };
        
        
        class SectionIntegratedBendingStiffnessMatrix: public MAST::FieldFunction<DenseRealMatrix > {
        public:
            SectionIntegratedBendingStiffnessMatrix(MAST::FieldFunction<DenseRealMatrix > *mat,
                                                    MAST::FieldFunction<Real> *h,
                                                    MAST::FieldFunction<Real> *off);
            
            SectionIntegratedBendingStiffnessMatrix(const MAST::Solid2DSectionElementPropertyCard::SectionIntegratedBendingStiffnessMatrix& f):
            MAST::FieldFunction<DenseRealMatrix >(f),
            _material_stiffness(f._material_stiffness->clone().release()),
            _h(f._h->clone().release()),
            _off(f._off->clone().release()) {
                _functions.insert(_material_stiffness);
                _functions.insert(_h);
                _functions.insert(_off);
            }
            

            virtual ~SectionIntegratedBendingStiffnessMatrix() {
                delete _material_stiffness;
                delete _h;
                delete _off;
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > > clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > >
                (new MAST::Solid2DSectionElementPropertyCard::SectionIntegratedBendingStiffnessMatrix(*this));
            }

            virtual void operator() (const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
        protected:
            
            MAST::FieldFunction<DenseRealMatrix > *_material_stiffness;
            MAST::FieldFunction<Real> *_h, *_off;
        };
        
        
        
        
        class SectionIntegratedTransverseStiffnessMatrix: public MAST::FieldFunction<DenseRealMatrix > {
        public:
            SectionIntegratedTransverseStiffnessMatrix(MAST::FieldFunction<DenseRealMatrix > *mat,
                                                       MAST::FieldFunction<Real>* h):
            MAST::FieldFunction<DenseRealMatrix >("SectionIntegratedTransverseStiffnessMatrix2D"),
            _material_stiffness(mat),
            _h(h) {
                _functions.insert(mat);
                _functions.insert(h);
            }
            
            
            SectionIntegratedTransverseStiffnessMatrix(const MAST::Solid2DSectionElementPropertyCard::SectionIntegratedTransverseStiffnessMatrix &f):
            MAST::FieldFunction<DenseRealMatrix >(f),
            _material_stiffness(f._material_stiffness->clone().release()),
            _h(f._h->clone().release()) {
                _functions.insert(_material_stiffness);
                _functions.insert(_h);
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > > clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > >
                (new MAST::Solid2DSectionElementPropertyCard::SectionIntegratedTransverseStiffnessMatrix(*this));
            }
            
            virtual ~SectionIntegratedTransverseStiffnessMatrix() {
                delete _material_stiffness;
                delete _h;
            }
            
            virtual void operator() (const libMesh::Point& p, const Real t, DenseRealMatrix& m) const {
                Real h;
                (*_h)(p, t, h);
                (*_material_stiffness)(p, t, m);
                m.scale(h);
            }
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const libMesh::Point& p, const Real t, DenseRealMatrix& m) const {
                DenseRealMatrix dm;
                Real h, dh;
                (*_h)(p, t, h); _h->partial(f, p, t, dh);
                (*_material_stiffness)(p, t, m); _material_stiffness->partial(f, p, t, dm);
                
                m.scale(dh);
                m.add(h, dm);
            }
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const libMesh::Point& p, const Real t, DenseRealMatrix& m) const {
                DenseRealMatrix dm;
                Real h, dh;
                (*_h)(p, t, h); _h->total(f, p, t, dh);
                (*_material_stiffness)(p, t, m); _material_stiffness->total(f, p, t, dm);
                
                m.scale(dh);
                m.add(h, dm);
            }
            
        protected:
            
            MAST::FieldFunction<DenseRealMatrix > *_material_stiffness;
            MAST::FieldFunction<Real> *_h;
        };

        
        class SectionIntegratedInertiaMatrix: public MAST::FieldFunction<DenseRealMatrix > {
        public:
            SectionIntegratedInertiaMatrix(MAST::FieldFunction<Real> *rho,
                                           MAST::FieldFunction<Real> *h,
                                           MAST::FieldFunction<Real> *off);
            
            SectionIntegratedInertiaMatrix(const MAST::Solid2DSectionElementPropertyCard::SectionIntegratedInertiaMatrix& f):
            MAST::FieldFunction<DenseRealMatrix >(f),
            _rho(f._rho->clone().release()),
            _h(f._h->clone().release()),
            _off(f._off->clone().release()) {
                _functions.insert(_rho);
                _functions.insert(_h);
                _functions.insert(_off);
            }
            

            virtual ~SectionIntegratedInertiaMatrix() {
                delete _rho;
                delete _h;
                delete _off;
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > > clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > >
                (new MAST::Solid2DSectionElementPropertyCard::SectionIntegratedInertiaMatrix(*this));
            }

            virtual void operator() (const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
        protected:
            
            MAST::FieldFunction<Real> *_rho, *_h, *_off;
        };
        
        
        
        class SectionIntegratedThermalExpansionAMatrix: public MAST::FieldFunction<DenseRealMatrix > {
        public:
            SectionIntegratedThermalExpansionAMatrix(MAST::FieldFunction<DenseRealMatrix > *mat_stiff,
                                                     MAST::FieldFunction<DenseRealMatrix > *mat_expansion,
                                                     MAST::FieldFunction<Real> *h);
            
            SectionIntegratedThermalExpansionAMatrix(const MAST::Solid2DSectionElementPropertyCard::SectionIntegratedThermalExpansionAMatrix& f):
            MAST::FieldFunction<DenseRealMatrix >(f),
            _material_stiffness(f._material_stiffness->clone().release()),
            _material_expansion(f._material_expansion->clone().release()),
            _h(f._h->clone().release()) {
                _functions.insert(_material_stiffness);
                _functions.insert(_material_expansion);
                _functions.insert(_h);
            }
            
            
            
            virtual ~SectionIntegratedThermalExpansionAMatrix() {
                delete _material_stiffness;
                delete _material_expansion;
                delete _h;
            }
            
            
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > > clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > >
                (new MAST::Solid2DSectionElementPropertyCard::SectionIntegratedThermalExpansionAMatrix(*this));
            }
            
            
            virtual void operator() (const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
        protected:
            
            MAST::FieldFunction<DenseRealMatrix > *_material_stiffness;
            MAST::FieldFunction<DenseRealMatrix > *_material_expansion;
            MAST::FieldFunction<Real> *_h;
        };

        
        
        class SectionIntegratedThermalExpansionBMatrix: public MAST::FieldFunction<DenseRealMatrix > {
        public:
            SectionIntegratedThermalExpansionBMatrix(MAST::FieldFunction<DenseRealMatrix > *mat_stiff,
                                                    MAST::FieldFunction<DenseRealMatrix > *mat_expansion,
                                                     MAST::FieldFunction<Real> *h,
                                                     MAST::FieldFunction<Real> *off);
            
            SectionIntegratedThermalExpansionBMatrix(const MAST::Solid2DSectionElementPropertyCard::SectionIntegratedThermalExpansionBMatrix& f):
            MAST::FieldFunction<DenseRealMatrix >(f),
            _material_stiffness(f._material_stiffness->clone().release()),
            _material_expansion(f._material_expansion->clone().release()),
            _h(f._h->clone().release()),
            _off(f._off->clone().release()) {
                _functions.insert(_material_stiffness);
                _functions.insert(_material_expansion);
                _functions.insert(_h);
                _functions.insert(_off);
            }
            
            
            virtual ~SectionIntegratedThermalExpansionBMatrix() {
                delete _material_stiffness;
                delete _material_expansion;
                delete _h;
                delete _off;
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > > clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > >
                (new MAST::Solid2DSectionElementPropertyCard::SectionIntegratedThermalExpansionBMatrix(*this));
            }
            
            
            virtual void operator() (const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
            
        protected:
            
            MAST::FieldFunction<DenseRealMatrix > *_material_stiffness;
            MAST::FieldFunction<DenseRealMatrix > *_material_expansion;
            MAST::FieldFunction<Real> *_h, *_off;
        };

        
        
        
        class SectionIntegratedPrestressAMatrix: public MAST::SectionIntegratedPrestressMatrixBase {
        public:
            SectionIntegratedPrestressAMatrix(MAST::FieldFunction<DenseRealMatrix > *prestress,
                                              MAST::FieldFunction<DenseRealMatrix > *T,
                                              MAST::FieldFunction<Real> *h);

            SectionIntegratedPrestressAMatrix(const MAST::Solid2DSectionElementPropertyCard::SectionIntegratedPrestressAMatrix& f):
            MAST::SectionIntegratedPrestressMatrixBase(f),
            _prestress(f._prestress->clone().release()),
            _T(f._T->clone().release()),
            _h(f._h->clone().release()) {
                _functions.insert(_prestress);
                _functions.insert(_T);
                _functions.insert(_h);
            }
            
            
            virtual ~SectionIntegratedPrestressAMatrix() {
                delete _prestress;
                delete _T;
                delete _h;
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > > clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > >
                (new MAST::Solid2DSectionElementPropertyCard::SectionIntegratedPrestressAMatrix(*this));
            }

            virtual void operator() (const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;

            virtual void convert_to_vector(const DenseRealMatrix& m, DenseRealVector& v) const;

        protected:
            
            MAST::FieldFunction<DenseRealMatrix > *_prestress, *_T;
            MAST::FieldFunction<Real> *_h;
        };

        
        
        class SectionIntegratedPrestressBMatrix: public MAST::SectionIntegratedPrestressMatrixBase {
        public:
            SectionIntegratedPrestressBMatrix(MAST::FieldFunction<DenseRealMatrix > *prestress,
                                              MAST::FieldFunction<DenseRealMatrix > *T,
                                              MAST::FieldFunction<Real> *h,
                                              MAST::FieldFunction<Real> *off);

            SectionIntegratedPrestressBMatrix(const MAST::Solid2DSectionElementPropertyCard::SectionIntegratedPrestressBMatrix& f):
            MAST::SectionIntegratedPrestressMatrixBase(f),
            _prestress(f._prestress->clone().release()),
            _T(f._T->clone().release()),
            _h(f._h->clone().release()),
            _off(f._off->clone().release()) {
                _functions.insert(_prestress);
                _functions.insert(_T);
                _functions.insert(_h);
                _functions.insert(_off);
            }
            

            virtual ~SectionIntegratedPrestressBMatrix() {
                delete _prestress;
                delete _T;
                delete _h;
                delete _off;
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > > clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > >
                (new MAST::Solid2DSectionElementPropertyCard::SectionIntegratedPrestressBMatrix(*this));
            }

            virtual void operator() (const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
            virtual void convert_to_vector(const DenseRealMatrix& m, DenseRealVector& v) const;
            
        protected:
            
            MAST::FieldFunction<DenseRealMatrix > *_prestress, *_T;
            MAST::FieldFunction<Real> *_h, *_off;
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
        virtual std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > >
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
