//
//  solid_1d_section_element_property_card.h
//  MAST
//
//  Created by Manav Bhatia on 1/30/14.
//  Copyright (c) 2014 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_solid_1d_section_element_property_card_h__
#define __MAST_solid_1d_section_element_property_card_h__


// MAST includes
#include "PropertyCards/element_property_card_1D.h"


namespace MAST {
    class Solid1DSectionElementPropertyCard : public MAST::ElementPropertyCard1D {
        
    public:
        
        Solid1DSectionElementPropertyCard(unsigned int pid):
        MAST::ElementPropertyCard1D(pid),
        _material(NULL)
        { }
        
        
        /*!
         *   virtual destructor
         */
        virtual ~Solid1DSectionElementPropertyCard() { }
        
        
        
        class Area: public MAST::FieldFunction<Real> {
        public:
            Area(MAST::FieldFunction<Real> *hy,
                 MAST::FieldFunction<Real>* hz):
            MAST::FieldFunction<Real>("Area"),
            _hy(hy),
            _hz(hz) {
                _functions.insert(hy);
                _functions.insert(hz);
            }
            
            Area(const MAST::Solid1DSectionElementPropertyCard::Area &f):
            MAST::FieldFunction<Real>(f),
            _hy(f._hy->clone().release()),
            _hz(f._hz->clone().release())
            { }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<Real> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<Real> >
                (new MAST::Solid1DSectionElementPropertyCard::Area(*this));
            }
            
            virtual ~Area() {
                delete _hy;
                delete _hz;
            }
            
            virtual void operator() (const Point& p, const Real t, Real& m) const {
                Real hy, hz;
                (*_hy)(p, t, hy);
                (*_hz)(p, t, hz);
                
                m = hy*hz;
            }
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const Point& p, const Real t, Real& m) const {
                Real hy, hz, dhy, dhz;
                (*_hy)(p, t, hy); _hy->partial(f, p, t, dhy);
                (*_hz)(p, t, hz); _hz->partial(f, p, t, dhz);
                
                m = dhy*hz + hy*dhz;
            }
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const Point& p, const Real t, Real& m) const {
                Real hy, hz, dhy, dhz;
                (*_hy)(p, t, hy); _hy->total(f, p, t, dhy);
                (*_hz)(p, t, hz); _hz->total(f, p, t, dhz);
                
                m = dhy*hz + hy*dhz;
            }
            
        protected:
            
            MAST::FieldFunction<Real> *_hy, *_hz;
        };
        
        
        
        class TorsionalConstant: public MAST::FieldFunction<Real> {
        public:
            TorsionalConstant(MAST::FieldFunction<Real> *hy,
                              MAST::FieldFunction<Real>* hz):
            MAST::FieldFunction<Real>("TorsionalConstant"),
            _hy(hy),
            _hz(hz) {
                _functions.insert(hy);
                _functions.insert(hz);
            }
            
            TorsionalConstant(const MAST::Solid1DSectionElementPropertyCard::TorsionalConstant &f):
            MAST::FieldFunction<Real>(f),
            _hy(f._hy->clone().release()),
            _hz(f._hz->clone().release())
            { }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<Real> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<Real> >
                (new MAST::Solid1DSectionElementPropertyCard::TorsionalConstant(*this));
            }
            
            virtual ~TorsionalConstant() {
                delete _hy;
                delete _hz;
            }
            
            virtual void operator() (const Point& p, const Real t, Real& m) const {
                Real hy, hz;
                (*_hy)(p, t, hy);
                (*_hz)(p, t, hz);
                
                m = hy*hz;
            }
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const Point& p, const Real t, Real& m) const {
                Real hy, hz, dhy, dhz;
                (*_hy)(p, t, hy); _hy->partial(f, p, t, dhy);
                (*_hz)(p, t, hz); _hz->partial(f, p, t, dhz);
                
                m = dhy*hz + hy*dhz;
            }
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const Point& p, const Real t, Real& m) const {
                Real hy, hz, dhy, dhz;
                (*_hy)(p, t, hy); _hy->total(f, p, t, dhy);
                (*_hz)(p, t, hz); _hz->total(f, p, t, dhz);
                
                m = dhy*hz + hy*dhz;
            }
            
        protected:
            
            MAST::FieldFunction<Real> *_hy, *_hz;
        };

        
        
        class PolarInertia: public MAST::FieldFunction<Real> {
        public:
            PolarInertia(MAST::FieldFunction<Real> *hy,
                              MAST::FieldFunction<Real>* hz):
            MAST::FieldFunction<Real>("PolarInertia"),
            _hy(hy),
            _hz(hz) {
                _functions.insert(hy);
                _functions.insert(hz);
            }
            
            PolarInertia(const MAST::Solid1DSectionElementPropertyCard::PolarInertia &f):
            MAST::FieldFunction<Real>(f),
            _hy(f._hy->clone().release()),
            _hz(f._hz->clone().release())
            { }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<Real> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<Real> >
                (new MAST::Solid1DSectionElementPropertyCard::PolarInertia(*this));
            }
            
            virtual ~PolarInertia() {
                delete _hy;
                delete _hz;
            }
            
            virtual void operator() (const Point& p, const Real t, Real& m) const {
                Real hy, hz;
                (*_hy)(p, t, hy);
                (*_hz)(p, t, hz);
                
                m = hy*hz;
            }
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const Point& p, const Real t, Real& m) const {
                Real hy, hz, dhy, dhz;
                (*_hy)(p, t, hy); _hy->partial(f, p, t, dhy);
                (*_hz)(p, t, hz); _hz->partial(f, p, t, dhz);
                
                m = dhy*hz + hy*dhz;
            }
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const Point& p, const Real t, Real& m) const {
                Real hy, hz, dhy, dhz;
                (*_hy)(p, t, hy); _hy->total(f, p, t, dhy);
                (*_hz)(p, t, hz); _hz->total(f, p, t, dhz);
                
                m = dhy*hz + hy*dhz;
            }
            
        protected:
            
            MAST::FieldFunction<Real> *_hy, *_hz;
        };

        
        
        
        class AreaYMoment: public MAST::FieldFunction<Real> {
        public:
            AreaYMoment(MAST::FieldFunction<Real>* hy,
                        MAST::FieldFunction<Real>* hz,
                        MAST::FieldFunction<Real>* hz_offset):
            MAST::FieldFunction<Real>("AreaYMoment"),
            _hy(hy),
            _hz(hz),
            _hz_offset(hz_offset) {
                _functions.insert(hy);
                _functions.insert(hz);
                _functions.insert(hz_offset);
            }
            
            AreaYMoment(const MAST::Solid1DSectionElementPropertyCard::AreaYMoment &f):
            MAST::FieldFunction<Real>(f),
            _hy(f._hy->clone().release()),
            _hz(f._hz->clone().release()),
            _hz_offset(f._hz_offset->clone().release())
            { }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<Real> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<Real> >
                (new MAST::Solid1DSectionElementPropertyCard::AreaYMoment(*this));
            }
            
            virtual ~AreaYMoment() {
                delete _hy;
                delete _hz;
                delete _hz_offset;
            }
            
            virtual void operator() (const Point& p, const Real t, Real& m) const {
                Real hy, hz, off;
                (*_hy)(p, t, hy);
                (*_hz)(p, t, hz);
                (*_hz_offset)(p, t, off);
                
                m = hy*hz*off;
            }
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const Point& p, const Real t, Real& m) const {
                Real hy, hz, off, dhy, dhz, doff;
                (*_hy)(p, t, hy); _hy->partial(f, p, t, dhy);
                (*_hz)(p, t, hz); _hz->partial(f, p, t, dhz);
                (*_hz_offset)(p, t, off); _hz_offset->partial(f, p, t, doff);
                
                m = dhy*hz*off + hy*dhz*off + hy*hz*doff;
            }
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const Point& p, const Real t, Real& m) const {
                Real hy, hz, off, dhy, dhz, doff;
                (*_hy)(p, t, hy); _hy->total(f, p, t, dhy);
                (*_hz)(p, t, hz); _hz->total(f, p, t, dhz);
                (*_hz_offset)(p, t, off); _hz_offset->total(f, p, t, doff);
                
                m = dhy*hz*off + hy*dhz*off + hy*hz*doff;
            }
            
        protected:
            
            MAST::FieldFunction<Real> *_hy, *_hz, *_hz_offset;
        };
        
        
        
        class AreaZMoment: public MAST::FieldFunction<Real> {
        public:
            AreaZMoment(MAST::FieldFunction<Real>* hy,
                        MAST::FieldFunction<Real>* hz,
                        MAST::FieldFunction<Real>* hy_offset):
            MAST::FieldFunction<Real>("AreaZMoment"),
            _hy(hy),
            _hz(hz),
            _hy_offset(hy_offset) {
                _functions.insert(hy);
                _functions.insert(hz);
                _functions.insert(hy_offset);
            }
            
            AreaZMoment(const MAST::Solid1DSectionElementPropertyCard::AreaZMoment &f):
            MAST::FieldFunction<Real>(f),
            _hy(f._hy->clone().release()),
            _hz(f._hz->clone().release()),
            _hy_offset(f._hy_offset->clone().release())
            { }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<Real> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<Real> >
                (new MAST::Solid1DSectionElementPropertyCard::AreaZMoment(*this));
            }
            
            virtual ~AreaZMoment() {
                delete _hy;
                delete _hz;
                delete _hy_offset;
            }
            
            virtual void operator() (const Point& p, const Real t, Real& m) const {
                Real hy, hz, off;
                (*_hy)(p, t, hy);
                (*_hz)(p, t, hz);
                (*_hy_offset)(p, t, off);
                
                m = hy*hz*off;
            }
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const Point& p, const Real t, Real& m) const {
                Real hy, hz, off, dhy, dhz, doff;
                (*_hy)(p, t, hy); _hy->partial(f, p, t, dhy);
                (*_hz)(p, t, hz); _hz->partial(f, p, t, dhz);
                (*_hy_offset)(p, t, off); _hy_offset->partial(f, p, t, doff);
                
                m = dhy*hz*off + hy*dhz*off + hy*hz*doff;
            }
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const Point& p, const Real t, Real& m) const {
                Real hy, hz, off, dhy, dhz, doff;
                (*_hy)(p, t, hy); _hy->total(f, p, t, dhy);
                (*_hz)(p, t, hz); _hz->total(f, p, t, dhz);
                (*_hy_offset)(p, t, off); _hy_offset->total(f, p, t, doff);
                
                m = dhy*hz*off + hy*dhz*off + hy*hz*doff;
            }
            
        protected:
            
            MAST::FieldFunction<Real> *_hy, *_hz, *_hy_offset;
        };
        
        
        
        
        class AreaInertiaMatrix: public MAST::FieldFunction<DenseMatrix<Real> > {
        public:
            AreaInertiaMatrix(MAST::FieldFunction<Real>* hy,
                              MAST::FieldFunction<Real>* hz,
                              MAST::FieldFunction<Real>* hy_offset,
                              MAST::FieldFunction<Real>* hz_offset):
            MAST::FieldFunction<DenseMatrix<Real> >("AreaInertiaMatrix"),
            _hy(hy),
            _hz(hz),
            _hy_offset(hy_offset),
            _hz_offset(hz_offset) {
                _functions.insert(hy);
                _functions.insert(hz);
                _functions.insert(hy_offset);
                _functions.insert(hz_offset);
            }
            
            AreaInertiaMatrix(const MAST::Solid1DSectionElementPropertyCard::AreaInertiaMatrix &f):
            MAST::FieldFunction<DenseMatrix<Real> >(f),
            _hy(f._hy->clone().release()),
            _hz(f._hz->clone().release()),
            _hy_offset(f._hy_offset->clone().release()),
            _hz_offset(f._hz_offset->clone().release())
            { }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>> clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>>
                (new MAST::Solid1DSectionElementPropertyCard::AreaInertiaMatrix(*this));
            }
            
            virtual ~AreaInertiaMatrix() {
                delete _hy;
                delete _hz;
                delete _hy_offset;
                delete _hz_offset;
            }
            
            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& m) const {
                Real hy, hz, offy, offz;
                m.resize(2,2);
                (*_hy)(p, t, hy);
                (*_hz)(p, t, hz);
                (*_hy_offset)(p, t, offy);
                (*_hz_offset)(p, t, offz);
                
                m(0,0) = hz*pow(hy,3)/12. + hy*hz*pow(offy,2) ; // Izz for v-bending
                m(0,1) = hy*hz*offy*offz;
                m(1,0) = m(0,1);
                m(0,0) = hy*pow(hz,3)/12. + hy*hz*pow(offz,2) ; // Iyy for w-bending
            }
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const Point& p, const Real t, DenseMatrix<Real>& m) const {
                Real hy, hz, offy, offz, dhy, dhz, doffy, doffz;
                m.resize(2,2);
                (*_hy)(p, t, hy); _hy->partial(f, p, t, dhy);
                (*_hz)(p, t, hz); _hz->partial(f, p, t, dhz);
                (*_hy_offset)(p, t, offy); _hy_offset->partial(f, p, t, doffy);
                (*_hz_offset)(p, t, offz); _hz_offset->partial(f, p, t, doffz);
                
                
                m(0,0) = dhz*pow(hy,3)/12. + hz*pow(hy,2)/4.*dhy +
                dhy*hz*pow(offy,2) + hy*dhz*pow(offy,2) + 2.*hy*hz*offy*doffy ;
                m(0,1) = dhy*hz*offy*offz + hy*dhz*offy*offz +
                hy*hz*doffy*offz + hy*hz*offy*doffz;
                m(1,0) = m(0,1);
                m(0,0) = dhy*pow(hz,3)/12. + hy*pow(hz,2)/4.*dhz +
                dhy*hz*pow(offz,2) + hy*dhz*pow(offz,2) + 2.*hy*hz*offz*doffz ;
            }
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const Point& p, const Real t, DenseMatrix<Real>& m) const {
                Real hy, hz, offy, offz, dhy, dhz, doffy, doffz;
                m.resize(2,2);
                (*_hy)(p, t, hy); _hy->total(f, p, t, dhy);
                (*_hz)(p, t, hz); _hz->total(f, p, t, dhz);
                (*_hy_offset)(p, t, offy); _hy_offset->total(f, p, t, doffy);
                (*_hz_offset)(p, t, offz); _hz_offset->total(f, p, t, doffz);
                
                
                m(0,0) = dhz*pow(hy,3)/12. + hz*pow(hy,2)/4.*dhy +
                dhy*hz*pow(offy,2) + hy*dhz*pow(offy,2) + 2.*hy*hz*offy*doffy ;
                m(0,1) = dhy*hz*offy*offz + hy*dhz*offy*offz +
                hy*hz*doffy*offz + hy*hz*offy*doffz;
                m(1,0) = m(0,1);
                m(0,0) = dhy*pow(hz,3)/12. + hy*pow(hz,2)/4.*dhz +
                dhy*hz*pow(offz,2) + hy*dhz*pow(offz,2) + 2.*hy*hz*offz*doffz ;
            }
            
        protected:
            
            MAST::FieldFunction<Real> *_hy, *_hz, *_hy_offset, *_hz_offset;
        };
        
        
        class SectionIntegratedExtensionStiffnessMatrix: public MAST::FieldFunction<DenseMatrix<Real> > {
        public:
            SectionIntegratedExtensionStiffnessMatrix(MAST::FieldFunction<DenseMatrix<Real> > *mat,
                                                      MAST::FieldFunction<Real>* A,
                                                      MAST::FieldFunction<Real>* J);
            
            SectionIntegratedExtensionStiffnessMatrix(const MAST::Solid1DSectionElementPropertyCard::SectionIntegratedExtensionStiffnessMatrix &f):
            MAST::FieldFunction<DenseMatrix<Real> >(f),
            _material_stiffness(f._material_stiffness->clone().release()),
            _A(f._A->clone().release()),
            _J(f._J->clone().release())
            { }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>> clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>>
                (new MAST::Solid1DSectionElementPropertyCard::SectionIntegratedExtensionStiffnessMatrix(*this));
            }
            
            virtual ~SectionIntegratedExtensionStiffnessMatrix() {
                delete _material_stiffness;
                delete _A;
                delete _J;
            }
            
            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
        protected:
            
            MAST::FieldFunction<DenseMatrix<Real> > *_material_stiffness;
            MAST::FieldFunction<Real> *_A, *_J;
        };
        
        
        
        class SectionIntegratedExtensionBendingStiffnessMatrix: public MAST::FieldFunction<DenseMatrix<Real> > {
        public:
            SectionIntegratedExtensionBendingStiffnessMatrix(MAST::FieldFunction<DenseMatrix<Real> > *mat,
                                                             MAST::FieldFunction<Real>* A_y_moment,
                                                             MAST::FieldFunction<Real>* A_z_moment);
            
            SectionIntegratedExtensionBendingStiffnessMatrix(const MAST::Solid1DSectionElementPropertyCard::SectionIntegratedExtensionBendingStiffnessMatrix &f):
            MAST::FieldFunction<DenseMatrix<Real> >(f),
            _material_stiffness(f._material_stiffness->clone().release()),
            _A_y_moment(f._A_y_moment->clone().release()),
            _A_z_moment(f._A_z_moment->clone().release())
            { }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>> clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>>
                (new MAST::Solid1DSectionElementPropertyCard::SectionIntegratedExtensionBendingStiffnessMatrix(*this));
            }
            
            virtual ~SectionIntegratedExtensionBendingStiffnessMatrix() {
                delete _material_stiffness;
                delete _A_y_moment;
                delete _A_z_moment;
            }
            
            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
        protected:
            
            MAST::FieldFunction<DenseMatrix<Real> > *_material_stiffness;
            MAST::FieldFunction<Real> *_A_y_moment, *_A_z_moment;
        };
        
        
        class SectionIntegratedBendingStiffnessMatrix: public MAST::FieldFunction<DenseMatrix<Real> > {
        public:
            SectionIntegratedBendingStiffnessMatrix(MAST::FieldFunction<DenseMatrix<Real> > *mat,
                                                    MAST::FieldFunction<DenseMatrix<Real> > *I);
            
            SectionIntegratedBendingStiffnessMatrix(const MAST::Solid1DSectionElementPropertyCard::SectionIntegratedBendingStiffnessMatrix &f):
            MAST::FieldFunction<DenseMatrix<Real> >(f),
            _material_stiffness(f._material_stiffness->clone().release()),
            _I(f._I->clone().release())
            { }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>> clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>>
                (new MAST::Solid1DSectionElementPropertyCard::SectionIntegratedBendingStiffnessMatrix(*this));
            }
            
            virtual ~SectionIntegratedBendingStiffnessMatrix() {
                delete _material_stiffness;
                delete _I;
            }
            
            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
        protected:
            
            MAST::FieldFunction<DenseMatrix<Real> > *_material_stiffness;
            MAST::FieldFunction<DenseMatrix<Real> > *_I;
        };
        
        
        class SectionIntegratedTransverseStiffnessMatrix: public MAST::FieldFunction<DenseMatrix<Real> > {
        public:
            SectionIntegratedTransverseStiffnessMatrix(MAST::FieldFunction<DenseMatrix<Real> > *mat,
                                                       MAST::FieldFunction<Real>* A):
            MAST::FieldFunction<DenseMatrix<Real> >("SectionIntegratedTransverseStiffnessMatrix1D"),
            _material_stiffness(mat),
            _A(A) {
                _functions.insert(mat);
                _functions.insert(A);
            }
            
            
            SectionIntegratedTransverseStiffnessMatrix(const MAST::Solid1DSectionElementPropertyCard::SectionIntegratedTransverseStiffnessMatrix &f):
            MAST::FieldFunction<DenseMatrix<Real> >(f),
            _material_stiffness(f._material_stiffness->clone().release()),
            _A(f._A->clone().release())
            { }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>> clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>>
                (new MAST::Solid1DSectionElementPropertyCard::SectionIntegratedTransverseStiffnessMatrix(*this));
            }
            
            virtual ~SectionIntegratedTransverseStiffnessMatrix() {
                delete _material_stiffness;
                delete _A;
            }
            
            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& m) const {
                Real A;
                (*_A)(p, t, A);
                (*_material_stiffness)(p, t, m);
                m.scale(A);
            }
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const Point& p, const Real t, DenseMatrix<Real>& m) const {
                DenseMatrix<Real> dm;
                Real A, dA;
                (*_A)(p, t, A); _A->partial(f, p, t, dA);
                (*_material_stiffness)(p, t, m); _material_stiffness->partial(f, p, t, dm);
                
                m.scale(dA);
                m.add(A, dm);
            }
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const Point& p, const Real t, DenseMatrix<Real>& m) const {
                DenseMatrix<Real> dm;
                Real A, dA;
                (*_A)(p, t, A); _A->total(f, p, t, dA);
                (*_material_stiffness)(p, t, m); _material_stiffness->total(f, p, t, dm);
                
                m.scale(dA);
                m.add(A, dm);
            }
            
        protected:
            
            MAST::FieldFunction<DenseMatrix<Real> > *_material_stiffness;
            MAST::FieldFunction<Real> *_A;
        };

        
        class SectionIntegratedInertiaMatrix: public MAST::FieldFunction<DenseMatrix<Real> > {
        public:
            SectionIntegratedInertiaMatrix(MAST::FieldFunction<Real>* rho,
                                           MAST::FieldFunction<Real>* A,
                                           MAST::FieldFunction<Real>* A_y_moment,
                                           MAST::FieldFunction<Real>* A_z_moment,
                                           MAST::FieldFunction<Real>* Ip,
                                           MAST::FieldFunction<DenseMatrix<Real> >* I);
            
            SectionIntegratedInertiaMatrix(const MAST::Solid1DSectionElementPropertyCard::SectionIntegratedInertiaMatrix &f):
            MAST::FieldFunction<DenseMatrix<Real> >(f),
            _rho(f._rho->clone().release()),
            _A(f._A->clone().release()),
            _A_y_moment(f._A_y_moment->clone().release()),
            _A_z_moment(f._A_z_moment->clone().release()),
            _Ip(f._Ip->clone().release()),
            _I(f._I->clone().release())
            { }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>> clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>>
                (new MAST::Solid1DSectionElementPropertyCard::SectionIntegratedInertiaMatrix(*this));
            }
            
            virtual ~SectionIntegratedInertiaMatrix() {
                delete _rho;
                delete _A;
                delete _A_y_moment;
                delete _A_z_moment;
                delete _Ip;
                delete _I;
            }
            
            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
        protected:
            
            MAST::FieldFunction<Real> *_rho, *_A, *_A_y_moment, *_A_z_moment, *_Ip;
            MAST::FieldFunction<DenseMatrix<Real> > *_I;
        };
        
        
        
        class SectionIntegratedThermalExpansionMatrix: public MAST::FieldFunction<DenseMatrix<Real> > {
        public:
            SectionIntegratedThermalExpansionMatrix(MAST::FieldFunction<DenseMatrix<Real> > *mat_stiff,
                                                    MAST::FieldFunction<DenseMatrix<Real> > *mat_expansion,
                                                    MAST::FieldFunction<Real> *A,
                                                    MAST::FieldFunction<Real> *A_y_moment,
                                                    MAST::FieldFunction<Real> *A_z_moment);
            
            SectionIntegratedThermalExpansionMatrix(const MAST::Solid1DSectionElementPropertyCard::SectionIntegratedThermalExpansionMatrix &f):
            MAST::FieldFunction<DenseMatrix<Real> >(f),
            _material_stiffness(f._material_stiffness->clone().release()),
            _material_expansion(f._material_expansion->clone().release()),
            _A(f._A->clone().release()),
            _A_y_moment(f._A_y_moment->clone().release()),
            _A_z_moment(f._A_z_moment->clone().release())
            { }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>> clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>>
                (new MAST::Solid1DSectionElementPropertyCard::SectionIntegratedThermalExpansionMatrix(*this));
            }
            
            virtual ~SectionIntegratedThermalExpansionMatrix() {
                delete _material_stiffness;
                delete _material_expansion;
                delete _A;
                delete _A_y_moment;
                delete _A_z_moment;
            }
            
            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
        protected:
            
            MAST::FieldFunction<DenseMatrix<Real> > *_material_stiffness;
            MAST::FieldFunction<DenseMatrix<Real> > *_material_expansion;
            MAST::FieldFunction<Real> *_A, *_A_y_moment, *_A_z_moment;
        };
        
        
        
        
        class SectionIntegratedPrestressAMatrix: public MAST::SectionIntegratedPrestressMatrixBase {
        public:
            SectionIntegratedPrestressAMatrix(MAST::FieldFunction<DenseMatrix<Real> > *prestress,
                                              MAST::FieldFunction<DenseMatrix<Real> > *T,
                                              MAST::FieldFunction<Real> *A);
            
            SectionIntegratedPrestressAMatrix(const MAST::Solid1DSectionElementPropertyCard::SectionIntegratedPrestressAMatrix &f):
            MAST::SectionIntegratedPrestressMatrixBase(f),
            _prestress(f._prestress->clone().release()),
            _T(f._T->clone().release()),
            _A(f._A->clone().release())
            { }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>> clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>>
                (new MAST::Solid1DSectionElementPropertyCard::SectionIntegratedPrestressAMatrix(*this));
            }
            
            virtual ~SectionIntegratedPrestressAMatrix() {
                delete _prestress;
                delete _T;
                delete _A;
            }
            
            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void convert_to_vector(const DenseMatrix<Real>& m, DenseVector<Real>& v) const;
            
        protected:
            
            MAST::FieldFunction<DenseMatrix<Real> > *_prestress, *_T;
            MAST::FieldFunction<Real> *_A;
        };
        
        
        
        class SectionIntegratedPrestressBMatrix: public MAST::SectionIntegratedPrestressMatrixBase {
        public:
            SectionIntegratedPrestressBMatrix(MAST::FieldFunction<DenseMatrix<Real> > *prestress,
                                              MAST::FieldFunction<DenseMatrix<Real> > *T,
                                              MAST::FieldFunction<Real> *A_y_moment,
                                              MAST::FieldFunction<Real> *A_z_moment);
            
            SectionIntegratedPrestressBMatrix(const MAST::Solid1DSectionElementPropertyCard::SectionIntegratedPrestressBMatrix &f):
            MAST::SectionIntegratedPrestressMatrixBase(f),
            _prestress(f._prestress->clone().release()),
            _T(f._T->clone().release()),
            _A_y_moment(f._A_y_moment->clone().release()),
            _A_z_moment(f._A_z_moment->clone().release())
            { }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>> clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>>
                (new MAST::Solid1DSectionElementPropertyCard::SectionIntegratedPrestressBMatrix(*this));
            }
            
            virtual ~SectionIntegratedPrestressBMatrix() {
                delete _prestress;
                delete _T;
                delete _A_y_moment;
                delete _A_z_moment;
            }
            
            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void convert_to_vector(const DenseMatrix<Real>& m, DenseVector<Real>& v) const;
            
        protected:
            
            MAST::FieldFunction<DenseMatrix<Real> > *_prestress, *_T;
            MAST::FieldFunction<Real> *_A_y_moment, *_A_z_moment;
        };
        
        
        /*!
         *    sets the material card
         */
        void set_material(MAST::MaterialPropertyCardBase& mat) {
            _material = &mat;
        }
        
        
        /*!
         *   return true if the property is isotropic
         */
        virtual bool if_isotropic() const {
            return true;
        }
        
        
        /*!
         *    returns a reference to the material
         */
        virtual const MAST::MaterialPropertyCardBase& get_material() const {
            libmesh_assert(_material); // make sure it has already been set
            return *_material;
        }
        
        /*!
         *   returns value of the property \par val. The string values for
         *   \par val are IYY, IZZ, IYZ
         */
        virtual Real value(const std::string& val) const {libmesh_error();}
        
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



#endif // __MAST_solid_1d_section_element_property_card_h__
