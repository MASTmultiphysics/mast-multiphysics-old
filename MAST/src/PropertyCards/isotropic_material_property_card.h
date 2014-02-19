//
//  isotropic_material_property_card.h
//  MAST
//
//  Created by Manav Bhatia on 1/30/14.
//  Copyright (c) 2014 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_isotropic_material_property_card_h__
#define __MAST_isotropic_material_property_card_h__


// C++ functional
#include <functional>

// MAST includes
#include "PropertyCards/material_property_card_base.h"


// libMesh includes
#include "libmesh/dense_matrix.h"


namespace MAST {
    
    // Forward declerations
    class ElementPropertyCardBase;
    
    class IsotropicMaterialPropertyCard: public MAST::MaterialPropertyCardBase {
        
    public:
        
        IsotropicMaterialPropertyCard(unsigned int pid):
        MAST::MaterialPropertyCardBase (pid)
        { }
        
        
        class StiffnessMatrix1D: public MAST::FieldFunction<DenseMatrix<Real> > {
        public:
            
            StiffnessMatrix1D( MAST::FieldFunction<Real>* E,
                              MAST::FieldFunction<Real>* nu);

            StiffnessMatrix1D(const MAST::IsotropicMaterialPropertyCard::StiffnessMatrix1D& f):
            MAST::FieldFunction<DenseMatrix<Real> >(f),
            _E(f._E->clone().release()),
            _nu(f._nu->clone().release()),
            _m(new std::function<void(const Real, const Real, DenseMatrix<Real>&)>(*f._m)),
            _parm_parE(new std::function<void(const Real, const Real, DenseMatrix<Real>&)>(*f._parm_parE)),
            _parm_parnu(new std::function<void(const Real, const Real, DenseMatrix<Real>&)>(*f._parm_parnu)){
                _functions.insert(_E);
                _functions.insert(_nu);
            }

            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>> clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>>
                (new MAST::IsotropicMaterialPropertyCard::StiffnessMatrix1D(*this));
            }

            virtual ~StiffnessMatrix1D();
            
            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            
        protected:
            
            MAST::FieldFunction<Real>* _E;
            MAST::FieldFunction<Real>* _nu;
            std::function<void(const Real, const Real, DenseMatrix<Real>&)>
            *_m, *_parm_parE, *_parm_parnu;
        };
        
        
        
        class TransverseShearStiffnessMatrix: public MAST::FieldFunction<DenseMatrix<Real> > {
        public:
            TransverseShearStiffnessMatrix( MAST::FieldFunction<Real>* E,
                                           MAST::FieldFunction<Real>* nu,
                                           MAST::FieldFunction<Real>* kappa);

            TransverseShearStiffnessMatrix(const MAST::IsotropicMaterialPropertyCard::TransverseShearStiffnessMatrix& f):
            MAST::FieldFunction<DenseMatrix<Real> >(f),
            _E(f._E->clone().release()),
            _nu(f._nu->clone().release()),
            _m(new std::function<void(const Real, const Real, const Real, DenseMatrix<Real>&)>(*f._m)),
            _parm_parE(new std::function<void(const Real, const Real, const Real, DenseMatrix<Real>&)>(*f._parm_parE)),
            _parm_parnu(new std::function<void(const Real, const Real, const Real, DenseMatrix<Real>&)>(*f._parm_parnu)),
            _parm_parkappa(new std::function<void(const Real, const Real, const Real, DenseMatrix<Real>&)>(*f._parm_parkappa)) {
                _functions.insert(_E);
                _functions.insert(_nu);
            }

            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>> clone() const  {
                return std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>>
                (new MAST::IsotropicMaterialPropertyCard::TransverseShearStiffnessMatrix(*this));
            }

            virtual ~TransverseShearStiffnessMatrix();
            
            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            
        protected:
            
            MAST::FieldFunction<Real>* _E;
            MAST::FieldFunction<Real>* _nu;
            MAST::FieldFunction<Real>* _kappa;
            std::function<void(const Real, const Real, const Real, DenseMatrix<Real>&)>
            *_m, *_parm_parE, *_parm_parnu, *_parm_parkappa;
        };
        
        
        class StiffnessMatrix2D: public MAST::FieldFunction<DenseMatrix<Real> > {
        public:
            StiffnessMatrix2D(MAST::FieldFunction<Real>* E,
                              MAST::FieldFunction<Real>* nu,
                              bool plane_stress);

            StiffnessMatrix2D(const MAST::IsotropicMaterialPropertyCard::StiffnessMatrix2D& f):
            MAST::FieldFunction<DenseMatrix<Real> >(f),
            _E(f._E->clone().release()),
            _nu(f._nu->clone().release()),
            _m(new std::function<void(const Real, const Real, DenseMatrix<Real>&)>(*f._m)),
            _parm_parE(new std::function<void(const Real, const Real, DenseMatrix<Real>&)>(*f._parm_parE)),
            _parm_parnu(new std::function<void(const Real, const Real, DenseMatrix<Real>&)>(*f._parm_parnu)) {
                _functions.insert(_E);
                _functions.insert(_nu);
            }

            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>> clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>>
                (new MAST::IsotropicMaterialPropertyCard::StiffnessMatrix2D(*this));
            }

            virtual ~StiffnessMatrix2D();
            
            
            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            
        protected:
            
            MAST::FieldFunction<Real>* _E;
            MAST::FieldFunction<Real>* _nu;
            bool _plane_stress;
            std::function<void(const Real, const Real, DenseMatrix<Real>&)>
            *_m, *_parm_parE, *_parm_parnu;
        };
        
        
        
        class StiffnessMatrix3D: public MAST::FieldFunction<DenseMatrix<Real> > {
        public:
            StiffnessMatrix3D(MAST::FieldFunction<Real>* E,
                              MAST::FieldFunction<Real>* nu);

            StiffnessMatrix3D(const MAST::IsotropicMaterialPropertyCard::StiffnessMatrix3D &f):
            MAST::FieldFunction<DenseMatrix<Real> >(f),
            _E(f._E->clone().release()),
            _nu(f._nu->clone().release()),
            _m(new std::function<void(const Real, const Real, DenseMatrix<Real>&)>(*f._m)),
            _parm_parE(new std::function<void(const Real, const Real, DenseMatrix<Real>&)>(*f._parm_parE)),
            _parm_parnu(new std::function<void(const Real, const Real, DenseMatrix<Real>&)>(*f._parm_parnu)) {
                _functions.insert(_E);
                _functions.insert(_nu);
            }

            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>> clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>>
                (new MAST::IsotropicMaterialPropertyCard::StiffnessMatrix3D(*this));
            }

            virtual ~StiffnessMatrix3D();
            
            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            
        protected:
            
            MAST::FieldFunction<Real>* _E;
            MAST::FieldFunction<Real>* _nu;
            std::function<void(const Real, const Real, DenseMatrix<Real>&)>
            *_m, *_parm_parE, *_parm_parnu;
        };
        
        
        
        /*!
         *   @returns the function object to calculate the requested quantity
         *   \par t, for an element of dimension \par dim
         */
        virtual std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>>
        get_property(MAST::MaterialPropertyMatrixType t,
                     const MAST::ElementPropertyCardBase& p,
                     const unsigned int dim) const;
        
    protected:
        
    };
    
}





#endif // __MAST_isotropic_material_property_card_h__
