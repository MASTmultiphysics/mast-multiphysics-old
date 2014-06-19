//
//  isotropic_material_property_card.h
//  MAST
//
//  Created by Manav Bhatia on 1/30/14.
//  Copyright (c) 2014 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_isotropic_material_property_card_h__
#define __MAST_isotropic_material_property_card_h__

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
        
        
        class StiffnessMatrix1D: public MAST::FieldFunction<DenseRealMatrix > {
        public:
            
            StiffnessMatrix1D( MAST::FieldFunction<Real>* E,
                              MAST::FieldFunction<Real>* nu);

            StiffnessMatrix1D(const MAST::IsotropicMaterialPropertyCard::StiffnessMatrix1D& f):
            MAST::FieldFunction<DenseRealMatrix >(f),
            _E(f._E->clone().release()),
            _nu(f._nu->clone().release()){
                _functions.insert(_E);
                _functions.insert(_nu);
            }

            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > > clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > >
                (new MAST::IsotropicMaterialPropertyCard::StiffnessMatrix1D(*this));
            }

            virtual ~StiffnessMatrix1D();
            
            virtual void operator() (const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
            
        protected:
            
            MAST::FieldFunction<Real>* _E;
            MAST::FieldFunction<Real>* _nu;
        };
        
        
        
        class TransverseShearStiffnessMatrix: public MAST::FieldFunction<DenseRealMatrix > {
        public:
            TransverseShearStiffnessMatrix( MAST::FieldFunction<Real>* E,
                                           MAST::FieldFunction<Real>* nu,
                                           MAST::FieldFunction<Real>* kappa);

            TransverseShearStiffnessMatrix(const MAST::IsotropicMaterialPropertyCard::TransverseShearStiffnessMatrix& f):
            MAST::FieldFunction<DenseRealMatrix >(f),
            _E(f._E->clone().release()),
            _nu(f._nu->clone().release()) {
                _functions.insert(_E);
                _functions.insert(_nu);
            }

            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > > clone() const  {
                return std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > >
                (new MAST::IsotropicMaterialPropertyCard::TransverseShearStiffnessMatrix(*this));
            }

            virtual ~TransverseShearStiffnessMatrix();
            
            virtual void operator() (const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
            
        protected:
            
            MAST::FieldFunction<Real>* _E;
            MAST::FieldFunction<Real>* _nu;
            MAST::FieldFunction<Real>* _kappa;
        };
        
        
        class StiffnessMatrix2D: public MAST::FieldFunction<DenseRealMatrix > {
        public:
            StiffnessMatrix2D(MAST::FieldFunction<Real>* E,
                              MAST::FieldFunction<Real>* nu,
                              bool plane_stress);

            StiffnessMatrix2D(const MAST::IsotropicMaterialPropertyCard::StiffnessMatrix2D& f):
            MAST::FieldFunction<DenseRealMatrix >(f),
            _E(f._E->clone().release()),
            _nu(f._nu->clone().release()),
            _plane_stress(f._plane_stress) {
                _functions.insert(_E);
                _functions.insert(_nu);
            }

            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > > clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > >
                (new MAST::IsotropicMaterialPropertyCard::StiffnessMatrix2D(*this));
            }

            virtual ~StiffnessMatrix2D();
            
            
            virtual void operator() (const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
            
        protected:
            
            MAST::FieldFunction<Real>* _E;
            MAST::FieldFunction<Real>* _nu;
            bool _plane_stress;
        };
        
        
        
        class StiffnessMatrix3D: public MAST::FieldFunction<DenseRealMatrix > {
        public:
            StiffnessMatrix3D(MAST::FieldFunction<Real>* E,
                              MAST::FieldFunction<Real>* nu);

            StiffnessMatrix3D(const MAST::IsotropicMaterialPropertyCard::StiffnessMatrix3D &f):
            MAST::FieldFunction<DenseRealMatrix >(f),
            _E(f._E->clone().release()),
            _nu(f._nu->clone().release()) {
                _functions.insert(_E);
                _functions.insert(_nu);
            }

            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > > clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > >
                (new MAST::IsotropicMaterialPropertyCard::StiffnessMatrix3D(*this));
            }

            virtual ~StiffnessMatrix3D();
            
            virtual void operator() (const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
            
        protected:
            
            MAST::FieldFunction<Real>* _E;
            MAST::FieldFunction<Real>* _nu;
        };
        
        
        
        class ThermalExpansionMatrix: public MAST::FieldFunction<DenseRealMatrix > {
        public:
            
            ThermalExpansionMatrix(unsigned int dim,
                                   MAST::FieldFunction<Real>* alpha):
            MAST::FieldFunction<DenseRealMatrix >("ThermalExpansionMatrix"),
            _dim(dim),
            _alpha(alpha) {
                _functions.insert(_alpha);
            }
            
            ThermalExpansionMatrix(const MAST::IsotropicMaterialPropertyCard::ThermalExpansionMatrix& f):
            MAST::FieldFunction<DenseRealMatrix >(f),
            _dim(f._dim),
            _alpha(f._alpha->clone().release()) {
                _functions.insert(_alpha);
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > > clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > >
                (new MAST::IsotropicMaterialPropertyCard::ThermalExpansionMatrix(*this));
            }
            
            virtual ~ThermalExpansionMatrix() {
                delete _alpha;
            }
            
            virtual void operator() (const libMesh::Point& p, const Real t, DenseRealMatrix& m) const {
                
                Real alpha;
                (*_alpha)(p, t, alpha);
                switch (_dim) {
                    case 1:
                        m.resize(2,1);
                        break;
                        
                    case 2:
                        m.resize(3,1);
                        break;
                        
                    case 3:
                        m.resize(3,1);
                        break;
                }
                
                for (unsigned int i=0; i<_dim; i++)
                    m(i,0) = alpha;
            }
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const libMesh::Point& p, const Real t, DenseRealMatrix& m) const {
                Real alpha;
                _alpha->partial(f, p, t, alpha);
                switch (_dim) {
                    case 1:
                        m.resize(2,1);
                        break;
                        
                    case 2:
                        m.resize(3,1);
                        break;
                        
                    case 3:
                        m.resize(3,1);
                        break;
                }
                
                for (unsigned int i=0; i<_dim; i++)
                    m(i,0) = alpha;

            }
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const libMesh::Point& p, const Real t, DenseRealMatrix& m) const {
                Real alpha;
                _alpha->total(f, p, t, alpha);
                switch (_dim) {
                    case 1:
                        m.resize(2,1);
                        break;
                        
                    case 2:
                        m.resize(3,1);
                        break;
                        
                    case 3:
                        m.resize(3,1);
                        break;
                }
                
                for (unsigned int i=0; i<_dim; i++)
                    m(i,0) = alpha;
            }
            
            
        protected:

            const unsigned int _dim;
            MAST::FieldFunction<Real>* _alpha;
        };

        
        
        /*!
         *   @returns the function object to calculate the requested quantity
         *   \par t, for an element of dimension \par dim
         */
        virtual std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > >
        get_property(MAST::MaterialPropertyMatrixType t,
                     const MAST::ElementPropertyCardBase& p,
                     const unsigned int dim) const;
        
    protected:
        
    };
    
}





#endif // __MAST_isotropic_material_property_card_h__
