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
    
    class IsotropicMaterialPropertyCard: public MAST::MaterialPropertyCardBase {
        
    public:
        
        IsotropicMaterialPropertyCard(unsigned int pid):
        MAST::MaterialPropertyCardBase (pid)
        { }
        
        
        class StiffnessMatrix1D: public MAST::FieldFunction<DenseMatrix<Real> > {
        public:
            
            StiffnessMatrix1D( MAST::FieldFunction<Real>& E,
                              MAST::FieldFunction<Real>& nu );
            
            virtual ~StiffnessMatrix1D() { }
            
            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void partial_derivative (const MAST::SensitivityParameters& par,
                                             const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void total_derivative (const MAST::SensitivityParameters& par,
                                           const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            
        protected:
            
            MAST::FieldFunction<Real>& _E;
            MAST::FieldFunction<Real>& _nu;
        };
        
        
        
        class TransverseShearStiffnessMatrix: public MAST::FieldFunction<DenseMatrix<Real> > {
        public:
            TransverseShearStiffnessMatrix( MAST::FieldFunction<Real>& E,
                                           MAST::FieldFunction<Real>& nu,
                                           MAST::FieldFunction<Real>& kappa);
            
            virtual ~TransverseShearStiffnessMatrix() { }
            
            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void partial_derivative (const MAST::SensitivityParameters& par,
                                             const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            
            virtual void total_derivative (const MAST::SensitivityParameters& par,
                                           const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            
        protected:
            
            MAST::FieldFunction<Real>& _E;
            MAST::FieldFunction<Real>& _nu;
            MAST::FieldFunction<Real>& _kappa;
        };
        
        
        
        
        class InertiaMatrix: public MAST::FieldFunction<DenseMatrix<Real> > {
        public:
            InertiaMatrix( MAST::FieldFunction<Real>& rho);
            
            virtual ~InertiaMatrix() { }
            
            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void partial_derivative (const MAST::SensitivityParameters& par,
                                             const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            
            virtual void total_derivative (const MAST::SensitivityParameters& par,
                                           const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            
        protected:
            
            MAST::FieldFunction<Real>& _rho;
        };
        
        
        
        
        class StiffnessMatrix2D: public MAST::FieldFunction<DenseMatrix<Real> > {
        public:
            StiffnessMatrix2D(MAST::FieldFunction<Real>& E,
                              MAST::FieldFunction<Real>& nu,
                              bool plane_stress = true);
            
            virtual ~StiffnessMatrix2D() { }
            
            
            void set_plane_stress(bool val);
            
            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void partial_derivative (const MAST::SensitivityParameters& par,
                                             const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            
            virtual void total_derivative (const MAST::SensitivityParameters& par,
                                           const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            
        protected:
            
            MAST::FieldFunction<Real>& _E;
            MAST::FieldFunction<Real>& _nu;
            bool _plane_stress;
        };
        
        
        
        class StiffnessMatrix3D: public MAST::FieldFunction<DenseMatrix<Real> > {
        public:
            StiffnessMatrix3D(MAST::FieldFunction<Real>& E,
                              MAST::FieldFunction<Real>& nu);
            
            virtual ~StiffnessMatrix3D() { }
            
            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            virtual void partial_derivative (const MAST::SensitivityParameters& par,
                                             const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            
            virtual void total_derivative (const MAST::SensitivityParameters& par,
                                           const Point& p, const Real t, DenseMatrix<Real>& m) const;
            
            
        protected:
            
            MAST::FieldFunction<Real>& _E;
            MAST::FieldFunction<Real>& _nu;
        };
        
        
        
        /*!
         *   @returns the function object to calculate the requested quantity
         *   \par t, for an element of dimension \par dim
         */
        virtual
        std::auto_ptr<MAST::FunctionBase>
        get_property(MAST::MaterialPropertyMatrixType t, const unsigned int dim);
        
    protected:
        
    };
    
}





#endif // __MAST_isotropic_material_property_card_h__
