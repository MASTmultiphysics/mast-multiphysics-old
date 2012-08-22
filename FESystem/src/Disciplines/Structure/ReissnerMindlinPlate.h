//
//  ReissnerMindlinPlate.h
//  FESystem
//
//  Created by Manav Bhatia on 8/13/12.
//
//

#ifndef FESystem_ReissnerMindlinPlate_h
#define FESystem_ReissnerMindlinPlate_h

// FESystem includes
#include "Disciplines/Structure/Structural2DElementBase.h"


namespace FESystem
{
    namespace Geometry {class Point;}
    
    namespace Structures
    {
        class ReissnerMindlinPlate: public FESystem::Structures::Structural2DElementBase
        {
        public:
            ReissnerMindlinPlate();
            
            virtual ~ReissnerMindlinPlate();
            
            virtual void clear();
            
            virtual void initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, const FESystem::Quadrature::QuadratureBase& q_bend,
                                    const FESystem::Quadrature::QuadratureBase& q_shear, FESystemDouble E, FESystemDouble nu, FESystemDouble rho, FESystemDouble th);
            
            virtual void calculateConsistentMassMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat);

            virtual void calculateDiagonalMassMatrix(FESystem::Numerics::VectorBase<FESystemDouble>& vec);

            virtual void calculateStiffnessMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat);
            
            virtual void transformMatrixToGlobalCoordinate(const std::vector<FESystem::Structures::StructuralVariable>& vars,
                                                           const FESystem::Numerics::MatrixBase<FESystemDouble>& elem_cs_mat,
                                                           FESystem::Numerics::MatrixBase<FESystemDouble>& global_cs_mat);
            
            
            virtual void getStressTensor(const FESystem::Numerics::VectorBase<FESystemDouble>& pt, const FESystem::Numerics::VectorBase<FESystemDouble>& sol,
                                         FESystem::Numerics::MatrixBase<FESystemDouble>& mat) ;

            
        protected:
            
            void getMaterialMassMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat);
            
            void getMaterialComplianceMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& bend_mat, FESystem::Numerics::MatrixBase<FESystemDouble>& shear_mat);
            
            void calculateInertiaOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat);
            
            void calculateBendingOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat);
            
            void calculateShearOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat);
            
            const FESystem::Quadrature::QuadratureBase*  quadrature_bending;
            
            const FESystem::Quadrature::QuadratureBase*  quadrature_shear;
            
            FESystemDouble kappa;
        };
    }
}


#endif
