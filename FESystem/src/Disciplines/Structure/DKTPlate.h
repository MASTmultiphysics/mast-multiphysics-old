//
//  DKTPlate.h
//  FESystem
//
//  Created by Manav Bhatia on 8/13/12.
//
//

#ifndef FESystem_DKTPlate_h
#define FESystem_DKTPlate_h

// FESystem includes
#include "Disciplines/Structure/LinearPlateElementBase.h"


namespace FESystem
{
    namespace Geometry {class Point;}
    namespace Mesh {class Tri6;}
    namespace Mesh {class Node;}
    
    namespace Structures
    {
        class DKTPlate: public FESystem::Structures::LinearPlateElementBase
        {
        public:
            DKTPlate();
            
            virtual ~DKTPlate();
            
            virtual void clear();
            

            virtual void initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe_tri3, FESystem::FiniteElement::FiniteElementBase& fe_tri6,
                                    const FESystem::Quadrature::QuadratureBase& q_tri3, const FESystem::Quadrature::QuadratureBase& q_tri6,
                                    FESystemDouble E, FESystemDouble nu, FESystemDouble rho, FESystemDouble th);

            virtual void calculateConsistentMassMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat);
            
            virtual void calculateDiagonalMassMatrix(FESystem::Numerics::VectorBase<FESystemDouble>& vec);
            
            virtual void calculateStiffnessMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat);            
            
            virtual void getStressTensor(const FESystem::Numerics::VectorBase<FESystemDouble>& pt, const FESystem::Numerics::VectorBase<FESystemDouble>& sol,
                                         FESystem::Numerics::MatrixBase<FESystemDouble>& mat) ;

            
        protected:
            
            void getMaterialMassMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat);
            
            void getMaterialComplianceMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat);
            
            void calculateInertiaOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat);
            
            void calculateBendingOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat);
            
            FESystemBoolean if_first_call;
            
            FESystem::FiniteElement::FiniteElementBase* finite_element_tri6;
            
            const FESystem::Quadrature::QuadratureBase* quadrature_tri6;
            
            FESystem::Geometry::Point* origin;
            
            FESystem::Mesh::Tri6* tri6_elem;
            
            FESystem::Mesh::Node* node0;
            
            FESystem::Mesh::Node* node1;
            
            FESystem::Mesh::Node* node2;
            
            FESystem::Mesh::Node* node3;
            
            FESystem::Mesh::Node* node4;
            
            FESystem::Mesh::Node* node5;
        };
    }
}


#endif
