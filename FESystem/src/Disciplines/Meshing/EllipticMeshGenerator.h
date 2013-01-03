//
//  EllipticMeshGenerator.h
//  FESystem
//
//  Created by Manav Bhatia on 12/23/12.
//
//

#ifndef __FESystem__EllipticMeshGenerator__
#define __FESystem__EllipticMeshGenerator__

// FESystem includes
#include "Base/FESystemTypes.h"

namespace FESystem
{
    // Forward declerations
    namespace Numerics {template<typename ValType> class VectorBase;}
    namespace Numerics {template<typename ValType> class MatrixBase;}
    namespace Quadrature {class QuadratureBase;}
    namespace FiniteElement {class FiniteElementBase;}
    namespace Mesh {class ElemBase;}
    namespace Geometry {class Point;}

    namespace Meshing
    {
        class EllipticMeshGenerator
        {
        public:
            
            EllipticMeshGenerator();
            
            virtual ~EllipticMeshGenerator();
            
            virtual void clear();
            
            void initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, const FESystem::Quadrature::QuadratureBase& q_rule, const FESystem::Numerics::VectorBase<FESystemDouble>& sol,
                            const FESystemDouble p_val, const FESystemDouble q_val);
            
            void calculateResidual(FESystem::Numerics::VectorBase<FESystemDouble>& res);
            
            void calculateTangentMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat);
            
        protected:

            FESystemBoolean if_initialized;
            
            const FESystem::Mesh::ElemBase* geometric_elem;
            
            const FESystem::Quadrature::QuadratureBase* quadrature;
            
            const FESystem::FiniteElement::FiniteElementBase* finite_element;
            
            const FESystem::Numerics::VectorBase<FESystemDouble>* solution;

            FESystemDouble P, Q;
        };
    }
}



#endif /* defined(__FESystem__EllipticMeshGenerator__) */
