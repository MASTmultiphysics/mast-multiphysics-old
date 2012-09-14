//
//  PistonTheory2D.h
//  FESystem
//
//  Created by Manav Bhatia on 9/13/12.
//
//

#ifndef __FESystem__PistonTheory2D__
#define __FESystem__PistonTheory2D__

// FESystem includes
#include "Disciplines/Structure/Structural2DElementBase.h"


namespace FESystem
{
    namespace Geometry {class Point;}
    
    namespace Structures
    {
        class PistonTheory2D: public FESystem::Structures::Structural2DElementBase
        {
        public:
            PistonTheory2D();
            
            virtual ~PistonTheory2D();
            
            virtual void clear();
            
            virtual FESystemUInt getNElemDofs() const;
            
            /*!
             *    Returns the indices for the matrix entries that this element will contribute to. For instance, a bar element will only have the extensional stiffness
             *    matrix, and the indices will correspond to the location of the stiffness terms corresponding to the u-dofs.
             */
            virtual void getActiveElementMatrixIndices(std::vector<FESystemUInt>& vec);
            
            virtual void initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, const FESystem::Quadrature::QuadratureBase& q_rule,
                                    FESystemUInt order, FESystemDouble m, FESystemDouble a, FESystemDouble g);
            
            virtual void calculateForceVector(const FESystem::Numerics::VectorBase<FESystemDouble>& sol, const FESystem::Numerics::VectorBase<FESystemDouble>& vel, FESystem::Numerics::VectorBase<FESystemDouble>& force);
            
            virtual void calculateTangentMatrix(const FESystem::Numerics::VectorBase<FESystemDouble>& sol, const FESystem::Numerics::VectorBase<FESystemDouble>& vel,
                                                FESystem::Numerics::MatrixBase<FESystemDouble>& dfdw, FESystem::Numerics::MatrixBase<FESystemDouble>& dfdwdot);
            
            virtual void getStressTensor(const FESystem::Geometry::Point& pt, const FESystem::Numerics::VectorBase<FESystemDouble>& sol,
                                         FESystem::Numerics::MatrixBase<FESystemDouble>& mat) ;
            
        protected:
            
            FESystemUInt theory_order;
            
            FESystemDouble mach;
            
            FESystemDouble a_inf;
            
            FESystemDouble gamma;
            
            FESystemDouble u_inf;
        };
    }
}


#endif /* defined(__FESystem__PistonTheory2D__) */
