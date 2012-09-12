//
//  VonKarmanStrain1D.h
//  FESystem
//
//  Created by Manav Bhatia on 8/13/12.
//
//

#ifndef FESystem_VonKarmanStrain1D_h
#define FESystem_VonKarmanStrain1D_h

// FESystem includes
#include "Disciplines/Structure/Structural1DElementBase.h"


namespace FESystem
{
    // Forward declerations
    namespace Geometry {class Point;}
    
    namespace Structures
    {
        // Forward declerations
        class ExtensionBar;
        class LinearBeamElementBase;
        
        class VonKarmanStrain1D: public FESystem::Structures::Structural1DElementBase
        {
        public:
            VonKarmanStrain1D();
            
            virtual ~VonKarmanStrain1D();
            
            virtual void clear();
            
            virtual FESystemUInt getNElemDofs() const;
            
            /*!
             *    Returns the indices for the matrix entries that this element will contribute to. For instance, a bar element will only have the extensional stiffness
             *    matrix, and the indices will correspond to the location of the stiffness terms corresponding to the u-dofs.
             */
            virtual void getActiveElementMatrixIndices(std::vector<FESystemUInt>& vec);

            /*!
             *    Method to initialize the vector when the bar element is used to calculate the extensional strain. In this case, the element has stiffness components 
             *    associated with the longitudinal degree of freedom as well.
             */
            virtual void initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, const FESystem::Quadrature::QuadratureBase& q_rule,
                                    FESystem::Structures::ExtensionBar& bar, FESystem::Structures::LinearBeamElementBase& beam);

            /*!
             *    Method to initialize the vector when the element has a given constant extensional stress. In this case, the element does not have stiffness components
             *    associated with the longitudinal degree of freedom as well.
             */
            virtual void initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, const FESystem::Quadrature::QuadratureBase& q_rule,
                                    FESystemDouble stress, FESystem::Structures::LinearBeamElementBase& beam);

            virtual void calculateInternalForceVector(const FESystem::Numerics::VectorBase<FESystemDouble>& sol, FESystem::Numerics::VectorBase<FESystemDouble>& vec);
            
            virtual void calculateTangentStiffnessMatrix(const FESystem::Numerics::VectorBase<FESystemDouble>& sol, FESystem::Numerics::MatrixBase<FESystemDouble>& mat);
            
            virtual void getStressTensor(const FESystem::Geometry::Point& pt, const FESystem::Numerics::VectorBase<FESystemDouble>& sol,
                                         FESystem::Numerics::MatrixBase<FESystemDouble>& mat);
 
        protected:
            
            void calculateTransverseDisplacementOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat);

            FESystemBoolean if_constant_extension_stress;
            
            FESystemDouble extension_stress;
                        
            FESystem::Structures::ExtensionBar* bar_elem;

            FESystem::Structures::LinearBeamElementBase* beam_elem;
        
        };
    }
}


#endif
