//
//  StructuralElementBase.h
//  FESystem
//
//  Created by Manav Bhatia on 8/12/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem_structural_element_base_h__
#define __fesystem_structural_element_base_h__

// C++ includes
#include <string>
#include <vector>

// FESystem includes
#include "Base/FESystemTypes.h"


namespace FESystem
{

    namespace Numerics {template <typename ValType> class MatrixBase;}
    namespace Numerics {template <typename ValType> class VectorBase;}
    namespace Quadrature {class QuadratureBase;}
    namespace FiniteElement {class FiniteElementBase;}
    namespace Mesh {class ElemBase;}
    namespace Geom {class Section;}

    
    namespace Structures
    {
        
        enum StructuralVariable
        {
            U_DISP,
            V_DISP,
            W_DISP,
            THETA_X,
            THETA_Y,
            THETA_Z
        };
        
        

        enum DiscreteFormulation
        {
            DISPLACEMENT_FORMULATION,
            MIXED_DISPLACEMENT_STRESS_FORMULATION
        };
                
        
        class StructuralElementBase
        {
        public:
            /*!
             *   Constructor
             */
            StructuralElementBase();
                        
            virtual ~StructuralElementBase();
            
            /*!
             *   Returns the number of active degrees of freedom for this element
             */
            virtual FESystemUInt getNElemDofs() const=0;
            
            /*!
             *   returns the name of the variable as a character
             */
            static std::string getVariableName(FESystem::Structures::StructuralVariable var);
            
            /*!
             *   returns the variable enumeration for the character name
             */
            static FESystem::Structures::StructuralVariable getVariableEnum(const std::string& var);
            
            /*!
             *    Returns the indices for the matrix entries that this element will contribute to. For instance, a bar element will only have the extensional stiffness
             *    matrix, and the indices will correspond to the location of the stiffness terms corresponding to the u-dofs.
             */
            virtual void getActiveElementMatrixIndices(std::vector<FESystemUInt>& vec) = 0;
            
            
            virtual void getStressTensor(const FESystem::Numerics::VectorBase<FESystemDouble>& pt, const FESystem::Numerics::VectorBase<FESystemDouble>& sol,
                                         FESystem::Numerics::MatrixBase<FESystemDouble>& mat) = 0;
            
            virtual void transformMatrixToGlobalSystem(const FESystem::Numerics::MatrixBase<FESystemDouble>& elem_mat, FESystem::Numerics::MatrixBase<FESystemDouble>& global_mat);

            virtual void transformVectorToGlobalSystem(const FESystem::Numerics::VectorBase<FESystemDouble>& elem_vec, FESystem::Numerics::VectorBase<FESystemDouble>& global_vec);

        protected:
            
            virtual void clear();
            

            virtual void initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, const FESystem::Quadrature::QuadratureBase& q_rule,
                                    FESystemDouble E, FESystemDouble nu, FESystemDouble rho);

            void calculateDeformationTransformationMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat);
            
            FESystemBoolean if_initialized;
            
            const FESystem::Mesh::ElemBase* geometric_elem;
            
            const FESystem::Quadrature::QuadratureBase* quadrature;
            
            const FESystem::FiniteElement::FiniteElementBase* finite_element;
            
            FESystemDouble E_val;
            
            FESystemDouble nu_val;
            
            FESystemDouble G_val;
            
            FESystemDouble rho_val;
        };
    }
}



#endif // __fesystem_structural_element_base_h__
