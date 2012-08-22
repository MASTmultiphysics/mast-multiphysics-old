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
             *   returns the name of the variable as a character
             */
            static std::string getVariableName(FESystem::Structures::StructuralVariable var);
            
            /*!
             *   returns the variable enumeration for the character name
             */
            static FESystem::Structures::StructuralVariable getVariableEnum(const std::string& var);
            
            
            virtual void transformMatrixToGlobalCoordinate(const std::vector<FESystem::Structures::StructuralVariable>& vars,
                                                           const FESystem::Numerics::MatrixBase<FESystemDouble>& elem_cs_mat,
                                                           FESystem::Numerics::MatrixBase<FESystemDouble>& global_cs_mat) = 0;
            
            
            virtual void getStressTensor(const FESystem::Numerics::VectorBase<FESystemDouble>& pt, const FESystem::Numerics::VectorBase<FESystemDouble>& sol,
                                         FESystem::Numerics::MatrixBase<FESystemDouble>& mat) = 0;
            
        protected:
            
            virtual void clear();
            

            virtual void initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, const FESystem::Quadrature::QuadratureBase& q_rule,
                                    FESystemDouble E, FESystemDouble nu, FESystemDouble rho);

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
