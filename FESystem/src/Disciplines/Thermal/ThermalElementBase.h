//
//  ThermalElementBase.h
//  FESystem
//
//  Created by Manav Bhatia on 8/12/12.
//
//

#ifndef __fesystem_thermal_element_base_h__
#define __fesystem_thermal_element_base_h__


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
    namespace Geometry {class Point;}
    namespace Mesh {class ElemBase;}
    
    
    namespace Thermal
    {
        
        enum ThermalVariable
        {
            TEMP
        };
        
        
        
        class ThermalElementBase
        {
        public:
            /*!
             *   Constructor
             */
            ThermalElementBase();
            
            virtual ~ThermalElementBase();
            
            /*!
             *   Returns the number of active degrees of freedom for this element
             */
            virtual FESystemUInt getNElemDofs() const=0;
            
            /*!
             *   returns the name of the variable as a character
             */
            static std::string getVariableName(FESystem::Thermal::ThermalVariable var);
            
            /*!
             *   returns the variable enumeration for the character name
             */
            static FESystem::Thermal::ThermalVariable getVariableEnum(const std::string& var);
            
            virtual void clear();
            
            
            virtual void initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, const FESystem::Quadrature::QuadratureBase& q_rule,
                                    FESystemDouble k, FESystemDouble cp, FESystemDouble rho, FESystemDouble cs = 0.0);

            
            
            virtual void getFluxVector(const FESystem::Numerics::VectorBase<FESystemDouble>& pt, const FESystem::Numerics::VectorBase<FESystemDouble>& sol,
                                       FESystem::Numerics::VectorBase<FESystemDouble>& vec);

            
            virtual void calculateConsistentCapacitanceMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat);
            
            void calculateDiagonalCapaciatanceMatrix(FESystem::Numerics::VectorBase<FESystemDouble>& vec);

            virtual void calculateConductanceMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat);
            
        protected:
            
            void calculateCapacitanceOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat);

            void calculateConductanceOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat);
            
            void getMaterialCapacitanceMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat);
            
            void getMaterialConductanceMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat);
                        
            
            FESystemBoolean if_initialized;
            
            const FESystem::Mesh::ElemBase* geometric_elem;
            
            const FESystem::Quadrature::QuadratureBase* quadrature;
            
            const FESystem::FiniteElement::FiniteElementBase* finite_element;
            
            FESystemDouble k_val;
            
            FESystemDouble cp_val;
            
            FESystemDouble rho_val;
            
            FESystemDouble cross_section_val;
        };
    }
}



#endif // __fesystem_thermal_element_base_h__

