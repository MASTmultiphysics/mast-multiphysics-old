//
//  FEInterpolationType.h
//  FESystem
//
//  Created by Manav Bhatia on 2/14/13.
//
//

#ifndef __FESystem__FEInterpolationType__
#define __FESystem__FEInterpolationType__

// C++ includes
#include <vector>
#include <string>


// FESystem includes
#include "FiniteElems/FiniteElementType.h"
#include "Quadrature/QuadratureType.h"
#include "Base/FESystemTypes.h"



namespace FESystem
{
    // Forward includes
    namespace Geometry {class Point;}
    
    namespace FiniteElement
    {
        class FEInterpolationType
        {
        public:
            /*!
             *   default constructor
             */
            FEInterpolationType();

            /*!
             *   Copy constructor
             */
            FEInterpolationType(const FESystem::FiniteElement::FEInterpolationType& fe);

            ~FEInterpolationType();

            void clear();
            
            void set(FESystem::FiniteElement::FiniteElementType fe, FESystem::FiniteElement::DofDistributionType dof, FESystem::FiniteElement::ContinuityType cont);
            
            void setNVariables(const FESystemUInt n);
            
            void setVariable(const FESystemUInt i, const FESystemUInt order);
            
            void setQRuleTypeForPointSpecification(FESystem::Quadrature::QuadratureType);
            
            const std::vector<FESystem::Geometry::Point*>& getDoFDistributionPoints() const;
            
            FESystemUInt getOrderForVariable(const FESystemUInt i) const;
            
            void incrementOrderForVariable(const FESystemUInt i);
            
            void decrementOrderForVariable(const FESystemUInt i);
            
            FESystem::FiniteElement::FiniteElementType getFEType() const;
            
            FESystem::FiniteElement::DofDistributionType getDofDistributionType() const;
            
            FESystem::FiniteElement::ContinuityType getContinuityType() const;
            
        protected:
            
            FESystemBoolean if_initialized;
            
            FESystem::FiniteElement::FiniteElementType fe_type;
            
            FESystem::FiniteElement::DofDistributionType dof_distribution_type;
            
            FESystem::FiniteElement::ContinuityType continuity_type;
            
            std::vector<FESystemUInt> interpolation_order_for_variable;
            
            std::vector<FESystem::Geometry::Point*> dof_points;
        };
    }
}


#endif /* defined(__FESystem__FEInterpolationType__) */
