//
//  ElemBoundaryBase.h
//  FESystem
//
//  Created by Manav Bhatia on 2/14/13.
//
//

#ifndef __FESystem__ElemBoundaryBase__
#define __FESystem__ElemBoundaryBase__

// C++ includes
#include <memory>

// FESystem includes
#include "Base/FESystemTypes.h"


namespace FESystem
{
    namespace Mesh
    {
        // Forward declerations
        class ElemBase;
        
        class ElemBoundaryBase
        {
        public:
            ElemBoundaryBase();
            
            virtual ~ElemBoundaryBase();
            
            const FESystem::Mesh::ElemBase& getParentElement() const;

            FESystem::Mesh::ElemBase& getParentElement();
            
            FESystemUInt getBoundaryNumberOfParentElement() const;
            
            FESystemUInt getDimension() const;

            std::auto_ptr<FESystem::Mesh::ElemBase*> createElementForBoundary() const;
            
        protected:
            
            FESystem::Mesh::ElemBase* parent_elem;
            
            FESystemUInt boundary_id;
        };
    }
}


#endif /* defined(__FESystem__ElemBoundaryBase__) */
