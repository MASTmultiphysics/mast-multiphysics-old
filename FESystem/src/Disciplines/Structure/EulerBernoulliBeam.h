//
//  EulerBernoulliBeam.h
//  FESystem
//
//  Created by Manav Bhatia on 8/13/12.
//
//

#ifndef FESystem_EulerBernoulliBeam_h
#define FESystem_EulerBernoulliBeam_h


// FESystem includes
#include "Disciplines/Structure/LinearBeamElementBase.h"


namespace FESystem
{
    namespace Geometry {class Point;}
    
    namespace Structures
    {
        class EulerBernoulliBeam: public FESystem::Structures::LinearBeamElementBase
        {
        public:
            EulerBernoulliBeam();
            
            virtual ~EulerBernoulliBeam();
            
            virtual void calculateConsistentMassMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat);
            
            virtual void calculateStiffnessMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat);
            
        protected:
                        
        };
    }
}


#endif
