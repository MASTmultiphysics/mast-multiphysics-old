//
//  PKFlutterSolver.h
//  FESystem
//
//  Created by Manav Bhatia on 7/24/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem_pk_flutter_solver_h__
#define __fesystem_pk_flutter_solver_h__


// FESystem includes
#include "Disciplines/Aeroelasticity/FrequencyDomainFlutterSolverBase.h"

namespace FESystem
{
    namespace Aeroelasticity
    {
        template <typename ValType>
        class PKFlutterSolver: public FESystem::Aeroelasticity::FrequencyDomainFlutterSolverBase<ValType> 
        {
        public:
            PKFlutterSolver();
            
            virtual ~PKFlutterSolver();
            
            virtual void solve();
            
        protected:
            
        };
    }
}


#endif // __fesystem_pk_flutter_solver_h__
