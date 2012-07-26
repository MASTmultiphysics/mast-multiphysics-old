//
//  UGFlutterSolver.h
//  FESystem
//
//  Created by Manav Bhatia on 7/24/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem_ug_flutter_solver_h__
#define __fesystem_ug_flutter_solver_h__


// FESystem includes
#include "Disciplines/Aeroelasticity/FrequencyDomainFlutterSolverBase.h"

namespace FESystem
{
    namespace Aeroelasticity
    {
        template <typename ValType>
        class UGFlutterSolver: public FESystem::Aeroelasticity::FrequencyDomainFlutterSolverBase<ValType>
        {
        public:
            UGFlutterSolver();
            
            virtual ~UGFlutterSolver();
            
            virtual void solve();
            
        protected:
            
        };
    }
}

#endif // __fesystem_ug_flutter_solver_h__
