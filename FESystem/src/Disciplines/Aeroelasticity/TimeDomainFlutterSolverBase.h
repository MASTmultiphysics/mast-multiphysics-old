//
//  TimeDomainFlutterSolverBase.h
//  FESystem
//
//  Created by Manav Bhatia on 7/24/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem_time_domain_flutter_solver_base_h__
#define __fesystem_time_domain_flutter_solver_base_h__

namespace FESystem
{
    namespace Aeroelasticity
    {
        class FlutterSolverBase
        {
        public:
            
            FlutterSolverBase();
            
            virtual ~FlutterSolverBase();
            
            virtual void solve()=0;
            
        protected:
            
        };
    }
}


#endif
