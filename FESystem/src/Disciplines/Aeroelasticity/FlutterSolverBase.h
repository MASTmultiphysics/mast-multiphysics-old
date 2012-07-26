//
//  FlutterSolverBase.h
//  FESystem
//
//  Created by Manav Bhatia on 7/24/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem_flutter_solver_base_h__
#define __fesystem_flutter_solver_base_h__


namespace FESystem
{
    namespace Aeroelasticity
    {
        /*!
         *   Base class for flutter solver
         */
        template <typename ValType>
        class FlutterSolverBase
        {
        public:
            /*!
             *   Constructor
             */
            FlutterSolverBase();
            
            virtual ~FlutterSolverBase();
            
            virtual void solve()=0;
            
        protected:
            
        };
    }
}



#endif // __fesystem_flutter_solver_base_h__
