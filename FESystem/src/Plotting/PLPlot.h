//
//  PLPlot.h
//  FESystem
//
//  Created by Manav Bhatia on 4/14/11.
//  Copyright 2011 . All rights reserved.
//


#ifndef __fesystem_plplot_h__
#define __fesystem_plplot_h__

// C++ includes
#include <string>


// PLPlot includes
#include "plplot/plstream.h"


// FESystem includes
#include "Plotting/Plot2DBase.h"


namespace FESystem
{
    namespace Plotting
    {
        template <typename ValType>
        class PLPlot: public FESystem::Plotting::Plot2DBase<ValType>
        {
        public:
            
            PLPlot(FESystem::Plotting::AxisKind a1,
                   FESystem::Plotting::AxisKind a2);
            
            virtual ~PLPlot();
                        
            virtual void plotData2D(const FESystem::Numerics::VectorBase<ValType>& x,
                                    const FESystem::Numerics::VectorBase<ValType>& y);

            virtual void plotData3DSurf(const FESystem::Numerics::VectorBase<ValType>& x,
                                        const FESystem::Numerics::VectorBase<ValType>& y,
                                        const FESystem::Numerics::MatrixBase<ValType>& z);

            virtual void plot2DMatrix(const FESystem::Numerics::MatrixBase<ValType>& z);

        protected:
            
            plstream* pl_stream;
        };
    }
}


#endif // __fesystem_plplot_h__

