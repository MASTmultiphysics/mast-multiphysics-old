//
//  Plot2DBase.h
//  FESystem
//
//  Created by Manav Bhatia on 4/14/11.
//  Copyright 2011 . All rights reserved.
//


#ifndef __fesystem_plot2d_base_h__
#define __fesystem_plot2d_base_h__

// C++ includes
#include <string>
#include <vector>


// FESystem includes
#include "Base/FESystemExceptions.h"
#include "Plotting/AxisBase.h"


namespace FESystem
{
    // Forward declerations
    namespace Numerics {template <typename ValType> class VectorBase;}
    namespace Numerics {template <typename ValType> class MatrixBase;}
    
    
    namespace Plotting
    {
        // Forward declerations
        template <typename ValType> class AxisBase;
        
        template <typename ValType>
        class Plot2DBase
        {
        public:
            
            Plot2DBase(FESystem::Plotting::AxisKind& a1,
                       FESystem::Plotting::AxisKind& a2);
            
            virtual ~Plot2DBase();
            
            FESystem::Plotting::AxisBase<ValType>& getAxis(FESystemUInt i);
            
            void setTitle(const std::string& s);
            
            virtual FESystemUInt getDimension();
            
            virtual void plotData2D(const FESystem::Numerics::VectorBase<ValType>& x,
                                    const FESystem::Numerics::VectorBase<ValType>& y) = 0;

            virtual void plotData3DSurf(const FESystem::Numerics::VectorBase<ValType>& x,
                                        const FESystem::Numerics::VectorBase<ValType>& y,
                                        const FESystem::Numerics::MatrixBase<ValType>& z) = 0;
                        
        protected:
            
            std::vector<FESystem::Plotting::AxisBase<ValType>* > axes;
            
            std::string title;
            
        };
        
        DeclareException2(IncompatibleDimensions, 
						  FESystemUInt,
						  FESystemUInt,
						  << "Dimension of X and Y data must be same. \n"
                          << "X : " << Arg1 
						  << "Y : " << Arg2 );

    }
}


#endif // __fesystem_plot2d_base_h__

