//
//  PLPlot.cpp
//  FESystem
//
//  Created by Manav Bhatia on 4/14/11.
//  Copyright 2011 . All rights reserved.
//

// C++ includes
#include <iostream>
#include <iomanip>

#ifdef HAVE_PLPLOT

// FESystem includes
#include "Plotting/PLPlot.h"
#include "Plotting/AxisBase.h"
#include "Numerics/VectorBase.h"
#include "Numerics/MatrixBase.h"


template <typename ValType>
FESystem::Plotting::PLPlot<ValType>::PLPlot(FESystem::Plotting::AxisKind a1,
                                                FESystem::Plotting::AxisKind a2):
FESystem::Plotting::Plot2DBase<ValType>(a1, a2)
{
    this->pl_stream = new plstream();
    this->pl_stream->scolbg(255,255,255); // set background to white
    this->pl_stream->scol0(1,0,0,0); // set axis colors to black
    this->pl_stream->start("aqt",1,1);
    this->pl_stream->env(-20, 20, -20, 20, 0, 2);
    //this->pl_stream->fontld( 1 );
    //this->pl_stream->init();
}


template <typename ValType>
FESystem::Plotting::PLPlot<ValType>::~PLPlot()
{
    if (this->pl_stream != NULL)
        delete this->pl_stream;
    this->pl_stream = NULL;
}
        

template <typename ValType>
void
FESystem::Plotting::PLPlot<ValType>::plotData2D(const FESystem::Numerics::VectorBase<ValType>& x,
                                                const FESystem::Numerics::VectorBase<ValType>& y)
{
    FESystemAssert2(x.getSize() == y.getSize(), 
                    FESystem::Plotting::IncompatibleDimensions,
                    x.getSize(), y.getSize());
    
    
    std::pair<ValType, ValType> r1 = this->getAxis(0).getRange(),
    r2 = this->getAxis(0).getRange();
    ValType xmin=x.getMinVal(), xmax=x.getMaxVal(), ymin=y.getMinVal(), ymax=y.getMaxVal();
    if (xmin == xmax)
        xmax = xmax+0.5;
    if (ymin == ymax)
        ymax = ymax+0.5;
    
    this->pl_stream->col0(1);
//    this->pl_stream->env(xmin, xmax, ymin, ymax, 0, 2);
//    this->pl_stream->adv(0);
//    this->pl_stream->vsta();
//    this->pl_stream->col0(2);
//    this->pl_stream->lab(this->getAxis(0).getLabel().c_str(),
//                         this->getAxis(1).getLabel().c_str(),
//                         this->title.c_str());
    
    //this->pl_stream->prec(1,1);
    //this->pl_stream->wind(0,20,0,20);
    this->pl_stream->col0(4);

    this->pl_stream->line(x.getSize(), 
                          const_cast<ValType*>(&(x.getVectorValues())[0]),
                          const_cast<ValType*>(&(y.getVectorValues())[0]));
}



template <typename ValType>
void
FESystem::Plotting::PLPlot<ValType>::plotData3DSurf(const FESystem::Numerics::VectorBase<ValType>& x,
                                                    const FESystem::Numerics::VectorBase<ValType>& y,
                                                    const FESystem::Numerics::MatrixBase<ValType>& z)
{
    std::pair<FESystemUInt, FESystemUInt> s=z.getSize();

    FESystemAssert2(x.getSize() == s.first, FESystem::Plotting::IncompatibleDimensions, x.getSize(), s.first);
    FESystemAssert2(y.getSize() == s.second, FESystem::Plotting::IncompatibleDimensions, y.getSize(), s.second);
    
    std::pair<ValType, ValType> r1 = this->getAxis(0).getRange(),
    r2 = this->getAxis(0).getRange();
    ValType xmin=x.getMinVal(), xmax=x.getMaxVal(), ymin=y.getMinVal(), ymax=y.getMaxVal(), zmin=z.getMinVal(), zmax=z.getMaxVal();
    if (xmin == xmax) xmax += 0.5;
    if (ymin == ymax) ymax += 0.5;
    if (zmin == zmax) zmax += 0.5;
    
    this->pl_stream->col0(1);
//    this->pl_stream->env(xmin, xmax, ymin, ymax, 0, 2);
//    this->pl_stream->adv(0);
//    this->pl_stream->vsta();
//    this->pl_stream->col0(2);
//    this->pl_stream->lab(this->getAxis(0).getLabel().c_str(),
//                         this->getAxis(1).getLabel().c_str(),
//                         this->title.c_str());
    
    //this->pl_stream->prec(1,1);
    this->pl_stream->adv( 0 );
    this->pl_stream->col0( 1 );
    this->pl_stream->vpor( 0.0, 1.0, 0.0, 0.9 );
    this->pl_stream->wind( -1.0, 1.0, -1.0, 1.5 );
    this->pl_stream->w3d(1.0, 1.0, 1.0, xmin, xmax, ymin, ymax, zmin, zmax, 30*xmax, 20.0);
    this->pl_stream->box3( "bnstu", "x axis", 0.0, 0,
                          "bnstu", "y axis", 0.0, 0,
                          "bcdmnstuv", "z axis", 0.0, 4 );
    this->pl_stream->col0(4);
    
    ValType** z_vals;
    FESystemUInt n_cont=10;
    std::vector<ValType> clevels(n_cont);
    for (FESystemUInt i=0; i<n_cont; i++)
        clevels[i] = zmin+(zmax-zmin)/(n_cont+1)*i;
    this->pl_stream->Alloc2dGrid(&z_vals, s.first, s.second);
    
    for (FESystemUInt i=0; i<s.first; i++)
        for (FESystemUInt j=0; j<s.second; j++)
            z_vals[i][j] = z.getVal(i,j);
    
    this->pl_stream->meshc(const_cast<ValType*>(&(x.getVectorValues())[0]),
                           const_cast<ValType*>(&(y.getVectorValues())[0]),
                           z_vals, x.getSize(), y.getSize(), DRAW_LINEXY|MAG_COLOR|BASE_CONT, &(clevels[0]), n_cont);

    
    this->pl_stream->Free2dGrid(z_vals, s.first, s.second);
}




template <typename ValType>
void
FESystem::Plotting::PLPlot<ValType>::plot2DMatrix(const FESystem::Numerics::MatrixBase<ValType>& z)
{
    std::pair<FESystemUInt, FESystemUInt> s=z.getSize();
    
    std::pair<ValType, ValType> r1 = this->getAxis(0).getRange(),
    r2 = this->getAxis(0).getRange();
    ValType zmin=z.getMinVal(), zmax=z.getMaxVal(), xmin=0.0, xmax=1.0, ymin=0.0, ymax=1.0;
    if (zmin == zmax) zmax += 0.5;
    
    this->pl_stream->col0(1);
    this->pl_stream->env(xmin, xmax, ymin, ymax, 0, 2);
    this->pl_stream->col0(4);
    
    ValType** z_vals;
    this->pl_stream->Alloc2dGrid(&z_vals, s.first, s.second);
    
    for (FESystemUInt i=0; i<s.first; i++)
        for (FESystemUInt j=0; j<s.second; j++)
            z_vals[i][j] = z.getVal(i,j);
    
    this->pl_stream->image(z_vals, s.first, s.second, xmin, xmax, ymin, ymax, zmin, zmax, xmin, xmax, ymin, ymax);
    
    
    this->pl_stream->Free2dGrid(z_vals, s.first, s.second);
}



/***************************************************************************************/
// Template instantiations for some generic classes

template class FESystem::Plotting::PLPlot<FESystemDouble>;


/***************************************************************************************/

#endif // HAVE_PLPLOT
