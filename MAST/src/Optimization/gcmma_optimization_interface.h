//
//  gcmma_optimization_interface.h
//  MAST
//
//  Created by Manav Bhatia on 11/12/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_gcmma_optimization_interface_h__
#define __MAST_gcmma_optimization_interface_h__

// MAST includes
#include "Optimization/optimization_interface.h"

extern "C" {
    extern void raasta_(int *M, int *N,
                        double *RAA0, double *RAA,
                        double *XMIN, double *XMAX,
                        double *DF0DX, double *DFDX);
    extern void asympg_(int *ITER, int *M, int *N,
                        double *XVAL, double *XMIN, double *XMAX,
                        double *XOLD1, double *XOLD2,
                        double *XLOW, double *XUPP,
                        double *ALFA, double *BETA);
    extern void mmasug_(int *ITER, int *M, int *N, double *GEPS, int *IYFREE,
                        double *XVAL, double *XMMA,
                        double *XMIN, double *XMAX,
                        double *XLOW, double *XUPP,
                        double *ALFA, double *BETA,
                        double *A, double *B, double *C, double *Y, double *Z,
                        double *RAA0, double *RAA, double *ULAM,
                        double *F0VAL, double *FVAL,
                        double *F0APP, double *FAPP,
                        double *FMAX, double *DF0DX, double *DFDX,
                        double *P, double *Q, double *P0, double *Q0,
                        double *UU, double *GRADF, double *DSRCH, double *HESSF);
    extern void conser_(int *M, int *ICONSE,
                        double *GEPS, double *F0NEW, double *F0APP,
                        double *FNEW, double *FAPP);
    extern void raaupd_(int *M, int *N, double *GEPS,
                        double *XMMA, double *XVAL,
                        double *XMIN, double *XMAX,
                        double *XLOW, double *XUPP,
                        double *F0NEW, double *FNEW,
                        double *F0APP, double *FAPP,
                        double *RAA0, double *RAA);
    extern void xupdat_(int *N, int *ITER,
                        double *XMMA, double *XVAL,
                        double *XOLD1, double *XOLD2);
    extern void fupdat_(int *M, double *F0NEW, double *FNEW,
                        double *F0VAL, double *FVAL);
}


namespace MAST {
    
    class GCMMAOptimizationInterface: public MAST::OptimizationInterface {
        
    public:
        
        GCMMAOptimizationInterface()
        { }
        
        virtual ~GCMMAOptimizationInterface()
        { }
        
        virtual void optimize();
        
        
    protected:
        
    };
}





#endif // __MAST_gcmma_optimization_interface_h__
