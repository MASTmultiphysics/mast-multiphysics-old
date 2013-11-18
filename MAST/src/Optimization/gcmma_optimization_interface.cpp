//
//  gcmma_optimization_interface.h
//  MAST
//
//  Created by Manav Bhatia on 11/12/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

// MAST includes
#include "Optimization/gcmma_optimization_interface.h"


void
MAST::GCMMAOptimizationInterface::optimize() {
    
    /*    C********+*********+*********+*********+*********+*********+*********+
     C
     C    MAIN PROGRAM for using the globally convergent version of MMA
     C    to solve a problem defined by the user's subroutines
     C    INITI, FUNC1 and FUNC2.
     C
     C       Version "August 2007".
     C    !-----------------------------------------!
     C    !  The author of this program is          !
     C    !  Krister Svanberg <krille@math.kth.se>  !
     C    !-----------------------------------------!
     C
     C    The problem is assumed to be on the following form:
     C
     C      minimize  f_0(x) + z + 0.05*z^2 + sum{c_i*y_i + 0.5*(y_i)^2}
     C
     C    subject to  f_i(x) - a_i*z - y_i <= fmax_i ,   i=1,..,M
     C                       xmin_j <= x_j <= xmax_j ,   j=1,..,N
     C                                 y_i >= 0 ,        i=1,..,M
     C                                   z >= 0 .
     C
     C    where both N and M are strictly positive integers (not zero).
     C    If the user wants to solve an unconstrained problem,
     C    i.e. a problem with M=0, then  a "dummy constraint" should
     C    be introduced so that M=1. An example of such a constraint is
     C    sum{x_j} <= S, where the constant S = 1 + sum{xmax_j}.
     C
     C    The functions f_0(x) and f_i(x) should be defined by the
     C    user's subroutines
     C    FUNC1, which defines only function values, and
     C    FUNC2, which defines both function values and derivatives.
     C    The problem sizes M and N, the tolerance parameter GEPS,
     C    the constants c_i, a_i, fmax_i, xmin_j, xmax_j,
     C    and the starting values on the variables x_j,
     C    should all be given in the user's subroutine INITI.
     C    Output is specified by the users subroutine OUTXOF.
     C
     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
     C
     DIMENSION XVAL(100),XOLD1(100),XOLD2(100),XMMA(100),
     1          XMIN(100), XMAX(100), XLOW(100),XUPP(100),
     2          ALFA(100), BETA(100),DF0DX(100),
     3          A(60),B(60),C(60),Y(60),RAA(60),ULAM(60),
     4          FVAL(60),FAPP(60),FNEW(60),FMAX(60),
     5          DFDX(6000),P(6000),Q(6000),P0(100),Q0(100),
     6          UU(60),GRADF(60),DSRCH(60),HESSF(1830)
     INTEGER IYFREE(60)
     C
     C********+*********+*********+*********+*********+*********+*********+
     C  The sizes of the above arrays must be at least as follows:
     C
     C    XVAL(N),XOLD1(N),XOLD2(N),XMMA(N),
     C    XMIN(N), XMAX(N), XLOW(N),XUPP(N),
     C    ALFA(N), BETA(N),DF0DX(N),
     C    A(M),B(M),C(M),Y(M),RAA(M),ULAM(M),
     C    FVAL(M),FAPP(M),FNEW(M),FMAX(M),
     C    DFDX(M*N),P(M*N),Q(M*N),P0(N),Q0(N),
     C    UU(M),GRADF(M),DSRCH(M),HESSF(M*(M+1)/2),
     C    IYFREE(M)
     C*/
    
    int N = _feval->n_vars(), M = _feval->n_eq() + _feval->n_ineq(),
    n_rel_change_iters = _feval->n_iters_relative_change();
    
    libmesh_assert_greater(M, 0);
    libmesh_assert_greater(N, 0);
    
    std::vector<Real>  XVAL(N, 0.), XOLD1(N, 0.), XOLD2(N, 0.),
    XMMA(N, 0.), XMIN(N, 0.), XMAX(N, 0.), XLOW(N, 0.), XUPP(N, 0.),
    ALFA(N, 0.), BETA(N, 0.), DF0DX(N, 0.),
    A(M, 0.), B(M, 0.), C(M, 0.), Y(M, 0.), RAA(M, 0.), ULAM(M, 0.),
    FVAL(M, 0.), FAPP(M, 0.), FNEW(M, 0.), FMAX(M, 0.),
    DFDX(M*N, 0.), P(M*N, 0.), Q(M*N, 0.), P0(N, 0.), Q0(N, 0.),
    UU(M, 0.), GRADF(M, 0.), DSRCH(M, 0.), HESSF(M*(M+1)/2, 0.),
    f0_iters(n_rel_change_iters);
    
    std::vector<int> IYFREE(M, 0);
    std::vector<bool> eval_grads(M, false);
    
    Real F0VAL, F0NEW, F0APP, RAA0, Z, GEPS=_feval->tolerance();
    
    
    /*C********+*********+*********+*********+*********+*********+*********+
     C
     C  The meaning of some of the scalars and vectors in the program:
     C
     C     N  = Number of variables x_j in the problem.
     C     M  = Number of constraints in the problem (not including
     C          the simple upper and lower bounds on the variables).
     C INNMAX = Maximal number of inner iterations within each outer iter.
     C          A reasonable choice is INNMAX=10.
     C  ITER  = Current outer iteration number ( =1 the first iteration).
     C  GEPS  = Tolerance parameter for the constraints.
     C          (Used in the termination criteria for the subproblem.)
     C
     C   XVAL(j) = Current value of the variable x_j.
     C  XOLD1(j) = Value of the variable x_j one iteration ago.
     C  XOLD2(j) = Value of the variable x_j two iterations ago.
     C   XMMA(j) = Optimal value of x_j in the MMA subproblem.
     C   XMIN(j) = Original lower bound for the variable x_j.
     C   XMAX(j) = Original upper bound for the variable x_j.
     C   XLOW(j) = Value of the lower asymptot l_j.
     C   XUPP(j) = Value of the upper asymptot u_j.
     C   ALFA(j) = Lower bound for x_j in the MMA subproblem.
     C   BETA(j) = Upper bound for x_j in the MMA subproblem.
     C    F0VAL  = Value of the objective function f_0(x)
     C   FVAL(i) = Value of the i:th constraint function f_i(x).
     C  DF0DX(j) = Derivative of f_0(x) with respect to x_j.
     C   FMAX(i) = Right hand side of the i:th constraint.
     C   DFDX(k) = Derivative of f_i(x) with respect to x_j,
     C             where k = (j-1)*M + i.
     C      P(k) = Coefficient p_ij in the MMA subproblem, where
     C             k = (j-1)*M + i.
     C      Q(k) = Coefficient q_ij in the MMA subproblem, where
     C             k = (j-1)*M + i.
     C     P0(j) = Coefficient p_0j in the MMA subproblem.
     C     Q0(j) = Coefficient q_0j in the MMA subproblem.
     C      B(i) = Right hand side b_i in the MMA subproblem.
     C    F0APP  = Value of the approximating objective function
     C             at the optimal soultion of the MMA subproblem.
     C   FAPP(i) = Value of the approximating i:th constraint function
     C             at the optimal soultion of the MMA subproblem.
     C    RAA0   = Parameter raa_0 in the MMA subproblem.
     C    RAA(i) = Parameter raa_i in the MMA subproblem.
     C      Y(i) = Value of the "artificial" variable y_i.
     C      Z    = Value of the "minimax" variable z.
     C      A(i) = Coefficient a_i for the variable z.
     C      C(i) = Coefficient c_i for the variable y_i.
     C   ULAM(i) = Value of the dual variable lambda_i.
     C  GRADF(i) = Gradient component of the dual objective function.
     C  DSRCH(i) = Search direction component in the dual subproblem.
     C  HESSF(k) = Hessian matrix component of the dual function.
     C IYFREE(i) = 0 for dual variables which are fixed to zero in
     C               the current subspace of the dual subproblem,
     C           = 1 for dual variables which are "free" in
     C               the current subspace of the dual subproblem.
     C
     C********+*********+*********+*********+*********+*********+*********+*/
    
    
    /*C
     C  The USER should now give values to the parameters
     C  M, N, GEPS, XVAL (starting point),
     C  XMIN, XMAX, FMAX, A and C.
     C*/
    // _initi(M,N,GEPS,XVAL,XMIN,XMAX,FMAX,A,C);
    // Assumed:  FMAX == A
    _feval->init_dvar(XVAL, XMIN, XMAX);
    // set the value of C[i] to be very large numbers
    Real max_x = 0.;
    for (unsigned int i=0; i<M; i++)
        if (max_x < fabs(XVAL[i]))
            max_x = fabs(XVAL[i]);
    std::fill(C.begin(), C.end(), std::max(1.e6*max_x, 1.e6));
    //IF(M.EQ.0) GOTO 100
    //IF(N.EQ.0) GOTO 100
    
    int INNMAX=15, ITER=0, ITE=0, INNER=0, ICONSE=0;
    /*C
     C  The USER should now calculate function values at XVAL.
     C  The result should be put in F0VAL and FVAL.
     C*/
    std::fill(eval_grads.begin(), eval_grads.end(), false);
    _feval->evaluate(XVAL,
                     F0VAL, false, DF0DX,
                     FVAL, eval_grads, DFDX);
    /*C
     C  The USER may now write the current (starting) solution.
     C*/
    _feval->output(ITER, XVAL, F0VAL, FVAL);
    /*C
     C  The outer iterative process starts.
     C*/
    bool terminate = false, inner_terminate=false;
    while (!terminate) {//30   CONTINUE
        
        ITER=ITER+1;
        ITE=ITE+1;
        /*C
         C  The USER should now calculate function values and gradients
         C  at XVAL. The result should be put in F0VAL,DF0DX,FVAL,DFDX.
         C*/
        std::fill(eval_grads.begin(), eval_grads.end(), true);
        _feval->evaluate(XVAL,
                         F0VAL, true, DF0DX,
                         FVAL, eval_grads, DFDX);
        
        /*C
         C  RAA0,RAA,XLOW,XUPP,ALFA and BETA are calculated.
         C*/
        raasta_(&M, &N, &RAA0, &RAA[0], &XMIN[0], &XMAX[0], &DF0DX[0], &DFDX[0]);
        asympg_(&ITER, &M, &N, &XVAL[0], &XMIN[0], &XMAX[0], &XOLD1[0], &XOLD2[0],
                &XLOW[0], &XUPP[0], &ALFA[0], &BETA[0]);
        /*C      CALL ASYMPG(ITER,M,N,XVAL,XMIN,XMAX,XOLD1,XOLD2,
         C     1            XLOW,XUPP,ALFA,BETA,RAA0,RAA)
         C
         C  The inner iterative process starts.
         C*/
        INNER=0;
        inner_terminate = false;
        while (!inner_terminate) {
            
            //40   CONTINUE
            /*C
             C  The subproblem is generated and solved.
             C*/
            mmasug_(&ITER, &M, &N, &GEPS, &IYFREE[0], &XVAL[0], &XMMA[0],
                    &XMIN[0], &XMAX[0], &XLOW[0], &XUPP[0], &ALFA[0], &BETA[0],
                    &A[0], &B[0], &C[0], &Y[0], &Z, &RAA0, &RAA[0], &ULAM[0],
                    &F0VAL, &FVAL[0], &F0APP, &FAPP[0], &FMAX[0], &DF0DX[0], &DFDX[0],
                    &P[0], &Q[0], &P0[0], &Q0[0], &UU[0], &GRADF[0], &DSRCH[0], &HESSF[0]);
            /*C
             C  The USER should now calculate function values at XMMA.
             C  The result should be put in F0NEW and FNEW.
             C*/
            std::fill(eval_grads.begin(), eval_grads.end(), false);
            _feval->evaluate(XMMA,
                             F0NEW, false, DF0DX,
                             FNEW, eval_grads, DFDX);
            
            if (INNER >= INNMAX) inner_terminate = true;
            /*C
             C  It is checked if the approximations were conservative.
             C*/
            if (!inner_terminate) {
                conser_( &M, &ICONSE, &GEPS, &F0NEW, &F0APP, &FNEW[0], &FAPP[0]);
                if (ICONSE == 1) inner_terminate = true;
            }
            /*C
             C  The approximations were not conservative, so RAA0 and RAA
             C  are updated and one more inner iteration is started.
             C*/
            if (!inner_terminate) {
                INNER=INNER+1;
                raaupd_( &M, &N, &GEPS, &XMMA[0], &XVAL[0],
                        &XMIN[0], &XMAX[0], &XLOW[0], &XUPP[0],
                        &F0NEW, &FNEW[0], &F0APP, &FAPP[0], &RAA0, &RAA[0]);
            }
            //GOTO 40
        } //60   CONTINUE
        
        /*C
         C  The inner iterative process has terminated, which means
         C  that an outer iteration has been completed.
         C  The variables are updated so that XVAL stands for the new
         C  outer iteration point. The fuction values are also updated.
         C*/
        xupdat_( &N, &ITER, &XMMA[0], &XVAL[0], &XOLD1[0], &XOLD2[0]);
        fupdat_( &M, &F0NEW, &FNEW[0], &F0VAL, &FVAL[0]);
        /*C
         C  The USER may now write the current solution.
         C*/
        _feval->output(ITER, XVAL, F0VAL, FVAL);
        f0_iters[(ITE-1)%n_rel_change_iters] = F0VAL;
        
        /*C
         C  One more outer iteration is started as long as
         C  ITE is less than MAXITE:
         C*/
        if (ITE == _feval->max_iters()) {
            libMesh::out
            << "GCMMA: Reached maximum iterations, terminating! "
            << std::endl;
            terminate = true;
        }
        
        // relative change in objective
        bool rel_change_conv = true;
        Real f0_curr = f0_iters[n_rel_change_iters-1];
        
        for (unsigned int i=0; i<n_rel_change_iters-1; i++) {
            if (f0_curr > sqrt(GEPS))
                rel_change_conv = (rel_change_conv &&
                                   fabs(f0_iters[i]-f0_curr)/fabs(f0_curr) < GEPS);
            else
                rel_change_conv = (rel_change_conv &&
                                   fabs(f0_iters[i]-f0_curr) < GEPS);
        }
        if (rel_change_conv) {
            libMesh::out
            << "GCMMA: Converged relative change tolerance, terminating! "
            << std::endl;
            terminate = true;
        }
        
    }//100  CONTINUE
}
