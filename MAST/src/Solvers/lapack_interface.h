//
//  lapack_interface.h
//  MAST
//
//  Created by Manav Bhatia on 8/8/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//


#ifndef MAST_lapack_interface_h
#define MAST_lapack_interface_h


// MAST includes
#include "Base/MAST_data_types.h"


extern "C" {
extern int zggev_(char* jobvl,
                  char* jobvr,
                  int* n,
                  std::complex<double>* a,
                  int* lda,
                  std::complex<double>* b,
                  int* ldb,
                  std::complex<double>* alpha,
                  std::complex<double>* beta,
                  std::complex<double>* vl,
                  int* ldvl,
                  std::complex<double>* vr,
                  int* ldvr,
                  std::complex<double>* work,
                  int* lwork,
                  double* rwork,
                  int* info);
/*
 *  Arguments
 *  =========
 *
 *  JOBVL   (input) CHARACTER*1
 *          = 'N':  do not compute the left generalized eigenvectors;
 *          = 'V':  compute the left generalized eigenvectors.
 *
 *  JOBVR   (input) CHARACTER*1
 *          = 'N':  do not compute the right generalized eigenvectors;
 *          = 'V':  compute the right generalized eigenvectors.
 *
 *  N       (input) INTEGER
 *          The order of the matrices A, B, VL, and VR.  N >= 0.
 *
 *  A       (input/output) COMPLEX array, dimension (LDA, N)
 *          On entry, the matrix A in the pair (A,B).
 *          On exit, A has been overwritten.
 *
 *  LDA     (input) INTEGER
 *          The leading dimension of A.  LDA >= max(1,N).
 *
 *  B       (input/output) COMPLEX array, dimension (LDB, N)
 *          On entry, the matrix B in the pair (A,B).
 *          On exit, B has been overwritten.
 *
 *  LDB     (input) INTEGER
 *          The leading dimension of B.  LDB >= max(1,N).
 *
 *  ALPHA   (output) COMPLEX array, dimension (N)
 *  BETA    (output) COMPLEX array, dimension (N)
 *          On exit, ALPHA(j)/BETA(j), j=1,...,N, will be the
 *          generalized eigenvalues.
 *
 *          Note: the quotients ALPHA(j)/BETA(j) may easily over- or
 *          underflow, and BETA(j) may even be zero.  Thus, the user
 *          should avoid naively computing the ratio alpha/beta.
 *          However, ALPHA will be always less than and usually
 *          comparable with norm(A) in magnitude, and BETA always less
 *          than and usually comparable with norm(B).
 *
 *  VL      (output) COMPLEX array, dimension (LDVL,N)
 *          If JOBVL = 'V', the left generalized eigenvectors u(j) are
 *          stored one after another in the columns of VL, in the same
 *          order as their eigenvalues.
 *          Each eigenvector is scaled so the largest component has
 *          abs(real part) + abs(imag. part) = 1.
 *          Not referenced if JOBVL = 'N'.
 *
 *  LDVL    (input) INTEGER
 *          The leading dimension of the matrix VL. LDVL >= 1, and
 *          if JOBVL = 'V', LDVL >= N.
 *
 *  VR      (output) COMPLEX array, dimension (LDVR,N)
 *          If JOBVR = 'V', the right generalized eigenvectors v(j) are
 *          stored one after another in the columns of VR, in the same
 *          order as their eigenvalues.
 *          Each eigenvector is scaled so the largest component has
 *          abs(real part) + abs(imag. part) = 1.
 *          Not referenced if JOBVR = 'N'.
 *
 *  LDVR    (input) INTEGER
 *          The leading dimension of the matrix VR. LDVR >= 1, and
 *          if JOBVR = 'V', LDVR >= N.
 *
 *  WORK    (workspace/output) COMPLEX array, dimension (MAX(1,LWORK))
 *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
 *
 *  LWORK   (input) INTEGER
 *          The dimension of the array WORK.  LWORK >= max(1,2*N).
 *          For good performance, LWORK must generally be larger.
 *
 *          If LWORK = -1, then a workspace query is assumed; the routine
 *          only calculates the optimal size of the WORK array, returns
 *          this value as the first entry of the WORK array, and no error
 *          message related to LWORK is issued by XERBLA.
 *
 *  RWORK   (workspace/output) REAL array, dimension (8*N)
 *
 *  INFO    (output) INTEGER
 *          = 0:  successful exit
 *          < 0:  if INFO = -i, the i-th argument had an illegal value.
 *          =1,...,N:
 *                The QZ iteration failed.  No eigenvectors have been
 *                calculated, but ALPHA(j) and BETA(j) should be
 *                correct for j=INFO+1,...,N.
 *          > N:  =N+1: other then QZ iteration failed in SHGEQZ,
 *                =N+2: error return from STGEVC.
 */
}

class LAPACK_ZGGEV{
    
public:
    
    LAPACK_ZGGEV():
    info_val(-1)
    { }
    
    /*!
     *    computes the eigensolution for A x = \lambda B x. A & B will be 
     *    overwritten
     */
    void compute(ComplexMatrixX& A, ComplexMatrixX& B,
                 bool computeEigenvectors = true);
  
    ComputationInfo info() const;
    
    const ComplexVectorX& alphas() const {
        libmesh_assert(info_val == 0);
        return this->alpha;
    }

    const ComplexVectorX& betas() const {
        libmesh_assert(info_val == 0);
        return this->beta;
    }
    
    const ComplexMatrixX& left_eigenvectors() const {
        libmesh_assert(info_val == 0);
        return this->VL;
    }

    const ComplexMatrixX& right_eigenvectors() const {
        libmesh_assert(info_val == 0);
        return this->VR;
    }

protected:
    
    ComplexMatrixX VL;
    
    ComplexMatrixX VR;
    
    ComplexVectorX alpha;
    
    ComplexVectorX beta;
    
    int info_val;
};


void
LAPACK_ZGGEV::compute(ComplexMatrixX &A, ComplexMatrixX &B,
                      bool computeEigenvectors)
{
    libmesh_assert(A.cols() == A.rows() &&
                   B.cols() == A.rows() &&
                   B.cols() == B.rows());

    int n = (int)A.cols();

    char L='N',R='N';
    
    if (computeEigenvectors)
    {
        L = 'V'; R = 'V';
        VL.setZero(n, n);
        VR.setZero(n, n);
    }

    int lwork=16*n, l_rwork=8*n; info_val=-1;

    alpha.setZero(n); beta.setZero(n);
    ComplexVectorX work; work.setZero(lwork);
    RealVectorX rwork; rwork.setZero(l_rwork);
    
    Complex *a_vals = A.data(),  *b_vals = B.data(),
    *alpha_v = alpha.data(), *beta_v = beta.data(),
    *VL_v = VL.data(), *VR_v = VR.data(),
    *work_v = work.data();
    
    Real *rwork_v = rwork.data();
    
    
    zggev_(&L, &R, &n,
           &(a_vals[0]), &n,
           &(b_vals[0]), &n,
           &(alpha_v[0]), &(beta_v[0]),
           &(VL_v[0]), &n, &(VR_v[0]), &n,
           &(work_v[0]), &lwork,
           &(rwork_v[0]),
           &info_val);
    
}



#endif
