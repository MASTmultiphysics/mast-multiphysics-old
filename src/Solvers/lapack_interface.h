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
     *  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
     *          On entry, the matrix A in the pair (A,B).
     *          On exit, A has been overwritten.
     *
     *  LDA     (input) INTEGER
     *          The leading dimension of A.  LDA >= max(1,N).
     *
     *  B       (input/output) DOUBLE PRECISION array, dimension (LDB, N)
     *          On entry, the matrix B in the pair (A,B).
     *          On exit, B has been overwritten.
     *
     *  LDB     (input) INTEGER
     *          The leading dimension of B.  LDB >= max(1,N).
     *
     *  ALPHAR  (output) DOUBLE PRECISION array, dimension (N)
     *  ALPHAI  (output) DOUBLE PRECISION array, dimension (N)
     *  BETA    (output) DOUBLE PRECISION array, dimension (N)
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
     *  VR      (output) DOUBLE PRECISION array, dimension (LDVR,N)
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
     *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
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
    extern int dggev_(char* jobvl,
                      char* jobvr,
                      int* n,
                      double* a,
                      int* lda,
                      double* b,
                      int* ldb,
                      double* alpha_r,
                      double* alpha_i,
                      double* beta,
                      double* vl,
                      int* ldvl,
                      double* vr,
                      int* ldvr,
                      double* work,
                      int* lwork,
                      int* info);
    
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
    
    const ComplexMatrixX& A() const {
        libmesh_assert(info_val == 0);
        return this->_A;
    }
    
    
    const ComplexMatrixX& B() const {
        libmesh_assert(info_val == 0);
        return this->_B;
    }
    
    
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
    
    /*!
     *    Scales the right eigenvector so that the inner product with respect
     *    to the B matrix is equal to an Identity matrix, i.e.
     *    VL* B * VR = I
     */
    void scale_eigenvectors_to_identity_innerproduct() {
        libmesh_assert(info_val == 0);
        
        // this product should be an identity matrix
        ComplexMatrixX r = this->VL.conjugate().transpose() * _B * this->VR;
        
        // scale the right eigenvectors by the inverse of the inner-product
        // diagonal
        Complex val;
        for (unsigned int i=0; i<_B.cols(); i++) {
            val = r(i,i);
            if (std::abs(val) > 0.)
                this->VR.col(i) *= (1./val);
        }
    }
    
    void print_inner_product(std::ostream& out) const {
        libmesh_assert(info_val == 0);
        ComplexMatrixX r;
        r = this->VL.conjugate().transpose() * _A * this->VR;
        out << "conj(VL)' * A * VR" << std::endl
        << r << std::endl;
        
        r = this->VL.conjugate().transpose() * _B * this->VR;
        out << "conj(VL)' * B * VR" << std::endl
        << r << std::endl;
        
    }
    
protected:
    
    ComplexMatrixX _A;
    
    ComplexMatrixX _B;
    
    ComplexMatrixX VL;
    
    ComplexMatrixX VR;
    
    ComplexVectorX alpha;
    
    ComplexVectorX beta;
    
    int info_val;
};



class LAPACK_DGGEV{
    
public:
    
    LAPACK_DGGEV():
    info_val(-1)
    { }
    
    /*!
     *    computes the eigensolution for A x = \lambda B x. A & B will be
     *    overwritten
     */
    void compute(RealMatrixX& A, RealMatrixX& B,
                 bool computeEigenvectors = true);
    
    ComputationInfo info() const;
    
    const RealMatrixX& A() const {
        libmesh_assert(info_val == 0);
        return this->_A;
    }
    
    
    const RealMatrixX& B() const {
        libmesh_assert(info_val == 0);
        return this->_B;
    }
    
    
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
    
    /*!
     *    Scales the right eigenvector so that the inner product with respect
     *    to the B matrix is equal to an Identity matrix, i.e.
     *    VL* B * VR = I
     */
    void scale_eigenvectors_to_identity_innerproduct() {
        libmesh_assert(info_val == 0);
        
        // this product should be an identity matrix
        ComplexMatrixX r = this->VL.conjugate().transpose() * _B * this->VR;
        
        // scale the right eigenvectors by the inverse of the inner-product
        // diagonal
        Complex val;
        for (unsigned int i=0; i<_B.cols(); i++) {
            val = r(i,i);
            if (std::abs(val) > 0.)
                this->VR.col(i) *= (1./val);
        }
    }
    
    void print_inner_product(std::ostream& out) const {
        libmesh_assert(info_val == 0);
        ComplexMatrixX r;
        r = this->VL.conjugate().transpose() * _A * this->VR;
        out << "conj(VL)' * A * VR" << std::endl
        << r << std::endl;
        
        r = this->VL.conjugate().transpose() * _B * this->VR;
        out << "conj(VL)' * B * VR" << std::endl
        << r << std::endl;
        
    }
    
protected:
    
    RealMatrixX _A;
    
    RealMatrixX _B;
    
    ComplexMatrixX VL;
    
    ComplexMatrixX VR;
    
    ComplexVectorX alpha;
    
    ComplexVectorX beta;
    
    int info_val;
};


inline void
LAPACK_ZGGEV::compute(ComplexMatrixX &A, ComplexMatrixX &B,
                      bool computeEigenvectors)
{
    libmesh_assert(A.cols() == A.rows() &&
                   B.cols() == A.rows() &&
                   B.cols() == B.rows());
    
    _A = A;
    _B = B;
    
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




inline void
LAPACK_DGGEV::compute(RealMatrixX &A, RealMatrixX &B,
                      bool computeEigenvectors)
{
    libmesh_assert(A.cols() == A.rows() &&
                   B.cols() == A.rows() &&
                   B.cols() == B.rows());
    
    _A = A;
    _B = B;
    
    int n = (int)A.cols();
    
    char L='N',R='N';
    
    if (computeEigenvectors)
    {
        L = 'V'; R = 'V';
        VL.setZero(n, n);
        VR.setZero(n, n);
    }
    
    int lwork=16*n; info_val=-1;
    
    alpha.setZero(n); beta.setZero(n);
    
    RealVectorX work, aval_r, aval_i, bval;
    RealMatrixX vecl, vecr;
    work.setZero(lwork);
    aval_r.setZero(n); aval_i.setZero(n); bval.setZero(n);
    vecl.setZero(n,n); vecr.setZero(n,n);
    
    Real *a_vals = A.data(),  *b_vals = B.data(),
    *alpha_r_v = aval_r.data(), *alpha_i_v = aval_i.data(),
    *beta_v = bval.data(), *vecl_v = vecl.data(), *vecr_v = vecr.data(),
    *work_v = work.data();
    
    
    dggev_(&L, &R, &n,
           &(a_vals[0]), &n,
           &(b_vals[0]), &n,
           &(alpha_r_v[0]), &(alpha_i_v[0]), &(beta_v[0]),
           &(vecl_v[0]), &n, &(vecr_v[0]), &n,
           &(work_v[0]), &lwork,
           &info_val);
    
    // now sort the eigenvalues for complex conjugates
    unsigned int n_located = 0;
    while (n_located < n) {
        // if the imaginary part of the eigenvalue is non-zero, it is a
        // complex conjugate
        if (aval_i(n_located) != 0.) { // complex conjugate
            alpha(  n_located) = std::complex<double>(aval_r(n_located),  aval_i(n_located));
            alpha(1+n_located) = std::complex<double>(aval_r(n_located), -aval_i(n_located));
            beta (  n_located) = bval(n_located);
            beta (1+n_located) = bval(n_located);
            // copy the eigenvectors if they were requested
            if (computeEigenvectors) {
                VL.col(  n_located) = (vecl.col(  n_located).cast<Complex>() +
                                       vecl.col(1+n_located).cast<Complex>() * std::complex<double>(0, 1.));
                VL.col(1+n_located) = (vecl.col(  n_located).cast<Complex>() -
                                       vecl.col(1+n_located).cast<Complex>() * std::complex<double>(0, 1.));
                VR.col(  n_located) = (vecr.col(  n_located).cast<Complex>() +
                                       vecr.col(1+n_located).cast<Complex>() * std::complex<double>(0, 1.));
                VR.col(1+n_located) = (vecr.col(  n_located).cast<Complex>() -
                                       vecr.col(1+n_located).cast<Complex>() * std::complex<double>(0, 1.));
                n_located +=2;
            }
            else {
                VL.col(n_located) = vecl.col(n_located).cast<Complex>();
                VR.col(n_located) = vecr.col(n_located).cast<Complex>();
                n_located++;
            }
        }
            
    }
}



#endif
