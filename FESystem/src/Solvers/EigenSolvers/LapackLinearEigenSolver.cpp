//
//  LapackLinearEigenSolver.cpp
//  FESystem
//
//  Created by Manav Bhatia on 5/3/12.
//  Copyright (c) 2012. All rights reserved.
//

// Lapack includes
//#include "clapack.h"

// FESystem includes
#include "Solvers/EigenSolvers/LapackLinearEigenSolver.h"
#include "Solvers/Factorizations/HouseHolderTriangulation.h"
#include "Solvers/Factorizations/TriangularBacksubstitution.h"
#include "Numerics/VectorBase.h"
#include "Numerics/MatrixBase.h"
#include "Base/macros.h"

// declaration of LAPACK routines
extern "C"
{
    // Generalized Hermitian
    extern int dggev_(char*, char*, int*, double*, int*, double*, int*, double*, double*,
                      double*, double*, int*, double*, int*, double*, int*, int*);
    // Hermitian
    extern int dgeev_(char*, char*, int*, double*, int*, double*, double*, double*,
                      int*, double*, int*, double*, int*, int*);
    // Generalized Hermitian
    extern int sggev_(char*, char*, int*, float*, int*, float*, int*, float*, float*,
                      float*, float*, int*, float*, int*, float*, int*, int*);
    // Hermitian
    extern int sgeev_(char*, char*, int*, float*, int*, float*, float*, float*,
                      int*, float*, int*, float*, int*, int*);
    extern int cggev_(char*, char*, int*, std::complex<float>*, int*, std::complex<float>*, int*,
                      std::complex<float>*, std::complex<float>*, std::complex<float>*, int*,
                      std::complex<float>*, int*, std::complex<float>*, int*, float*, int*);
}

template <typename ValType>
FESystem::EigenSolvers::LapackLinearEigenSolver<ValType>::LapackLinearEigenSolver():
FESystem::EigenSolvers::LinearEigenSolverBase<ValType>()
{
    
}


template <typename ValType>
FESystem::EigenSolvers::LapackLinearEigenSolver<ValType>::~LapackLinearEigenSolver()
{
    
}



template <>
void
FESystem::EigenSolvers::LapackLinearEigenSolver<FESystemDouble>::solve()
{
    FESystemAssert0(this->matrices_are_set, FESystem::EigenSolvers::MatrixNotSet);
    
    // LAPACK solver is only for local dense matrices
    std::auto_ptr<FESystem::Numerics::MatrixBase<FESystemDouble> >
    tmatA(FESystem::Numerics::MatrixCreate<FESystemDouble>(FESystem::Numerics::LOCAL_DENSE_MATRIX)),
    tmatB(FESystem::Numerics::MatrixCreate<FESystemDouble>(FESystem::Numerics::LOCAL_DENSE_MATRIX));
    
    std::pair<FESystemUInt, FESystemUInt> s = this->A_mat->getSize();
    
    tmatA->resize(s.first, s.second);
    tmatA->copyMatrix(*(this->A_mat));
    FESystemDouble* a_vals=tmatA->getMatrixValues();
    
    switch (this->getEigenProblemType()) {
        case FESystem::EigenSolvers::HERMITIAN:
        {
            FESystemInt lwork=16*s.second, info=0;
            std::vector<FESystemDouble> alpha_r(s.first), alpha_i(s.first), work(lwork);
            FESystemDouble* eig_vec = this->eig_vec_mat_right->getMatrixValues();
            char N='N',V='V';
            std::pair<FESystemInt, FESystemInt> s_int = s;
            
            dgeev_(&N, &V, &(s_int.second), &(a_vals[0]), &(s_int.first),
                   &(alpha_r[0]), &(alpha_i[0]), &(eig_vec[0]), &(s_int.first), &(eig_vec[0]), &(s_int.first), &(work[0]), &lwork, &info);
            
            // check the convergence
            if (info == 0)
            {
                std::auto_ptr<FESystem::Numerics::VectorBase<FESystemDouble> >
                tmp_vec1(FESystem::Numerics::VectorCreate<FESystemDouble>(FESystem::Numerics::LOCAL_VECTOR).release());
                tmp_vec1->resize(s.first);
                
                this->n_converged_eig_vals = s.first;
                for (FESystemUInt i=0; i<s.first; i++)
                {
                    // rescale the eigenvectors so that the x^T B x = 1
                    this->eig_vec_mat_right->getColumnVals(i, 0, s.first-1, *tmp_vec1);
                    tmp_vec1->scale(1.0/sqrt(tmp_vec1->dotProduct(*tmp_vec1)));
                    this->eig_vec_mat_right->setColumnVals(i, 0, s.first-1, *tmp_vec1);
                    
                    this->eig_val_vec->setVal(i, alpha_r[i]);
                    //FESystemAssert0(alpha_i[i]==0.0, FESystem::Exception::InvalidValue);
                }
            }
            else if (info < 0)
            {
                std::cout << "Error in argument no.: " << info << std::endl;
                FESystemAssert0(false, FESystem::Exception::InvalidValue);
            }
            else
            {
                std::cout << "Eigensolver did not converge" << std::endl;
                FESystemAssert0(false, FESystem::Exception::InvalidValue);
            }
            
            this->n_converged_eig_vals = s.first;
        }
            break;
            
            
            
        case FESystem::EigenSolvers::NONHERMITIAN:
        {
            FESystemInt lwork=16*s.second, info=0;
            std::vector<FESystemDouble> alpha_r(s.first), alpha_i(s.first), work(lwork), eig_vec_right(s.first*s.first), eig_vec_left(s.first*s.first);
            char N='V',V='V';
            std::pair<FESystemInt, FESystemInt> s_int = s;
            
            dgeev_(&N, &V, &(s_int.second), &(a_vals[0]), &(s_int.first),
                   &(alpha_r[0]), &(alpha_i[0]), &(eig_vec_left[0]), &(s_int.first), &(eig_vec_right[0]), &(s_int.first), &(work[0]), &lwork, &info);
            
            // check the convergence
            if (info == 0)
            {
                std::auto_ptr<FESystem::Numerics::VectorBase<FESystemComplexDouble> >
                tmp_vec1(FESystem::Numerics::VectorCreate<FESystemComplexDouble>(FESystem::Numerics::LOCAL_VECTOR).release()),
                tmp_vec2(FESystem::Numerics::VectorCreate<FESystemComplexDouble>(FESystem::Numerics::LOCAL_VECTOR).release());;
                tmp_vec1->resize(s.first); tmp_vec2->resize(s.first);
                
                // first copy the eigenvectors and eigenvalues to the solver data structure
                // any complex eigenvalue appears in conjugate pairs, and the associated eigenvector would
                // occupy two consecutive columns in the eigenvector matrices. Hence, one should look for conjugate pairs
                
                this->n_converged_eig_vals = s.first;
                FESystemBoolean if_conjugate = false;
                FESystemUInt eig_id = 0;
                while (eig_id < s.first)
                {
                    // if this is the last eigenvalue, then it is not conjugate
                    if (eig_id == s.first-1)
                        if_conjugate = false;
                    else
                    {
                        // if the real part for two consecutive eigenvalues is the same, and imaginary part
                        // is opposite in sign, then it is a complex conjugate pair
                        if ((alpha_r[eig_id] == alpha_r[eig_id+1]) && (alpha_i[eig_id] == -alpha_i[eig_id+1]))
                            if_conjugate = true;
                        else
                            if_conjugate = false;
                    }
                    
                    // look at the two
                    if (if_conjugate)
                    {
                        this->eig_val_vec_complex->setVal(  eig_id,   std::complex<FESystemDouble>(alpha_r[eig_id],   alpha_i[eig_id]));
                        this->eig_val_vec_complex->setVal(eig_id+1, std::complex<FESystemDouble>(alpha_r[eig_id+1], alpha_i[eig_id+1]));
                        // first eigenvector (right)
                        for (FESystemUInt i=0; i<s.first; i++)
                            tmp_vec1->setVal(i, std::complex<FESystemDouble>(eig_vec_right[eig_id*s.first+i], eig_vec_right[(eig_id+1)*s.first+i]));
                        // do the same for the left eigenvector
                        for (FESystemUInt i=0; i<s.first; i++)
                            tmp_vec2->setVal(i, std::complex<FESystemDouble>(eig_vec_left[eig_id*s.first+i], eig_vec_left[(eig_id+1)*s.first+i]));                        
                        this->eig_vec_mat_right_complex->setColumnVals(eig_id, 0, s.first-1, *tmp_vec1);
                        this->eig_vec_mat_left_complex->setColumnVals(eig_id, 0, s.first-1, *tmp_vec2);

                        // next place the complex conjugate of this eigenvector in the next column
                        tmp_vec1->complexConjugate();
                        tmp_vec2->complexConjugate();
                        this->eig_vec_mat_right_complex->setColumnVals(eig_id+1, 0, s.first-1, *tmp_vec1);
                        this->eig_vec_mat_left_complex->setColumnVals(eig_id+1, 0, s.first-1, *tmp_vec2);
                        
                        eig_id+=2;
                    }
                    else
                    {
                        this->eig_val_vec_complex->setVal(  eig_id,   std::complex<FESystemDouble>(alpha_r[eig_id],   alpha_i[eig_id]));
                        // first eigenvector (right)
                        for (FESystemUInt i=0; i<s.first; i++)
                            tmp_vec1->setVal(i, std::complex<FESystemDouble>(eig_vec_right[eig_id*s.first+i], 0.0));
                        this->eig_vec_mat_right_complex->setColumnVals(eig_id, 0, s.first-1, *tmp_vec1);
                        
                        // do the same for the left eigenvector
                        for (FESystemUInt i=0; i<s.first; i++)
                            tmp_vec2->setVal(i, std::complex<FESystemDouble>(eig_vec_left[eig_id*s.first+i], 0.0));
                        this->eig_vec_mat_left_complex->setColumnVals(eig_id, 0, s.first-1, *tmp_vec2);
                        
                        eig_id+=1;
                    }
                }
            }
            else if (info < 0)
            {
                std::cout << "Error in argument no.: " << info << std::endl;
                FESystemAssert0(false, FESystem::Exception::InvalidValue);
            }
            else
            {
                std::cout << "Eigensolver did not converge" << std::endl;
                FESystemAssert0(false, FESystem::Exception::InvalidValue);
            }
            
            this->n_converged_eig_vals = s.first;
        }
            break;
            

            
        case FESystem::EigenSolvers::GENERALIZED_HERMITIAN:
        {
            tmatB->resize(s.first, s.second);
            tmatB->copyMatrix(*(this->B_mat));
            FESystemDouble* b_vals=tmatB->getMatrixValues();
            
            FESystemInt lwork=16*s.second, info=0;
            std::vector<FESystemDouble> alpha_r(s.first), alpha_i(s.second), beta(s.second), work(lwork);
            FESystemDouble* eig_vec = this->eig_vec_mat_right->getMatrixValues();
            char N='N',V='V';
            std::pair<FESystemInt, FESystemInt> s_int = s;
            
            dggev_(&N, &V, &(s_int.second), &(a_vals[0]), &(s_int.first),
                   &(b_vals[0]), &(s_int.first), &(alpha_r[0]), &(alpha_i[0]), &(beta[0]),
                   &(eig_vec[0]), &(s_int.first), &(eig_vec[0]), &(s_int.first), &(work[0]), &lwork, &info);
            
            // check the convergence
            if (info == 0)
            {
                std::auto_ptr<FESystem::Numerics::VectorBase<FESystemDouble> >
                tmp_vec1(FESystem::Numerics::VectorCreate<FESystemDouble>(FESystem::Numerics::LOCAL_VECTOR).release()),
                tmp_vec2(FESystem::Numerics::VectorCreate<FESystemDouble>(FESystem::Numerics::LOCAL_VECTOR).release());
                tmp_vec1->resize(s.first);
                tmp_vec2->resize(s.first);
                
                this->n_converged_eig_vals = s.first;
                for (FESystemUInt i=0; i<s.first; i++)
                {
                    // rescale the eigenvectors so that the x^T B x = 1
                    this->eig_vec_mat_right->getColumnVals(i, 0, s.first-1, *tmp_vec1);
                    this->getBMatrix().rightVectorMultiply(*tmp_vec1, *tmp_vec2);
                    tmp_vec1->scale(1.0/sqrt(tmp_vec1->dotProduct(*tmp_vec2)));
                    this->eig_vec_mat_right->setColumnVals(i, 0, s.first-1, *tmp_vec1);
                    
                    this->eig_val_vec->setVal(i, alpha_r[i]/beta[i]);
                    //FESystemAssert0(alpha_i[i]==0.0, FESystem::Exception::InvalidValue);
                }
            }
            else if (info < 0)
            {
                std::cout << "Error in argument no.: " << info << std::endl;
                FESystemAssert0(false, FESystem::Exception::InvalidValue);
            }
            else
            {
                std::cout << "Eigensolver did not converge" << std::endl;
                FESystemAssert0(false, FESystem::Exception::InvalidValue);
            }
        }
            break;
            
            
            
        case FESystem::EigenSolvers::GENERALIZED_NONHERMITIAN:
        {
            tmatB->resize(s.first, s.second);
            tmatB->copyMatrix(*(this->B_mat));
            FESystemDouble* b_vals=tmatB->getMatrixValues();
            
            FESystemInt lwork=16*s.second, info=0;
            std::vector<FESystemDouble> alpha_r(s.first), alpha_i(s.second), beta(s.second), work(lwork);
            FESystemDouble* eig_vec = this->eig_vec_mat_right->getMatrixValues();
            char N='N',V='V';
            std::pair<FESystemInt, FESystemInt> s_int = s;
            
            dggev_(&N, &V, &(s_int.second), &(a_vals[0]), &(s_int.first),
                   &(b_vals[0]), &(s_int.first), &(alpha_r[0]), &(alpha_i[0]), &(beta[0]),
                   &(eig_vec[0]), &(s_int.first), &(eig_vec[0]), &(s_int.first), &(work[0]), &lwork, &info);
            
            // check the convergence
            if (info == 0)
            {
                std::auto_ptr<FESystem::Numerics::VectorBase<FESystemDouble> >
                tmp_vec1(FESystem::Numerics::VectorCreate<FESystemDouble>(FESystem::Numerics::LOCAL_VECTOR).release()),
                tmp_vec2(FESystem::Numerics::VectorCreate<FESystemDouble>(FESystem::Numerics::LOCAL_VECTOR).release());
                tmp_vec1->resize(s.first);
                tmp_vec2->resize(s.first);
                
                this->n_converged_eig_vals = s.first;
                for (FESystemUInt i=0; i<s.first; i++)
                {
                    // rescale the eigenvectors so that the x^T B x = 1
                    this->eig_vec_mat_right->getColumnVals(i, 0, s.first-1, *tmp_vec1);
                    this->getBMatrix().rightVectorMultiply(*tmp_vec1, *tmp_vec2);
                    tmp_vec1->scale(1.0/sqrt(tmp_vec1->dotProduct(*tmp_vec2)));
                    this->eig_vec_mat_right->setColumnVals(i, 0, s.first-1, *tmp_vec1);
                    
                    this->eig_val_vec->setVal(i, alpha_r[i]/beta[i]);
                    //FESystemAssert0(alpha_i[i]==0.0, FESystem::Exception::InvalidValue);
                }
            }
            else if (info < 0)
            {
                std::cout << "Error in argument no.: " << info << std::endl;
                FESystemAssert0(false, FESystem::Exception::InvalidValue);
            }
            else
            {
                std::cout << "Eigensolver did not converge" << std::endl;
                FESystemAssert0(false, FESystem::Exception::InvalidValue);
            }
        }
            break;
            
            
            
        default:
            FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
            break;
    }
    
}




template <>
void
FESystem::EigenSolvers::LapackLinearEigenSolver<FESystemFloat>::solve()
{
    FESystemAssert0(this->matrices_are_set, FESystem::EigenSolvers::MatrixNotSet);
    
    std::auto_ptr<FESystem::Numerics::MatrixBase<FESystemFloat> >
    tmatA(FESystem::Numerics::MatrixCreate<FESystemFloat>(this->A_mat->getType())),
    tmatB(FESystem::Numerics::MatrixCreate<FESystemFloat>(this->A_mat->getType()));
    
    std::pair<FESystemUInt, FESystemUInt> s = this->A_mat->getSize();
    
    tmatA->resize(s.first, s.second);
    tmatA->copyMatrix(*(this->A_mat));
    FESystemFloat* a_vals=tmatA->getMatrixValues();
    
    switch (this->getEigenProblemType()) {
        case FESystem::EigenSolvers::HERMITIAN:
        {
            FESystemInt lwork=16*s.second, info=0;
            std::vector<FESystemFloat> alpha_r(s.first), alpha_i(s.second), beta(s.second), work(lwork);
            FESystemFloat* eig_vec = this->eig_vec_mat_right->getMatrixValues();
            char N='N',V='V';
            std::pair<FESystemInt, FESystemInt> s_int = s;
            
            sgeev_(&N, &V, &(s_int.second), &(a_vals[0]), &(s_int.first),
                   &(alpha_r[0]), &(alpha_i[0]), &(eig_vec[0]), &(s_int.first), &(eig_vec[0]), &(s_int.first), &(work[0]), &lwork, &info);
            
            // check the convergence
            if (info == 0)
            {
                std::auto_ptr<FESystem::Numerics::VectorBase<FESystemFloat> >
                tmp_vec1(FESystem::Numerics::VectorCreate<FESystemFloat>(FESystem::Numerics::LOCAL_VECTOR).release());
                tmp_vec1->resize(s.first);
                
                this->n_converged_eig_vals = s.first;
                for (FESystemUInt i=0; i<s.first; i++)
                {
                    // rescale the eigenvectors so that the x^T B x = 1
                    this->eig_vec_mat_right->getColumnVals(i, 0, s.first-1, *tmp_vec1);
                    tmp_vec1->scale(1.0/sqrt(tmp_vec1->dotProduct(*tmp_vec1)));
                    this->eig_vec_mat_right->setColumnVals(i, 0, s.first-1, *tmp_vec1);
                    
                    this->eig_val_vec->setVal(i, alpha_r[i]);
                    //FESystemAssert0(alpha_i[i]==0.0, FESystem::Exception::InvalidValue);
                }
            }
            else if (info < 0)
            {
                std::cout << "Error in argument no.: " << info << std::endl;
                FESystemAssert0(false, FESystem::Exception::InvalidValue);
            }
            else
            {
                std::cout << "Eigensolver did not converge" << std::endl;
                FESystemAssert0(false, FESystem::Exception::InvalidValue);
            }
            
            this->n_converged_eig_vals = s.first;
        }
            break;
            
        case FESystem::EigenSolvers::GENERALIZED_HERMITIAN:
        {
            tmatB->resize(s.first, s.second);
            tmatB->copyMatrix(*(this->B_mat));
            FESystemFloat* b_vals=tmatB->getMatrixValues();
            
            FESystemInt lwork=16*s.second, info=0;
            std::vector<FESystemFloat> alpha_r(s.first), alpha_i(s.second), beta(s.second), work(lwork);
            FESystemFloat* eig_vec = this->eig_vec_mat_right->getMatrixValues();
            char N='N',V='V';
            std::pair<FESystemInt, FESystemInt> s_int = s;
            
            sggev_(&N, &V, &(s_int.second), &(a_vals[0]), &(s_int.first),
                   &(b_vals[0]), &(s_int.first), &(alpha_r[0]), &(alpha_i[0]), &(beta[0]),
                   &(eig_vec[0]), &(s_int.first), &(eig_vec[0]), &(s_int.first), &(work[0]), &lwork, &info);
            
            // check the convergence
            if (info == 0)
            {
                std::auto_ptr<FESystem::Numerics::VectorBase<FESystemFloat> >
                tmp_vec1(FESystem::Numerics::VectorCreate<FESystemFloat>(FESystem::Numerics::LOCAL_VECTOR).release()),
                tmp_vec2(FESystem::Numerics::VectorCreate<FESystemFloat>(FESystem::Numerics::LOCAL_VECTOR).release());
                tmp_vec1->resize(s.first);
                tmp_vec2->resize(s.first);
                
                this->n_converged_eig_vals = s.first;
                for (FESystemUInt i=0; i<s.first; i++)
                {
                    // rescale the eigenvectors so that the x^T B x = 1
                    this->eig_vec_mat_right->getColumnVals(i, 0, s.first-1, *tmp_vec1);
                    this->getBMatrix().rightVectorMultiply(*tmp_vec1, *tmp_vec2);
                    tmp_vec1->scale(1.0/sqrt(tmp_vec1->dotProduct(*tmp_vec2)));
                    this->eig_vec_mat_right->setColumnVals(i, 0, s.first-1, *tmp_vec1);
                    
                    this->eig_val_vec->setVal(i, alpha_r[i]/beta[i]);
                    //FESystemAssert0(alpha_i[i]==0.0, FESystem::Exception::InvalidValue);
                }
                
            }
            else if (info < 0)
            {
                std::cout << "Error in argument no.: " << info << std::endl;
                FESystemAssert0(false, FESystem::Exception::InvalidValue);
            }
            else
            {
                std::cout << "Eigensolver did not converge" << std::endl;
                FESystemAssert0(false, FESystem::Exception::InvalidValue);
            }
        }
            break;
            
        case FESystem::EigenSolvers::NONHERMITIAN:
        case FESystem::EigenSolvers::GENERALIZED_NONHERMITIAN:
        default:
            FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
            break;
    }
    
}


template <>
void
FESystem::EigenSolvers::LapackLinearEigenSolver<FESystemComplexDouble>::solve()
{
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}



template <>
void
FESystem::EigenSolvers::LapackLinearEigenSolver<FESystemComplexFloat>::solve()
{
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}


/***************************************************************************************/
// Template instantiations for some generic classes

INSTANTIATE_CLASS_FOR_ALL_DATA_TYPES(FESystem::EigenSolvers::LapackLinearEigenSolver);

/***************************************************************************************/
