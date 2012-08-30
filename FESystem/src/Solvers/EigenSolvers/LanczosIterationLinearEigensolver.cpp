
//
//  LanczosIterationLinearEigensolver.cpp
//  FESystem
//
//  Created by Manav Bhatia on 7/14/11.
//  Copyright 2011 . All rights reserved.
//

// C++ includes
#include <cmath>
#include <stdlib.h>

// FESystem includes 
#include "Solvers/EigenSolvers/LanczosIterationLinearEigensolver.h"
#include "Solvers/EigenSolvers/QRMethodLinearEigensolver.h"
#include "Solvers/Factorizations/ModifiedQRFactorization.h"
#include "Numerics/VectorBase.h"
#include "Numerics/MatrixBase.h"
#include "Base/macros.h"


template <typename ValType> 
FESystem::EigenSolvers::LanczosIterationLinearEigenSolver<ValType>::LanczosIterationLinearEigenSolver():
FESystem::EigenSolvers::LinearEigenSolverBase<ValType>()
{
    
}


template <typename ValType> 
FESystem::EigenSolvers::LanczosIterationLinearEigenSolver<ValType>::~LanczosIterationLinearEigenSolver()
{
    
}



//template <typename ValType> 
//void
//FESystem::EigenSolvers::LanczosIterationLinearEigenSolver<ValType>::setShift(ValType v)
//{
//    this->solver_shift = v;
//}


template <typename ValType>
void
FESystem::EigenSolvers::LanczosIterationLinearEigenSolver<ValType>::initializeMatrices()
{
    
    FESystemAssert0(this->A_mat != NULL, FESystem::Exception::NULLQuantity);
    
    this->krylov_basis_mat.reset(FESystem::Numerics::MatrixCreate<ValType>(this->getAMatrix().getType()).release());
    this->hessenberg_mat.reset(FESystem::Numerics::MatrixCreate<ValType>(this->getAMatrix().getType()).release());
    
    std::pair<FESystemUInt, FESystemUInt> s = this->getAMatrix().getSize();
    
    this->krylov_basis_mat->resize(s.first, s.second);
    this->hessenberg_mat->resize(s.first, s.second);
}



template <typename ValType> 
void
FESystem::EigenSolvers::LanczosIterationLinearEigenSolver<ValType>::solve()
{
    FESystemAssert0(this->matrices_are_set, FESystem::EigenSolvers::MatrixNotSet);
    
    this->initializeMatrices();
        
    std::auto_ptr<FESystem::Numerics::MatrixBase<ValType> > 
    tri_diag_mat(FESystem::Numerics::MatrixCreate<ValType>(FESystem::Numerics::LOCAL_DENSE_MATRIX)),
    orthogonalization_factors(FESystem::Numerics::MatrixCreate<ValType>(FESystem::Numerics::LOCAL_DENSE_MATRIX)); // used to store the orthogonalizaiton factor approximation

    std::auto_ptr<FESystem::Numerics::VectorBase<ValType> > 
    q_vec_old(FESystem::Numerics::VectorCreate<ValType>(FESystem::Numerics::LOCAL_VECTOR)),
    q_vec(FESystem::Numerics::VectorCreate<ValType>(FESystem::Numerics::LOCAL_VECTOR)),
    r_vec(FESystem::Numerics::VectorCreate<ValType>(FESystem::Numerics::LOCAL_VECTOR)),
    eig_val_convergence(FESystem::Numerics::VectorCreate<ValType>(FESystem::Numerics::LOCAL_VECTOR));
    
    std::pair<FESystemUInt, FESystemUInt> s = this->A_mat->getSize();

    orthogonalization_factors->resize(s.first, s.second); orthogonalization_factors->zero();
    
    q_vec_old->resize(s.first);
    q_vec->resize(s.first);
    r_vec->resize(s.first);
        
    ValType alpha = 0.0, beta = 0.0, omega = 0.0;
    
    FESystemUInt n_iters=0;
        
    this->krylov_basis_mat->zero();
    this->eig_val_vec->zero();
    
    q_vec->initializeToRandomUnitVector(0.0, 1.0);
    q_vec_old->zero();

    FESystem::EigenSolvers::QRMethodLinearEigenSolver<ValType> qr_eigen_solver;
    qr_eigen_solver.setEigenProblemType(this->getEigenProblemType());

    switch (this->getEigenProblemType()) {
        case FESystem::EigenSolvers::HERMITIAN:
        {    
                        
            FESystemBoolean convergence = false;
            
            while (!convergence)
            {
                std::cout << "Iteration #: " << n_iters << std::endl;
                // place the unit vector in the Krylov basis matrix
                this->krylov_basis_mat->setColumnVals(n_iters, 0, s.first-1, *q_vec);
//                if (n_iters > 1)
//                {
//                    FESystem::EigenSolvers::ModifiedQRFactorization<FESystemDouble> mqr;
//                    mqr.setMatrix(this->krylov_basis_mat.get());
//                    mqr.factorize();
//                    this->krylov_basis_mat->copyMatrix(mqr.getQMatrix());
//                    this->krylov_basis_mat->getColumnVals(n_iters, 0, s.first-1, *q_vec);
//                }
                
                this->A_mat->rightVectorMultiply(*q_vec, *r_vec); // r = A q 
                alpha = r_vec->dotProduct(*q_vec);  // alpha = r^T q
                
                // copy this value to the hessenberg matrix
                this->hessenberg_mat->setVal(n_iters, n_iters, alpha);

                // calculate the eigenvalue approximation and look for convergence
                tri_diag_mat->resize(n_iters+1, n_iters+1);
                tri_diag_mat->setSubMatrixVals(0, n_iters, 0, n_iters,
                                               0, n_iters, 0, n_iters,
                                               *(this->hessenberg_mat));
                
                std::cout << "**** H **** " << std::endl;
                tri_diag_mat->write(std::cout);
                
                // get eigenvalue of tridiagonal matrix
                qr_eigen_solver.setMatrix(tri_diag_mat.get());
                qr_eigen_solver.solve();

//                // calculate the convergence in the eigenvalues
//                eig_val_convergence->resize(n_iters+1);
//                eig_val_convergence->setSubVectorVals(0, n_iters, 
//                                                      0, n_iters,
//                                                      *(this->eig_val_vec));
//                
//                eig_val_convergence->add(-1.0, qr_eigen_solver.getEigenValues());
//                eig_val_convergence->setVal(n_iters, 0.0); /* the last eigenvalue is set to zero because the previous 
//                                                            * iteration had only n_iters-1 eigenvalue estimates */

                // copy the current eigenvalue estimates for further reference
                this->eig_val_vec->setSubVectorVals(0, n_iters, 
                                                    0, n_iters, 
                                                    qr_eigen_solver.getEigenValues());
                
                std::cout << "**** EIG **** " << std::endl;
                this->eig_val_vec->write(std::cout);
                
                if (n_iters == (s.first-1))
                    break;
                
                r_vec->add(-alpha, *q_vec);
                r_vec->add(-beta, *q_vec_old);
                
                beta = r_vec->getL2Norm(); std::cout << beta << std::endl;
                
                // set the value in the hessenberg matrix
                this->hessenberg_mat->setVal(n_iters, n_iters+1, beta);
                this->hessenberg_mat->setVal(n_iters+1, n_iters, beta);

                r_vec->scale(1.0/beta);  // q = r / beta
                q_vec_old->copyVector(*q_vec);
                q_vec->copyVector(*r_vec);
                r_vec->zero();

                // now conduct the various checks on the eigenvalue convergence, and purity of basis
                // check for converged eigenvalues
                for (FESystemUInt i=0; i<n_iters; i++)
                {
                    
                }
                
                // value of beta should be non-singular
                if (fabs(beta) < MACHINE_EPSILON)
                {
                    // restart the Lanczos process with a new starting vector
                }
                    
                // repeated eigenvalues
                
                // orthogonality loss
                // update the orthogonalizaiton error estimates 
                orthogonalization_factors->setVal(n_iters,n_iters,1.0);
                if (n_iters > 1) // lower than this is always orthogonal as it contains only 3 vectors
                    if (true)
                    {
                        krylov_basis_mat->matrixTransposeRightMultiply (1.0,
                                                                        *krylov_basis_mat,
                                                                        *orthogonalization_factors);
                    }
                else
                {
                    // the recursion formula is 
                    // beta_{j+1,j} omega_{k,j+1} = beta_{k+1,k} omega_{j,k+1} + (alpha_{k,k}-alpha_{j,j}) omega_{k,j} + beta_{k-1,k} omega_{j,k-1} - beta_{j-1,j} omega_{k,j-1} - errors
                    for (FESystemUInt k=0;k<(n_iters-1);k++)
                    {
                        omega = hessenberg_mat->getVal(k+1,k) * orthogonalization_factors->getVal(n_iters,k+1); // beta_{k+1,k} omega_{j,k+1} 
                        if (k>0)
                            omega += hessenberg_mat->getVal(k-1,k) * orthogonalization_factors->getVal(n_iters,k-1); // beta_{k-1,k} omega_{j,k-1} 
                        omega -= hessenberg_mat->getVal(n_iters-1,n_iters) * orthogonalization_factors->getVal(n_iters-1,k); // beta_{j-1,j} omega_{j-1,k}
                        omega += (hessenberg_mat->getVal(k,k)-hessenberg_mat->getVal(n_iters,n_iters)) * orthogonalization_factors->getVal(n_iters,k); // (alpha_{k,k}-alpha_{j,j}) omega_{j,k} 
                        
                        omega /= beta; // divide by omega_{j+1,j}
                        omega += ((1.0 * rand()) / (1.0 * RAND_MAX) * 2.0-1.0)*MACHINE_EPSILON; // add the random perturbation which is a uniformly distributed parameter between -epsilon and +epsilon
                        orthogonalization_factors->setVal(n_iters+1,k,omega); 
                    }
                }
                std::cout << "**** OMEGA **** " << std::endl;
                orthogonalization_factors->write(std::cout);
                                
                n_iters++; // increment the iteration number 
                
            }
            
            this->krylov_basis_mat->matrixRightMultiply(1.0,qr_eigen_solver.getEigenVectorMatrix(), *(this->eig_vec_mat));
            this->n_krylov_basis = n_iters; 
        }
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
            break;
    }
    
}




/***************************************************************************************/
// Template instantiations for some generic classes

INSTANTIATE_CLASS_FOR_ONLY_REAL_DATA_TYPES(FESystem::EigenSolvers::LanczosIterationLinearEigenSolver);

/***************************************************************************************/

