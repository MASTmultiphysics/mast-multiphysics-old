//
//  ArpackLinearEigenSolver.cpp
//  FESystem
//
//  Created by Manav Bhatia on 5/31/12.
//  Copyright (c) 2012. All rights reserved.
//

// FESystem includes
#include "Solvers/EigenSolvers/ArpackLinearEigenSolver.h"
#include "Solvers/LinearSolvers/LinearSolverBase.h"
#include "Numerics/MatrixBase.h"
#include "Numerics/LocalVector.h"




template <typename ValType>
FESystem::Solvers::ArpackLinearEigenSolver<ValType>::ArpackLinearEigenSolver():
FESystem::Solvers::LinearEigenSolverBase<ValType>::LinearEigenSolverBase(),
solution_completed(false),
if_initialized(false),
n_converged_eigen_pairs(0),
ido(0),
n(0),
nev(0),
tol(0.0),
ncv(0),
ldv(0),
lworkl(0),
info(0),
rvec(0),
ldz(0),
sigmar(0.0),
sigmai(0.0)
{
    
}




template <typename ValType>
FESystem::Solvers::ArpackLinearEigenSolver<ValType>::~ArpackLinearEigenSolver()
{
    this->clear();
}





template <typename ValType>
void
FESystem::Solvers::ArpackLinearEigenSolver<ValType>::clear()
{
    this->solution_completed = false;
    this->if_initialized = false;
    
    /// number of converged eigenpairs
    this->n_converged_eigen_pairs = 0;
    
    this->ido = 0;
    this->n = 0;
    this->nev = 0;
    this->tol = 0.0;
    this->ncv = 0;
    this->ldv = 0;
    this->lworkl = 0;
    this->info = 0;
    this->rvec = 0;
    this->ldz = 0;
    this->sigmar = 0.0;
    this->sigmai = 0.0;
    this->bmat.clear();
    this->which.clear();
    this->select.clear();
    this->HowMny.clear();
    this->Z.clear();
    this->dr.clear();
    this->di.clear();
    this->V.clear();
    this->resid.clear();
    this->iparam.clear();
    this->ipntr.clear();
    this->workd.clear();
    this->workl.clear();
    this->workev.clear();
    this->operator_matrix.reset();
    this->linear_solver = NULL;
}


template <typename ValType> 
void 
FESystem::Solvers::ArpackLinearEigenSolver<ValType>::setLinearSolver(FESystem::Solvers::LinearSolverBase<ValType>& solver)
{
    this->linear_solver = &solver;
}



template <typename ValType>
void 
FESystem::Solvers::ArpackLinearEigenSolver<ValType>::init(FESystemUInt n_eig_to_compute, FESystemBoolean if_calculate_eig_vec)
{
    // make sure that the solution has not been computed
    FESystemAssert0(!this->if_initialized, FESystem::Exception::InvalidState);
    FESystemAssert0(!this->solution_completed, FESystem::Exception::InvalidState);
    FESystemAssert0(this->matrices_are_set, FESystem::Solvers::MatrixNotSet);
    
    std::pair<FESystemUInt, FESystemUInt> s = this->A_mat->getSize();

    this->ido = 0;
    this->n = s.first;
    this->nev = n_eig_to_compute;
    this->tol = this->getConvergenceTolerance();
    // a higher than recommended value of ncv is kept. 
    this->ncv = 3 * this->nev+1;
    if (this->ncv > this->n) this->ncv = this->n;
    this->ldv = this->n;
    this->info = 0;
    if (if_calculate_eig_vec)
    {
        this->rvec = 1;
        this->HowMny = "A";
        // the size for hermitian is nev*n + 1, and for non-hermitian it 
        // is (nev+1)*n+1. Since both have the same number of rows, it is 
        // safe to use just the bigger of the two number of columns. 
        this->Z.resize((this->nev+1) * this->n + 1);
        std::fill(this->Z.begin(), this->Z.end(), 0.0);
    }
    else
        this->rvec = 0;
    // select should contain the selected eigen vectors to be computed, but we are computing 
    // all by default
    this->select.resize(this->ncv + 1);
    std::fill(this->select.begin(), this->select.end(), 0);
    this->ldz = this->n;
    
    
    switch (this->getEigenProblemType())
    {
        case FESystem::Solvers::HERMITIAN:
        case FESystem::Solvers::NONHERMITIAN:
            this->bmat = "I";
            break;
            
        case FESystem::Solvers::GENERALIZED_HERMITIAN:
        case FESystem::Solvers::GENERALIZED_NONHERMITIAN:
            this->bmat = "G";
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
            break;
    }
    
    
    switch (this->getEigenProblemType())
    {
        case FESystem::Solvers::HERMITIAN:
        case FESystem::Solvers::GENERALIZED_HERMITIAN:
        {
            this->lworkl = this->ncv * (this->ncv + 8);
            this->ipntr.resize(12);
            std::fill(this->ipntr.begin(), this->ipntr.end(), 0);
            
            
            switch(this->getEigenSpectrumType())
            {
                case FESystem::Solvers::LARGEST_MAGNITUDE:
                    this->which = "LM";
                    break;
                    
                case FESystem::Solvers::SMALLEST_MAGNITUDE:
                    this->which = "SM";
                    break;
                    
                case FESystem::Solvers::LARGEST_REAL:
                    this->which = "LA";
                    break;
                    
                case FESystem::Solvers::SMALLEST_REAL:
                    this->which = "SA";
                    break;
                    
                case FESystem::Solvers::LARGEST_IMAGINARY:
                case FESystem::Solvers::SMALLEST_IMAGINARY:
                default: 
                    FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
            }
        }
            break;
            
            
        case FESystem::Solvers::NONHERMITIAN:
        case FESystem::Solvers::GENERALIZED_NONHERMITIAN:
        {
            this->lworkl = 2*( this->ncv * (3 * this->ncv + 6));
            this->ipntr.resize(15);
            std::fill(this->ipntr.begin(), this->ipntr.end(), 0);
            this->workev.resize(3* this->ncv+1);
            std::fill(this->workev.begin(), this->workev.end(), 0.0);
            
            switch(this->getEigenSpectrumType())
            {
                case FESystem::Solvers::LARGEST_MAGNITUDE:
                    this->which = "LM";
                    break;
                    
                case FESystem::Solvers::SMALLEST_MAGNITUDE:
                    this->which = "SM";
                    break;
                    
                case FESystem::Solvers::LARGEST_REAL:
                    this->which = "LR";
                    break;
                    
                case FESystem::Solvers::SMALLEST_REAL:
                    this->which = "SR";
                    break;
                    
                case FESystem::Solvers::LARGEST_IMAGINARY:
                    this->which = "LI";
                    break;
                    
                case FESystem::Solvers::SMALLEST_IMAGINARY:
                    this->which = "SI";
                    break;
                    
                default: 
                    FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
            }
        }
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
            break;
    }
    
    
    this->dr.resize(this->nev+1+1);
    std::fill(this->dr.begin(), this->dr.end(), 0.0);
    this->di.resize(this->nev+1+1);
    std::fill(this->di.begin(), this->di.end(), 0.0);
    this->V.resize(this->ncv * this->n + 1);
    std::fill(this->V.begin(), this->V.end(), 0.0);
    this->resid.resize(this->n+1);
    std::fill(this->resid.begin(), this->resid.end(), 0.0);
    this->iparam.resize(12);
    std::fill(this->iparam.begin(), this->iparam.end(), 0);
    this->workd.resize(3 * this->n + 1);
    std::fill(this->workd.begin(), this->workd.end(), 0.0);
    this->workl.resize(this->lworkl + 1);   
    std::fill(this->workl.begin(), this->workl.end(), 0.0);
    
    //  tell the solver that exact shifts are used
    this->iparam[1] = 1;
    
    // set the max allowable iterations
    this->iparam[3] = this->getMaxIterations();
    
    // now set the problem shift type
    switch (this->getEigenShiftType())
    {
        case FESystem::Solvers::NO_SHIFT:
            if (this->bmat == "I")
                this->iparam[7] = 1;
            else 
                this->iparam[7] = 2;
            break;
            
        case FESystem::Solvers::SHIFT_AND_INVERT:
            this->iparam[7] = 3;
            break;
            
        case FESystem::Solvers::CAYLEY_SHIFT:
        {
            switch (this->getEigenProblemType())
            {
                case FESystem::Solvers::HERMITIAN:
                case FESystem::Solvers::GENERALIZED_HERMITIAN:
                    this->iparam[7] = 5;
                    break;
                    
                case FESystem::Solvers::NONHERMITIAN:
                case FESystem::Solvers::GENERALIZED_NONHERMITIAN:
                default:
                    FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
                    break;
            }
        }
            break;
            
        case FESystem::Solvers::SPECTRUM_FOLD:
        case FESystem::Solvers::ORIGIN_SHIFT:
        default: 
            FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
    }
    
    // only sigmar is used, since the imaginary part of the shift is not allowed
    if (this->iparam[7] > 2)
        this->sigmar = this->getEigenShiftValue();

    
    // initialize the eigensolver
    FESystemAssert0(this->linear_solver != NULL, FESystem::Exception::NULLQuantity);
    FESystemAssert0(this->operator_matrix.get() == NULL, FESystem::Exception::InvalidState);

    this->operator_matrix.reset(FESystem::Numerics::MatrixCreate<ValType>(this->getAMatrix().getType()).release());
    
    switch (this->getEigenShiftType())
    {
        case FESystem::Solvers::NO_SHIFT:
            if (this->bmat == "G")
            {
                this->operator_matrix->copyMatrix(this->getBMatrix());
                this->linear_solver->setSystemMatrix(*this->operator_matrix);
            }
            break;
            
        case FESystem::Solvers::SHIFT_AND_INVERT:
            switch (this->getEigenProblemType())
        {
            case FESystem::Solvers::HERMITIAN:
            case FESystem::Solvers::NONHERMITIAN:
            {
                this->operator_matrix->copyMatrix(this->getAMatrix());
                this->operator_matrix->shiftDiagonal(- this->sigmar);
                this->linear_solver->setSystemMatrix(*this->operator_matrix);
            }
                break;
                
            case FESystem::Solvers::GENERALIZED_HERMITIAN:
            case FESystem::Solvers::GENERALIZED_NONHERMITIAN:
            {
                this->operator_matrix->copyMatrix(this->getAMatrix());
                this->operator_matrix->add(-this->sigmar, this->getBMatrix());
                this->linear_solver->setSystemMatrix(*this->operator_matrix);
            }
                break;
                
            default:
                FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
                break;
        }
            break;
            
        case FESystem::Solvers::CAYLEY_SHIFT:
        case FESystem::Solvers::SPECTRUM_FOLD:
        case FESystem::Solvers::ORIGIN_SHIFT:
        default: 
            FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
    }

    this->if_initialized = true;
}




template <>
void 
FESystem::Solvers::ArpackLinearEigenSolver<FESystemDouble>::solve()
{
    FESystemAssert0(this->if_initialized, FESystem::Solvers::MatrixNotSet);
    FESystemAssert0(this->matrices_are_set, FESystem::Solvers::MatrixNotSet);
    
    FESystem::Numerics::MatrixBase<FESystemDouble>* solver_A_mat = NULL, *solver_B_mat = NULL;
    
    std::pair<FESystemUInt, FESystemUInt> s = this->A_mat->getSize();
    
    switch(this->getEigenProblemType())
    {
        case FESystem::Solvers::HERMITIAN:
        case FESystem::Solvers::NONHERMITIAN:
            solver_A_mat = &(this->getAMatrix());
            break;
            
        case FESystem::Solvers::GENERALIZED_HERMITIAN:
        case FESystem::Solvers::GENERALIZED_NONHERMITIAN:
            solver_A_mat = &(this->getAMatrix());
            solver_B_mat = &(this->getBMatrix());
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
            break;
    }
    
    
    // create the vectors for matrix operations
    FESystem::Numerics::LocalVector<FESystemDouble> vec1, vec2, vec3;
    
    vec3.resize(s.first);
    
    // solve 
    while (this->ido != 99)
    {
        switch (this->getEigenProblemType())
        {
                // this is for real symmetric
            case FESystem::Solvers::HERMITIAN:
            case FESystem::Solvers::GENERALIZED_HERMITIAN:
                dsaupd_(&ido, const_cast<char*>(bmat.c_str()), 
                        &n, const_cast<char*>(which.c_str()), 
                        &nev, &tol,
                        &resid[1], &ncv,
                        &V[1], &ldv,
                        &iparam[1], &ipntr[1],
                        &workd[1], &workl[1],
                        &lworkl, &info);
                break;
                
                // this is for real un-symmetric
            case FESystem::Solvers::NONHERMITIAN:
            case FESystem::Solvers::GENERALIZED_NONHERMITIAN:
                dnaupd_(&ido, const_cast<char*>(bmat.c_str()), 
                        &n, const_cast<char*>(which.c_str()), 
                        &nev, &tol,
                        &resid[1], &ncv,
                        &V[1], &ldv,
                        &iparam[1], &ipntr[1],
                        &workd[1], &workl[1],
                        &lworkl, &info);
                break;
                
            default:
                FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
                break;
        }
        
        
        // check the error, if any
        this->checkError(this->info);
        
        // copy the vector into vec1
        vec1.resize(s.first, &(workd[ipntr[1]]));
        vec2.resize(s.first, &(this->workd[ipntr[2]]));
        vec3.zero();
        
        
        // this is for real symmetric problems only
        switch (this->ido)
        {
            case -1:
            case 1:
            {
                // compute  Y = OP * X  for case -1
                // compute  Y = OP * X without performing the B * X operation in 3, 4 and 5 for case 1
                switch (this->iparam[7])
                {
                    case 1:
                        // OP = A, B = I
                        solver_A_mat->rightVectorMultiply(vec1, vec2);
                        break;
                        
                    case 2:
                        // OP = M^(-1) A , B = M
                        solver_A_mat->rightVectorMultiply(vec1, vec3);
                        this->linear_solver->solve(vec3, vec2);
                        break;
                        
                    case 3:
                    {
                        // OP = (A - sigma I)^(-1) , B = I
                        // OP = (A - sigma M)^(-1) M , B = M  (for generalized problem)
                        switch (this->getEigenProblemType())
                        {
                            case FESystem::Solvers::HERMITIAN:
                            case FESystem::Solvers::NONHERMITIAN:
                                this->linear_solver->solve(vec1, vec2);
                                break;
                                
                            case FESystem::Solvers::GENERALIZED_HERMITIAN:
                            case FESystem::Solvers::GENERALIZED_NONHERMITIAN:
                            {
                                if (this->ido == -1)
                                {
                                    solver_B_mat->rightVectorMultiply(vec1, vec3);
                                    this->linear_solver->solve(vec3, vec2);
                                }
                                else 
                                {
                                    vec1.resize(s.first, &(workd[ipntr[3]]));
                                    this->linear_solver->solve(vec1, vec2);
                                }
                            }
                                break;
                                
                            default:
                                FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
                                break;
                        }
                    }
                        break;
                        
                    case 4: // not handled for now
                    case 5: // not handled for now
                    default:
                        FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
                }
            }
                break;
                
            case 2:
            {
                // compute  Y = B * X 
                switch (this->iparam[7])
                {
                    case 1:
                        // OP = A, B = I
                        vec2.copyVector(vec1);
                        break;
                        
                    case 2:
                        // OP = M^(-1) A , B = M
                        solver_B_mat->rightVectorMultiply(vec1, vec2);
                        break;
                        
                    case 3:
                    {
                        // OP = (A - sigma I)^(-1) , B = I
                        // OP = (A - sigma M)^(-1) M , B = M  (for generalized problem)
                        switch (this->getEigenProblemType())
                        {
                            case FESystem::Solvers::HERMITIAN:
                            case FESystem::Solvers::NONHERMITIAN:
                                vec2.copyVector(vec1);
                                break;
                                
                            case FESystem::Solvers::GENERALIZED_HERMITIAN:
                            case FESystem::Solvers::GENERALIZED_NONHERMITIAN:
                                solver_B_mat->rightVectorMultiply(vec1, vec2);
                                break;
                                
                            default:
                                FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
                                break;
                        }
                    }
                        break;
                        
                    case 4: // not handled for now
                    case 5: // not handled for now
                    default:
                        FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
                }
                
            }
                break;
                                
            case 99:
                break;
                
            case 3:
            default:
                FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
        }
    }
    
    
    // once the solution is done, calculate the eigenvectors and eigenvalues
    switch (this->getEigenProblemType())
    {
            // this is for real symmetric
        case FESystem::Solvers::HERMITIAN:
        case FESystem::Solvers::GENERALIZED_HERMITIAN:
            dseupd_(&rvec, const_cast<char*>(HowMny.c_str()), 
                    &select[1], &dr[1],
                    &Z[1], &ldz, 
                    &sigmar, const_cast<char*>(bmat.c_str()),
                    &n, const_cast<char*>(which.c_str()),
                    &nev, &tol, 
                    &resid[1], &ncv, 
                    &V[1], &ldv, 
                    &iparam[1],
                    &ipntr[1], &workd[1], 
                    &workl[1], &lworkl, &info );
            break;
            
            // this is for real un-symmetric
        case FESystem::Solvers::NONHERMITIAN:
        case FESystem::Solvers::GENERALIZED_NONHERMITIAN:
        {
            // the Z vector is overwritten with the first NEV+1 vectors of 
            // V, since we need the ritz vectors of the problem, instead of the
            // the Arnoldi basis
            //      for (int j=0; j<(this->nev+1); j++)
            //        for (int i=0; i < this->n; i++)
            //          this->Z[ j*this->n + i + 1] = this->V[ j*this->n + i + 1];
            
            dneupd_(&rvec, const_cast<char*>(HowMny.c_str()), 
                    &select[1], &dr[1], &di[1],
                    &Z[1], &ldz, &sigmar, &sigmai, &workev[1],
                    const_cast<char*>(bmat.c_str()),
                    &n, const_cast<char*>(which.c_str()),
                    &nev, &tol, 
                    &resid[1], &ncv, 
                    &V[1], &ldv, 
                    &iparam[1],
                    &ipntr[1], &workd[1], 
                    &workl[1], &lworkl, &info );
        }
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
            break;
    }
    
    
    // check the error, if any
    this->checkError(this->info);
    
    // if everything is done and is complete, then the solution will be available 
    // in the vectors
    this->solution_completed = true;
    
}




template <typename ValType>
void 
FESystem::Solvers::ArpackLinearEigenSolver<ValType>::checkError(const unsigned int ierr)
{
    switch (ierr)
    {
        case 0:
        {
            //    =  0   : Normal exit.
            // nothing to be done, just return
            if (this->ido == 99)
            {
                // number of converged eigenvalues 
                this->n_converged_eig_vals = this->iparam[5];
                
                switch (this->getEigenProblemType())
                {
                        // this is for real symmetric
                    case HERMITIAN:
                    case GENERALIZED_HERMITIAN:
                    {
                        this->eig_val_vec->resize(this->n_converged_eig_vals);
                        this->eig_vec_mat->resize(this->n , this->n_converged_eig_vals);
                        
                        for (FESystemUInt i=0; i<this->n_converged_eig_vals; i++)
                        {
                            this->eig_val_vec->setVal(i, dr[i+1]);
                            for (FESystemUInt j=0; j < this->n; j++)
                                this->eig_vec_mat->setVal(j, i, this->Z[(i)*this->n+j+1]);
                        }
                    }
                        break;
                        
                        // this is for real un-symmetric
                    case NONHERMITIAN:
                    case GENERALIZED_NONHERMITIAN:
                    {
                        this->eig_val_vec_complex->resize(this->n_converged_eig_vals);
                        this->eig_vec_mat_complex->resize(this->n , this->n_converged_eig_vals);
                        
                        FESystemUInt index = 0;
                        for (FESystemUInt i=0; i<this->n_converged_eig_vals; i++)
                        {
                            this->eig_val_vec_complex->setVal(i, std::complex<ValType>(dr[i+1], di[i+1]));
                            
                            if (i % 2 ==0) 
                                index = i / 2;
                            else 
                                index = (i-1) / 2;
                            
                            for (int j=0; j < this->n; j++)
                                this->eig_vec_mat_complex->setVal(i, j, std::complex<ValType>(this->V[(2*index) * this->n+j+1], this->V[(2*index+1) * this->n+j+1]));
                        }
                    }
                        break;
                        
                    default:
                        FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
                        break;
                }
            }

        }
            break;
            
        case 1:
            //    =  1   : Maximum number of iterations taken. All possible
            //    eigenvalues of OP has been found. iparam[5]
            //    returns the number of converged Ritz values.
            std::cout << "Eigensolver reached maximum iterations." << std::endl;
            std::cout << "Number of converged eigenvalues: "
            << this->iparam[5] << std::endl;
            break;
            
        case 3:
            //    =  3   : No shifts could be applied during a cycle of the
            //    Implicitly restarted Arnoldi iteration. One
            //    possibility is to increase the size of NCV relative
            //    to nev. See remark 4 below.
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
            
        case -1:
            //    = -1   : n must be positive.
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
            
        case -2:
            //    = -2   : nev must be positive.
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
            
        case -3:
            //    = -3   : ncv must satisfy nev < ncv <= n.
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
            
        case -4:
            //    = -4   : The maximum number of Arnoldi update iterations allowed
            //    must be greater than zero.
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
            
        case -5:
            //    = -5   : which must be one of 'LM', 'SM', 'LA', 'SA' or 'BE'.
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
            
        case -6:
            //    = -6   : bmat must be one of 'I' or 'G'.
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
            
        case -7:
            //    = -7   : Length of private work array workl is not sufficient.
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
            
        case -8:
            //    = -8   : Error return from trid. eigenvalue calculation;
            //    Informational error from LAPACK routine dsteqr.
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
            
        case -9:
            //      = -9   : Starting vector is zero.
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
            
        case -10:
            //      = -10  : iparam[7] must be 1,2,3,4,5.
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
            
        case -11:
            //      = -11  : iparam[7] = 1 and bmat = 'G' are incompatible.
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
            
        case -12:
            //      = -12  : iparam[1] must be equal to 0 or 1.
            //    = -12: nev and which = 'BE' are incompatible.
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
            
        case -13:
            //      = -13  : nev and which = 'BE' are incompatible.
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
            
        case -14:
            //    = -14: dsaupp did not find any eigenvalues to sufficient
            //    accuracy.
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
            
        case -15:
            //    = -15: HowMny must be one of 'A' or 'S' if rvec = true.
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
            
        case -16:
            //      = -16: HowMny = 'S' not yet implemented.
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
            
        case -9999:
            //      = -9999: Could not build an Arnoldi factorization. iparam[5]
            //      returns the size of the current Arnoldi factorization.
            //      The user is advised to check that enough workspace
            //      and array storage has been allocated.
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
            
            
        default:
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
    }
}


/***************************************************************************************/
// Template instantiations for some generic classes

template class FESystem::Solvers::ArpackLinearEigenSolver<FESystemDouble> ;
//INSTANTIATE_CLASS_FOR_ONLY_REAL_DATA_TYPES(FESystem::Solvers::ArpackLinearEigenSolver);


/***************************************************************************************/

