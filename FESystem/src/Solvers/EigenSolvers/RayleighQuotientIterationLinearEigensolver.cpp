//
//  RayleighQuotientIterationEigensolver.cpp
//  FESystem
//
//  Created by Manav Bhatia on 5/9/11.
//  Copyright 2011 . All rights reserved.
//

// C++ inlcudes
#include <stdlib.h>

// FESystem includes
#include "Solvers/EigenSolvers/RayleighQuotientIterationLinearEigensolver.h"
#include "Solvers/Factorizations/HouseHolderTriangulation.h"
#include "Solvers/Factorizations/TriangularBacksubstitution.h"
#include "Numerics/VectorBase.h"
#include "Numerics/MatrixBase.h"
#include "Base/macros.h"


template <typename ValType> 
FESystem::Solvers::RayleighQuotientIterationLinearEigenSolver<ValType>::RayleighQuotientIterationLinearEigenSolver():
FESystem::Solvers::LinearEigenSolverBase<ValType>(),
solver_shift(0.0)
{
    
}


template <typename ValType> 
FESystem::Solvers::RayleighQuotientIterationLinearEigenSolver<ValType>::~RayleighQuotientIterationLinearEigenSolver()
{
    
}



template <typename ValType> 
void
FESystem::Solvers::RayleighQuotientIterationLinearEigenSolver<ValType>::setShift(ValType v)
{
    this->solver_shift = v;
}


template <typename ValType> 
void
FESystem::Solvers::RayleighQuotientIterationLinearEigenSolver<ValType>::solve()
{
    FESystemAssert0(this->matrices_are_set, FESystem::Solvers::MatrixNotSet);
    
    
    switch (this->getEigenProblemType()) {
        case FESystem::Solvers::HERMITIAN:
        case FESystem::Solvers::NONHERMITIAN:    
        {    
            
            this->shiftAndInvertPowerIterations();
        }
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
            break;
    }
    
}


template <typename ValType> 
void
FESystem::Solvers::RayleighQuotientIterationLinearEigenSolver<ValType>::shiftAndInvertPowerIterations
()
{
    this->shifted_mat.reset(FESystem::Numerics::MatrixCreate<ValType>(this->getAMatrix().getType()).release());
    this->shifted_mat->resize(this->getAMatrix().getSize().first,
                              this->getAMatrix().getSize().second);
    
    FESystemUInt n = this->shifted_mat->getSize().second;
    
    std::auto_ptr<FESystem::Numerics::VectorBase<ValType> > 
    vec1(FESystem::Numerics::VectorCreate<ValType>(FESystem::Numerics::LOCAL_VECTOR).release()),
    vec2(FESystem::Numerics::VectorCreate<ValType>(FESystem::Numerics::LOCAL_VECTOR).release()),
    vec3(FESystem::Numerics::VectorCreate<ValType>(FESystem::Numerics::LOCAL_VECTOR).release());
    
    vec1->resize(n);
    vec2->resize(n);
    vec3->resize(n);
    
    // only the first eigenvalue is calculated
    // assuem an initial unit vector 
    for (FESystemUInt i=0; i<n; i++)
        vec1->setVal(i, 1.0);
    
    FESystemDouble conv = 1.0e6;
    ValType eig0 = this->solver_shift, eig1 = this->solver_shift, val;
    FESystem::Solvers::HouseholderTriangulation<ValType> qr_householder;
    FESystem::Solvers::TriangularBacksubstitution<ValType> back_substitute;

    back_substitute.setTriangularMatrixType(FESystem::Solvers::UPPER_TRIANGULAR);
    qr_householder.setMatrix(this->shifted_mat.get());
    
    while (fabs(conv) >= this->getConvergenceTolerance())
    {
        this->shifted_mat->copyMatrix(this->getAMatrix());
        this->shifted_mat->shiftDiagonal(-eig1);
                
        // calculate the QR Factorization of this matrix, and used that during the inversion operations
        qr_householder.clear();
        
        qr_householder.setMatrix(this->shifted_mat.get());
        qr_householder.factorize();
        
        FESystem::Numerics::MatrixBase<ValType>& q_mat = qr_householder.getQMatrix();
        FESystem::Numerics::MatrixBase<ValType>& r_mat = qr_householder.getRMatrix();
                
        back_substitute.setMatrix(r_mat);
        
        vec1->scaleToUnitLength();
        q_mat.leftVectorMultiply(*vec1, *vec2); // v2 = Q^T * v1;
        back_substitute.backSubstitute(*vec2, *vec3); // this solves Q R v3 = v1 
        
        val = 1.0; val /= vec3->dotProduct(*vec1);
        eig0 += val; // (v3' (A- lambda I)^(-1) v1)
        vec1->copyVector(*vec3);
        
        // convergence metric
        conv = fabs(abs(eig0 - eig1));
        eig1 = eig0; 
    }
    
    this->eig_val_vec->setVal(this->n_converged_eig_vals, eig1); // copy the eigenvalue
    vec1->scaleToUnitLength();
    this->eig_vec_mat->setColumnVals(this->n_converged_eig_vals, 0, n-1, *vec1); // copy the eigenvector 
    this->n_converged_eig_vals++;
}



/***************************************************************************************/
// Template instantiations for some generic classes

INSTANTIATE_CLASS_FOR_ONLY_REAL_DATA_TYPES(FESystem::Solvers::RayleighQuotientIterationLinearEigenSolver);

/***************************************************************************************/
