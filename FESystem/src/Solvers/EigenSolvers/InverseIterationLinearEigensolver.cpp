//
//  InverseIterationLinearEigensolver.cpp
//  FESystem
//
//  Created by Manav Bhatia on 3/27/11.
//  Copyright 2011 . All rights reserved.
//


// FESystem includes 
#include "Solvers/EigenSolvers/InverseIterationLinearEigensolver.h"
#include "Solvers/Factorizations/HouseHolderTriangulation.h"
#include "Solvers/Factorizations/TriangularBacksubstitution.h"
#include "Numerics/VectorBase.h"
#include "Numerics/MatrixBase.h"
#include "Base/macros.h"


template <typename ValType> 
FESystem::EigenSolvers::InverseIterationLinearEigenSolver<ValType>::InverseIterationLinearEigenSolver():
FESystem::EigenSolvers::LinearEigenSolverBase<ValType>(),
solver_shift(0.0)
{

}


template <typename ValType> 
FESystem::EigenSolvers::InverseIterationLinearEigenSolver<ValType>::~InverseIterationLinearEigenSolver()
{
    
}



template <typename ValType> 
void
FESystem::EigenSolvers::InverseIterationLinearEigenSolver<ValType>::setShift(ValType v)
{
    this->solver_shift = v;
}

    
template <typename ValType> 
void
FESystem::EigenSolvers::InverseIterationLinearEigenSolver<ValType>::solve()
{
    FESystemAssert0(this->matrices_are_set, FESystem::EigenSolvers::MatrixNotSet);

    std::pair<FESystemUInt, FESystemUInt> s = this->getAMatrix().getSize();

    this->shifted_mat.reset(FESystem::Numerics::MatrixCreate<ValType>(this->getAMatrix().getType()).release());
    this->shifted_mat->resize(s.first, s.second);
    
    switch (this->getEigenProblemType()) {
        case FESystem::EigenSolvers::HERMITIAN:
        case FESystem::EigenSolvers::NONHERMITIAN:    
        {    
            this->shifted_mat->copyMatrix(this->getAMatrix());
            this->shifted_mat->shiftDiagonal(-this->solver_shift);

            this->shiftAndInvertPowerIterations(*(this->shifted_mat));
        }
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
            break;
    }
    
}


template <typename ValType> 
void
FESystem::EigenSolvers::InverseIterationLinearEigenSolver<ValType>::shiftAndInvertPowerIterations
(FESystem::Numerics::MatrixBase<ValType>& mat)
{
    FESystemUInt n = mat.getSize().second;
    
    std::auto_ptr<FESystem::Numerics::VectorBase<ValType> > 
    vec1(FESystem::Numerics::VectorCreate<ValType>(FESystem::Numerics::LOCAL_VECTOR).release()),
    vec2(FESystem::Numerics::VectorCreate<ValType>(FESystem::Numerics::LOCAL_VECTOR).release()),
    vec3(FESystem::Numerics::VectorCreate<ValType>(FESystem::Numerics::LOCAL_VECTOR).release());
    
    vec1->resize(n);
    vec2->resize(n);
    vec3->resize(n);
    this->eig_val_vec->zero();

    // calculate the QR Factorization of this matrix, and used that during the inversion operations
    FESystem::FactorizationSolvers::HouseholderTriangulation<ValType> qr_householder;
    qr_householder.setMatrix(&mat);
    qr_householder.factorize();
    
    FESystem::Numerics::MatrixBase<ValType>& q_mat = qr_householder.getQMatrix();
    FESystem::Numerics::MatrixBase<ValType>& r_mat = qr_householder.getRMatrix();

    FESystem::FactorizationSolvers::TriangularBacksubstitution<ValType> back_substitute;
    back_substitute.setTriangularMatrixType(FESystem::FactorizationSolvers::UPPER_TRIANGULAR);
    back_substitute.setMatrix(r_mat);
    
    
    // only the first eigenvalue is calculated
    // assuem an initial unit vector 
    for (FESystemUInt i=0; i<n; i++)
        vec1->setVal(i, 1.0);
    
    typename RealOperationType(ValType) conv = 1.0e6;
    ValType eig0 = 0.0, eig1 = 0.0;
    
    while (fabs(conv) >= this->getConvergenceTolerance())
    {
        vec1->scaleToUnitLength();
        q_mat.leftVectorMultiply(*vec1, *vec2); // v2 = Q^T * v1;
        back_substitute.backSubstitute(*vec2, *vec3); // this solves Q R v3 = v1 
        
        eig0 = 1.0;
        eig0 /= vec3->dotProduct(*vec1)+this->solver_shift; // (v3' (A- lambda I)^(-1) v1)
        vec1->copyVector(*vec3);
        
        // convergence metric
        conv =  FESystem::Base::magnitude<ValType, typename RealOperationType(ValType)>(eig0 - eig1);
        eig1 = eig0;         
    }
    
    this->eig_val_vec->setVal(this->n_converged_eig_vals, eig1); // copy the eigenvalue
    vec1->scaleToUnitLength(); // scale the eigenvector to unit length before copying
    this->eig_vec_mat_right->setColumnVals(this->n_converged_eig_vals, 0, n-1, *vec1); // copy the eigenvector

    this->n_converged_eig_vals++; // increment the counter for the converged eigenvalue    
}



/***************************************************************************************/
// Template instantiations for some generic classes

INSTANTIATE_CLASS_FOR_ONLY_REAL_DATA_TYPES(FESystem::EigenSolvers::InverseIterationLinearEigenSolver);


/***************************************************************************************/
