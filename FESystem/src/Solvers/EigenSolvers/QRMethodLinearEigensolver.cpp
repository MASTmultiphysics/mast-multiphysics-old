//
//  QRMethodLinearEigensolver.cpp
//  FESystem
//
//  Created by Manav Bhatia on 7/10/11.
//  Copyright 2011 . All rights reserved.
//


// FESystem includes 
#include "Solvers/EigenSolvers/QRMethodLinearEigensolver.h"
#include "Solvers/Factorizations/HouseHolderTriangulation.h"
#include "Solvers/Factorizations/TriangularBacksubstitution.h"
#include "Numerics/VectorBase.h"
#include "Numerics/MatrixBase.h"
#include "Base/macros.h"


template <typename ValType> 
FESystem::EigenSolvers::QRMethodLinearEigenSolver<ValType>::QRMethodLinearEigenSolver():
FESystem::EigenSolvers::LinearEigenSolverBase<ValType>()
{
    
}


template <typename ValType> 
FESystem::EigenSolvers::QRMethodLinearEigenSolver<ValType>::~QRMethodLinearEigenSolver()
{
    
}



//template <typename ValType> 
//void
//FESystem::EigenSolvers::QRMethodLinearEigenSolver<ValType>::setShift(ValType v)
//{
//    this->solver_shift = v;
//}



template <typename ValType> 
void
FESystem::EigenSolvers::QRMethodLinearEigenSolver<ValType>::solve()
{
    FESystemAssert0(this->matrices_are_set, FESystem::EigenSolvers::MatrixNotSet);
    
    std::auto_ptr<FESystem::Numerics::MatrixBase<ValType> > 
    tmat1(FESystem::Numerics::MatrixCreate<ValType>(this->A_mat->getType())),
    tmat2(FESystem::Numerics::MatrixCreate<ValType>(this->A_mat->getType()));

    std::pair<FESystemUInt, FESystemUInt> s = this->A_mat->getSize();
    
    tmat1->resize(s.first, s.second);
    tmat2->resize(s.first, s.second);
    
    std::auto_ptr<FESystem::FactorizationSolvers::MatrixQRFactorizationBase<ValType> >
    qr_factorization(new FESystem::FactorizationSolvers::HouseholderTriangulation<ValType>());
    
    switch (this->getEigenProblemType()) {
        case FESystem::EigenSolvers::HERMITIAN:
        {    
            std::auto_ptr<FESystem::Numerics::VectorBase<ValType> > 
            vec1(FESystem::Numerics::VectorCreate<ValType>(FESystem::Numerics::LOCAL_VECTOR)),
            vec2(FESystem::Numerics::VectorCreate<ValType>(FESystem::Numerics::LOCAL_VECTOR));
            
            this->eig_vec_mat_right->setToIdentity();
            this->eig_val_vec->zero();

            vec1->resize(s.first);
            vec2->resize(s.first);

            tmat1->copyMatrix(*(this->A_mat));
            
            FESystemBoolean convergence = false;
            
            while (!convergence)
            {                
                qr_factorization->clear();
                qr_factorization->setMatrix(tmat1.get());
                qr_factorization->factorize();
                
                qr_factorization->getRMatrix().matrixRightMultiply(1.0, qr_factorization->getQMatrix(),
                                                                   *tmat1); // tmat1 = R Q 
                
                this->eig_vec_mat_right->matrixRightMultiply(1.0, qr_factorization->getQMatrix(),
                                                       *tmat2);  // H = T3  Q 
                this->eig_vec_mat_right->copyMatrix(*tmat2);
                
                // check for convergence
//                std::cout << "**** RQ mat ****" << std::endl; tmat1->write(std::cout);
                tmat1->getDiagonal(*vec1);
                vec2->copyVector(*vec1);
                vec2->add(-1.0, *(this->eig_val_vec));
                std::cout << vec2->getL2Norm() << std::endl;
                //if (vec2->getL2Norm() < this->getConvergenceTolerance())
                if (vec2->getL2Norm() < 50.0)
                    convergence = true;

                // copy the final eigenvalue to vector
                this->eig_val_vec->copyVector(*vec1);
            }
            this->n_converged_eig_vals = s.first;
        }
            break;
            
        case FESystem::EigenSolvers::NONHERMITIAN:
        {    
            this->eig_vec_mat_right_complex->setToIdentity();
            this->eig_val_vec_complex->zero();
            
            std::auto_ptr<FESystem::Numerics::VectorBase<typename ComplexOperationType(ValType)> > 
            vec1(FESystem::Numerics::VectorCreate<typename ComplexOperationType(ValType)>(FESystem::Numerics::LOCAL_VECTOR)),
            vec2(FESystem::Numerics::VectorCreate<typename ComplexOperationType(ValType)>(FESystem::Numerics::LOCAL_VECTOR));

            std::auto_ptr<FESystem::Numerics::VectorBase<ValType> > 
            vec3(FESystem::Numerics::VectorCreate<ValType>(FESystem::Numerics::LOCAL_VECTOR));
            
            vec1->resize(s.first);
            vec2->resize(s.first);
            vec3->resize(s.first);
            
            tmat1->copyMatrix(*(this->A_mat));
            
            FESystemBoolean convergence = false;
            std::pair<FESystemUInt,FESystemUInt> s = this->A_mat->getSize();
            
            typename RealOperationType(ValType) threshold = this->A_mat->getFrobeniusNorm()*sqrt(FESystem::Base::getMachineEpsilon<typename RealOperationType(ValType)>());
            typename ComplexOperationType(ValType) eig1, eig2;
            
            while (!convergence)
            {                
                qr_factorization->clear();
                qr_factorization->setMatrix(tmat1.get());  // Q R = T1
                qr_factorization->factorize();
                
                qr_factorization->getRMatrix().matrixRightMultiply(1.0, qr_factorization->getQMatrix(),
                                                                   *tmat1); // T1 = R Q 
                
//                this->eig_vec_mat_right->matrixRightMultiply(1.0, qr_factorization->getQMatrix(),
//                                                       *tmat2);  // T3 = T3  Q 
//                this->eig_vec_mat_right->copyMatrix(*tmat2);
                
                // check for convergence
                std::cout << "**** RQ mat ****" << std::endl; tmat1->write(std::cout);
                
                tmat1->getDiagonal(*vec3); // the eigenvalue approaximations from the diagonal, which will be updated using block-diagonals
                for (FESystemUInt i=0; i<s.first; i++) 
                    vec1->setVal(i, vec3->getVal(i));  // copy to the complex vector
                
                // look at the sub-diagonals for convergence in the eigenvalues of the matrix blocks
                for (FESystemUInt i=0; i<s.first; ) 
                {
                    if ((i+1) < s.second)  // check subdiagonal entries till the third last entry
                    {
                        // if the entry is smaller than the threshold, use the update the eigenvalue 
                        // the last 2x2 matrix block is always handled through its eigenvalues
                        if ((i == (s.second-2)) || (FESystem::Base::magnitude<ValType, typename RealOperationType(ValType)>(tmat1->getVal(i+2,i)) < threshold))
                        {                                                               
                            // if the subdiagonal entry is less than the threshold, then the matrix diagonal 
                            // is a valid eigenvalue approximation, else use the 2x2 diagonal-block
                            if ((i == (s.second-2)) || (FESystem::Base::magnitude<ValType, typename RealOperationType(ValType)>(tmat1->getVal(i+1,i)) > threshold))
                            {
                                tmat1->getDiagonalBlockEigenvalue(i, eig1, eig2);
                                vec1->setVal(i,eig1);
                                vec1->setVal(i+1,eig2);
                                i+=2; // 2x2 block has been tested, hence increment two
                            }
                            else
                                i++; // only the diagonal entry was used, hence increment one
                        }
                        else
                            i++; // go to the next row
                    }
                    else
                        i++; // go to the next row
                }
                
                vec2->copyVector(*vec1);
                vec2->add(-1.0, *(this->eig_val_vec_complex));
                std::cout << "**** Error Norm ****" << std::endl; std::cout << vec2->getL2Norm() << std::endl;
                if (vec2->getL2Norm() < this->getConvergenceTolerance())
                    convergence = true;
                
                // copy the final eigenvalue to vector
                this->eig_val_vec_complex->copyVector(*vec1); 
                std::cout << "**** EIG VAL ****" << std::endl; this->eig_val_vec_complex->write(std::cout);
            }
            this->n_converged_eig_vals = s.first;
        }
            break;

        default:
            FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
            break;
    }
    
}




/***************************************************************************************/
// Template instantiations for some generic classes

INSTANTIATE_CLASS_FOR_ONLY_REAL_DATA_TYPES(FESystem::EigenSolvers::QRMethodLinearEigenSolver);

/***************************************************************************************/
