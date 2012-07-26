//
//  EigenSolverBase.h
//  FESystem
//
//  Created by Manav Bhatia on 3/26/11.
//  Copyright 2011 . All rights reserved.
//

#ifndef __fesystem_linear_eigen_solver_base_h__
#define __fesystem_linear_eigen_solver_base_h__

// C++ includes
#include <memory>
#include <vector>

// FESystem includes
#include "Base/FESystemExceptions.h"



namespace FESystem
{
    // Forward declerations
    namespace Numerics {template <typename ValType> class VectorBase;}
    namespace Numerics {template <typename ValType> class MatrixBase;}

    namespace Solvers
    {

        enum LinearEigenProblemType
        {
            HERMITIAN,
            GENERALIZED_HERMITIAN, 
            NONHERMITIAN,
            GENERALIZED_NONHERMITIAN 
        };

        
        enum EigenSpectrumType
        {
            SMALLEST_MAGNITUDE,
            LARGEST_MAGNITUDE, 
            SMALLEST_REAL,
            LARGEST_REAL,
            SMALLEST_IMAGINARY,
            LARGEST_IMAGINARY
        };

        enum EigenShiftType
        {
            NO_SHIFT,
            SHIFT_AND_INVERT,
            CAYLEY_SHIFT,
            SPECTRUM_FOLD,
            ORIGIN_SHIFT
        };

	enum EigenValueSortingCriteria
	{
	  VALUE,
	  MAGNITUDE,
	  REAL_PART,
	  IMAGINARY_PART
	};

        
        template <typename ValType> 
        class LinearEigenSolverBase
        {
        public:
            
            LinearEigenSolverBase();
            
            virtual ~LinearEigenSolverBase();
            
            void setEigenProblemType(FESystem::Solvers::LinearEigenProblemType p);
            
            FESystem::Solvers::LinearEigenProblemType getEigenProblemType() const;

            void setEigenSpectrumType(FESystem::Solvers::EigenSpectrumType p);
            
            FESystem::Solvers::EigenSpectrumType getEigenSpectrumType() const;

            void setEigenShiftType(FESystem::Solvers::EigenShiftType p);
            
            FESystem::Solvers::EigenShiftType getEigenShiftType() const;

            void setEigenShiftValue(ValType v);
            
            ValType getEigenShiftValue() const;

            void setMaxIterations(FESystemUInt n);
            
            FESystemUInt getMaxIterations() const;

            void setMatrix(FESystem::Numerics::MatrixBase<ValType>* lhs, 
                           FESystem::Numerics::MatrixBase<ValType>* rhs = NULL);

            FESystem::Numerics::MatrixBase<ValType>& getAMatrix();

            FESystem::Numerics::MatrixBase<ValType>& getBMatrix();

            FESystemUInt getNConvergedEigenValues() const;
            
            void setConvergenceTolerance(FESystemDouble& val);

            FESystemDouble getConvergenceTolerance() const;

            /*!
             *  \brief returns the real valued eigenvector matrix for Hermitian problem
             */
            const FESystem::Numerics::MatrixBase<ValType>& getEigenVectorMatrix() const;
            
            const FESystem::Numerics::MatrixBase<FESystemComplexDouble>& getComplexEigenVectorMatrix() const;
            
            const FESystem::Numerics::VectorBase<ValType>& getEigenValues() const;

            const FESystem::Numerics::VectorBase<FESystemComplexDouble>& getComplexEigenValues() const;

	    /*!
	     *    Prepares the vector with ids of the eigenvalues in ascending magnitude. Can be used to access the 
	     *    lowest or highest eigenvalues in the spectrum
	     */
	    void  prepareSortingVector(FESystem::Solvers::EigenValueSortingCriteria sort, std::vector<FESystemUInt>& sorted_ids) const;

            virtual void solve() = 0;
            
        protected:

            void completePowerIterations(FESystem::Numerics::MatrixBase<ValType>& mat,
                                         FESystem::Numerics::VectorBase<ValType>& eig_vals,
                                         FESystem::Numerics::MatrixBase<ValType>& eig_vecs);

            FESystemUInt n_max_iterations;
            
            FESystem::Solvers::LinearEigenProblemType problem_type;

            FESystem::Solvers::EigenShiftType shift_type;
            
            ValType shift_value;

            FESystem::Solvers::EigenSpectrumType spectrum;
            
            FESystemBoolean matrices_are_set;
            
            /*!
             *  \brief stores the number of converged eigenvalues
             */
            FESystemUInt n_converged_eig_vals;
            
            /*!
             *  \brief required convergence tolerance, default value is 1.0e-12;
             */
            FESystemDouble convergence_tolerance; 

            FESystem::Numerics::MatrixBase<ValType>* A_mat;
            
            FESystem::Numerics::MatrixBase<ValType>* B_mat;
            
            /*!
             *  \brief storage for the eigenvectors once they are calculated after the iterations, this is used only for hermitian problems
             */
            std::auto_ptr<FESystem::Numerics::MatrixBase<ValType> > eig_vec_mat; 

            /*!
             *  \brief storage for the eigenvectors once they are calculated after the iterations, this is used only for nonhermitian problems
             */
            std::auto_ptr<FESystem::Numerics::MatrixBase<FESystemComplexDouble> > eig_vec_mat_complex; 

            /*!
             *  \brief storage for the eigenvalues as they are updated during the iterations, this is used only for hermitian problems
             */
            std::auto_ptr<FESystem::Numerics::VectorBase<ValType> > eig_val_vec;
            
            /*!
             *  \brief storage for the complex eigenvalues as they are updated during the iterations, this is specifically used for nonhermitian problems
             */
            std::auto_ptr<FESystem::Numerics::VectorBase<FESystemComplexDouble> > eig_val_vec_complex;

        };
        
        DeclareException0(GeneralizedEigenproblemDoesNotHaveRHSMatrix, 
                          << "RHS matrix does not exist for Hermitian and NonHermitian Problems\n");

        DeclareException0(MatrixNotSet, 
                          << "Matrices not set for eigensolver before use.\n");

    }
}


#endif // __fesystem_linear_eigen_solver_base_h__

