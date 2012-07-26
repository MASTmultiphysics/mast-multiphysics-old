//
//  QRMethodLinearEigensolver.h
//  FESystem
//
//  Created by Manav Bhatia on 7/10/11.
//  Copyright 2011 . All rights reserved.
//


#ifndef __fesystem_qr_method_linear_eigensolver_h__
#define __fesystem_qr_method_linear_eigensolver_h__

// C++ includes
#include <memory>

// FESystem includes
#include "Solvers/EigenSolvers/LinearEigenSolverBase.h"

namespace FESystem
{
    // Forward declerations
    namespace Numerics {template <typename ValType> class MatrixBase;}

    namespace Solvers
    {
        /*!
         *    This class uses the QR algorithm for solution of given eigenproblems. The algorithm 
         *    calculates all eigenvectors and eigenvalues for the system. 
         *
         *     \todo Improve the efficiency by application of shifts, and locking of converged eigenvalues.
         *     \todo The eigenvector and eigenvalues are currently only real numbers, and consideration of nonhermitian problems will need more careful implementation. 
         */
        
        template <typename ValType> 
        class QRMethodLinearEigenSolver: 
        public FESystem::Solvers::LinearEigenSolverBase<ValType>
        {
        public:
            
            QRMethodLinearEigenSolver();
            
            virtual ~QRMethodLinearEigenSolver();
            
            //void setShift(ValType v);
            
            virtual void solve();
            
        protected:
                                    
        };
        
    }
}



#endif // __fesystem_qr_method_linear_eigensolver_h__

