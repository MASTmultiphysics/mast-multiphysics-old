
// C++ includes
#include <iostream>
#include <memory>
#include <string>

// FESystem includes
#include "Base/FESystemBase.h"
#include "Base/FESystemTypes.h"
#include "Utils/CommandLineArguments.h"
#include "Numerics/DenseMatrix.h"
#include "Numerics/LocalVector.h"
#include "Solvers/EigenSolvers/RayleighQuotientIterationLinearEigensolver.h"
#include "Solvers/EigenSolvers/InverseIterationLinearEigensolver.h"
#include "Solvers/EigenSolvers/QRMethodLinearEigensolver.h"
#include "Solvers/EigenSolvers/PowerMethodLinearEigensolver.h"
#include "Solvers/EigenSolvers/LanczosIterationLinearEigensolver.h"
#include "Solvers/Factorizations/ModifiedQRFactorization.h"
#include "Solvers/Factorizations/ClassicalQRFactorization.h"
#include "Solvers/Factorizations/HouseHolderTriangulation.h"
#include "Solvers/Factorizations/HessenbergFormReduction.h"
#include "Plotting/PLPlot.h"

int eigen_main (int argc, char * const argv[]) {
    
	std::auto_ptr<FESystem::Base::InitializerBase> 
	fes_library(FESystem::Base::createInitializer(FESystem::Base::FESystem)); 
	
	
	FESystem::Utility::CommandLineArguments cargs;
	std::string str;
	str = "-np 1 ";
	cargs.appendArgument(str);
    
	fes_library->setCommandArguments(cargs);
	fes_library->initialize();

    FESystemUInt mat_size = 3;
    
    std::auto_ptr<FESystem::Numerics::MatrixBase<FESystemDouble> > 
    m(FESystem::Numerics::MatrixCreate<FESystemDouble>(FESystem::Numerics::LOCAL_DENSE_MATRIX).release()),
    tmat1(FESystem::Numerics::MatrixCreate<FESystemDouble>(FESystem::Numerics::LOCAL_DENSE_MATRIX).release()),
    tmat2(FESystem::Numerics::MatrixCreate<FESystemDouble>(FESystem::Numerics::LOCAL_DENSE_MATRIX).release());
    tmat1->resize(mat_size,mat_size);
    tmat2->resize(mat_size,mat_size);
    
    m->resize(mat_size,mat_size);
    FESystemUInt k=1;
    for (FESystemUInt i=0; i<m->getSize().first; i++)
        for (FESystemUInt j=0; j<m->getSize().second; j++)    
            //m->setVal(i, j, i+j+1);  // this is for symmetric matrix
            m->setVal(i,j,k++);
    for (FESystemUInt i=0; i<m->getSize().first; i++)
        for (FESystemUInt j=0; j<=i; j++) 
            //m->setVal(i, j, i+j+1);  // this is for symmetric matrix
            m->setVal(i,j,pow(k++,2));
    
//    m->zero(); for (FESystemUInt i=0; i<m->getSize().first; i++) m->setVal(i, i, pow(i,4)); 
//    m->resize(2,2); m->setVal(0,1,1); m->setVal(1, 0, 1); // this is a bad matrix, and the normal QR algorithm woudl never converge

    std::cout << "****  Matrix: ****" << std::endl;
    m->write(std::cout);
    std::cout << "Determinant: " << std::endl << m->getDeterminant() << std::endl;
    m->getInverse(*tmat1);
    std::cout << "Inverse: " << std::endl;
    tmat1->write(std::cout);
    
//    FESystem::EigenSolvers::QRMethodLinearEigenSolver<FESystemDouble> eig;
//    //eig.setShift(.000);
//    eig.setEigenProblemType(FESystem::EigenSolvers::NONHERMITIAN);
//    eig.setMatrix(m.get());
//    
//    eig.solve();
//    
//    std::cout << "****  Eigenvalues: ****" << std::endl;
//    eig.getEigenValues().write(std::cout);
//
//    std::cout << "****  Eigenvectors: ****" << std::endl;
//    eig.getEigenVectorMatrix().write(std::cout);

    
    return 0;
}

