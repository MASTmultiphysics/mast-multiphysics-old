//
//  VonKarmanStrain1D.cpp
//  FESystem
//
//  Created by Manav Bhatia on 8/13/12.
//
//

// FESystem includes
#include "Disciplines/Structure/VonKarmanStrain1D.h"
#include "Base/FESystemExceptions.h"
#include "Mesh/ElemBase.h"
#include "Numerics/DenseMatrix.h"
#include "Numerics/LocalVector.h"
#include "FiniteElems/FiniteElementBase.h"
#include "Quadrature/QuadratureBase.h"
#include "Geom/Point.h"
#include "Disciplines/Structure/ExtensionBar.h"
#include "Disciplines/Structure/LinearBeamElementBase.h"


FESystem::Structures::VonKarmanStrain1D::VonKarmanStrain1D():
FESystem::Structures::Structural1DElementBase(),
bar_elem(NULL),
beam_elem(NULL)
{
    
}


FESystem::Structures::VonKarmanStrain1D::~VonKarmanStrain1D()
{
    
}



void
FESystem::Structures::VonKarmanStrain1D::clear()
{
    FESystem::Structures::Structural1DElementBase::clear();
    this->bar_elem = NULL;
    this->beam_elem = NULL;
}



void
FESystem::Structures::VonKarmanStrain1D::initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, const FESystem::Quadrature::QuadratureBase& q_rule,
                                                    FESystem::Structures::ExtensionBar& bar, FESystem::Structures::LinearBeamElementBase& beam)
{
    FESystem::Structures::StructuralElementBase::initialize(elem, fe, q_rule, 0, 0, 0);
    this->bar_elem = &bar;
    this->beam_elem = &beam;
}



void
FESystem::Structures::VonKarmanStrain1D::calculateTangentStiffnessMatrix(const FESystem::Numerics::VectorBase<FESystemDouble>& bar_sol, FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    const FESystemUInt n = this->geometric_elem->getNNodes();;
    FESystemUInt dims=n;
    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    
    if (this->beam_elem->ifIncludeChordwiseStiffness())
        dims *= 2;

    FESystemAssert4(((s.first == dims) && (s.second== dims)), FESystem::Numerics::MatrixSizeMismatch, dims, dims, s.first, s.second);
    
    
    static FESystem::Numerics::DenseMatrix<FESystemDouble> B_mat, stress_mat, tmp_mat1, tmp_mat2, mat_tr_and_ch;
    static FESystem::Numerics::LocalVector<FESystemDouble> pt;
    stress_mat.resize(1,1); B_mat.resize(1, 2*n); tmp_mat1.resize(1, 2*n), tmp_mat2.resize(2*n, 2*n), mat_tr_and_ch.resize(2*n, 2*n);
    stress_mat.zero(); B_mat.zero(); tmp_mat1.zero(); tmp_mat2.zero();
    
    const std::vector<FESystem::Geometry::Point*>& q_pts = this->quadrature->getQuadraturePoints();
    const std::vector<FESystemDouble>& q_weight = this->quadrature->getQuadraturePointWeights();
    
    FESystemDouble jac=0.0;
    mat_tr_and_ch.zero(); mat.zero();
    this->bar_elem->getStressTensor(pt, bar_sol, stress_mat);
    
    for (FESystemUInt i=0; i<q_pts.size(); i++)
    {
        jac = this->finite_element->getJacobianValue(*(q_pts[i]));
        this->calculateOperatorMatrix(*(q_pts[i]), B_mat);
        
        stress_mat.matrixRightMultiply(1.0, B_mat, tmp_mat1); // sigma B
        B_mat.matrixTransposeRightMultiply(1.0, tmp_mat1, tmp_mat2); // B^T sigma B
        
        mat_tr_and_ch.add(q_weight[i]*jac, tmp_mat2);
    }
    
    mat.setSubMatrixVals(0, dims-1, 0, dims-1, 0, dims-1, 0, dims-1, mat_tr_and_ch);
}



void
FESystem::Structures::VonKarmanStrain1D::calculateOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat)
{
    const FESystemUInt n = this->geometric_elem->getNNodes();;
    const std::pair<FESystemUInt, FESystemUInt> s = B_mat.getSize();
    FESystemAssert4(((s.first == 1) && (s.second== 2*n)), FESystem::Numerics::MatrixSizeMismatch, 1, 2*n, s.first, s.second);
   
    static std::vector<FESystemUInt> derivatives_x(2),derivatives_y(2);
    derivatives_x[0] = 1; derivatives_x[1] = 0;
    derivatives_y[0] = 0; derivatives_y[1] = 1;
    static FESystem::Numerics::LocalVector<FESystemDouble> Nvec;
    Nvec.resize(n);
    B_mat.zero();
    
    Nvec.zero();
    this->finite_element->getShapeFunctionDerivativeForPhysicalCoordinates(derivatives_x, pt, Nvec);
    B_mat.setRowVals(0, 0,    n-1, Nvec); // dv/dx
    B_mat.setRowVals(2, n,  2*n-1, Nvec); // dw/dx
    
    Nvec.zero();
    this->finite_element->getShapeFunctionDerivativeForPhysicalCoordinates(derivatives_y, pt, Nvec);
    B_mat.setRowVals(1, 0,    n-1, Nvec); // dv/dy
    B_mat.setRowVals(3, n,  2*n-1, Nvec); // dw/dy
}


