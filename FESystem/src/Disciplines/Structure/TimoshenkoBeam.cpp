//
//  TimoshenkoBeam.cpp
//  FESystem
//
//  Created by Manav Bhatia on 8/26/12.
//
//

// FESystem includes
#include "Disciplines/Structure/TimoshenkoBeam.h"
#include "Base/FESystemExceptions.h"
#include "Mesh/ElemBase.h"
#include "Numerics/DenseMatrix.h"
#include "Numerics/LocalVector.h"
#include "FiniteElems/FiniteElementBase.h"
#include "Quadrature/QuadratureBase.h"
#include "Geom/Point.h"


FESystem::Structures::TimoshenkoBeam::TimoshenkoBeam():
FESystem::Structures::LinearBeamElementBase(),
quadrature_bending(NULL),
quadrature_shear(NULL),
kappa(5.0/6.0)
{
}


FESystem::Structures::TimoshenkoBeam::~TimoshenkoBeam()
{
    
}



void
FESystem::Structures::TimoshenkoBeam::clear()
{
    FESystem::Structures::LinearBeamElementBase::clear();
    this->quadrature_bending = NULL;
    this->quadrature_shear = NULL;
}





void
FESystem::Structures::TimoshenkoBeam::getStressTensor(const FESystem::Numerics::VectorBase<FESystemDouble>& pt, const FESystem::Numerics::VectorBase<FESystemDouble>& sol,
                                                            FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}




void
FESystem::Structures::TimoshenkoBeam::initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, const FESystem::Quadrature::QuadratureBase& q_bend,
                                                 const FESystem::Quadrature::QuadratureBase& q_shear, FESystemDouble E, FESystemDouble nu, FESystemDouble rho, FESystemDouble I_tr, FESystemDouble I_ch,
                                                 FESystemDouble A)
{
    FESystem::Structures::LinearBeamElementBase::initialize(elem, fe, q_bend, E, nu, rho, I_tr, I_ch, A);
    this->quadrature_bending = &q_bend;
    this->quadrature_shear = &q_shear;
}




void
FESystem::Structures::TimoshenkoBeam::calculateStiffnessMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);

    const FESystemUInt n = this->geometric_elem->getNNodes();
    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    
    static FESystem::Numerics::DenseMatrix<FESystemDouble> B_mat, B_mat_shear, C_mat_bend, C_mat_shear, tmp_mat1, tmp_mat2;
    
    FESystemAssert4(((s.first == 4*n) && (s.second == 4*n)), FESystem::Numerics::MatrixSizeMismatch, 4*n, 4*n, s.first, s.second);
    C_mat_bend.resize(2,2); C_mat_shear.resize(2,2); B_mat.resize(2, 4*n); B_mat_shear.resize(2, 4*n); tmp_mat1.resize(2, 4*n); tmp_mat2.resize(4*n, 4*n);
    
    C_mat_bend.zero(); C_mat_shear.zero(); B_mat.zero(); tmp_mat1.zero(); tmp_mat2.zero();
    
    
    const std::vector<FESystem::Geometry::Point*>& q_pts_bend = this->quadrature_bending->getQuadraturePoints();
    const std::vector<FESystemDouble>& q_weight_bend = this->quadrature_bending->getQuadraturePointWeights();
    
    const std::vector<FESystem::Geometry::Point*>& q_pts_shear = this->quadrature_shear->getQuadraturePoints();
    const std::vector<FESystemDouble>& q_weight_shear = this->quadrature_shear->getQuadraturePointWeights();
    
    FESystemDouble jac=0.0;
    mat.zero();
    this->getMaterialComplianceMatrix(C_mat_bend, C_mat_shear);
    
    // bending contribution
    for (FESystemUInt i=0; i<q_pts_bend.size(); i++)
    {
        jac = this->finite_element->getJacobianValue(*(q_pts_bend[i]));
        this->calculateBendingOperatorMatrix(*(q_pts_bend[i]), B_mat);
        
        C_mat_bend.matrixRightMultiply(1.0, B_mat, tmp_mat1);
        B_mat.matrixTransposeRightMultiply(1.0, tmp_mat1, tmp_mat2);
        
        mat.add(q_weight_bend[i]*jac, tmp_mat2);
    }
    
    // transverse shear contribution
    B_mat_shear.zero();
    for (FESystemUInt i=0; i<q_pts_shear.size(); i++)
    {
        jac = this->finite_element->getJacobianValue(*(q_pts_shear[i]));
        this->calculateShearOperatorMatrix(*(q_pts_shear[i]), B_mat_shear);
        
        C_mat_shear.matrixRightMultiply(1.0, B_mat_shear, tmp_mat1);
        B_mat_shear.matrixTransposeRightMultiply(1.0, tmp_mat1, tmp_mat2);
        
        mat.add(q_weight_shear[i]*jac, tmp_mat2);
    }
}





void
FESystem::Structures::TimoshenkoBeam::calculateBendingOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat)
{
    const FESystemUInt n = this->geometric_elem->getNNodes();
    const std::pair<FESystemUInt, FESystemUInt> s = B_mat.getSize();
    
    static std::vector<FESystemUInt> derivatives_x(1);
    derivatives_x[0] = 1;
    static FESystem::Numerics::LocalVector<FESystemDouble> Nvec;
    Nvec.resize(n);
    B_mat.zero();
    Nvec.zero();
    this->finite_element->getShapeFunctionDerivativeForPhysicalCoordinates(derivatives_x, pt, Nvec);
    
    FESystemAssert4(((s.first == 2) && (s.second == 4*n)), FESystem::Numerics::MatrixSizeMismatch, 2, 4*n, s.first, s.second);
    
    B_mat.setRowVals(0, 3*n,  4*n-1, Nvec); // epsilon-x: thetaz
    Nvec.scale(-1.0);
    B_mat.setRowVals(1, 2*n,  3*n-1, Nvec); // epsilon-x: thetay
}



void
FESystem::Structures::TimoshenkoBeam::calculateShearOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat)
{
    const FESystemUInt n = this->geometric_elem->getNNodes();
    const std::pair<FESystemUInt, FESystemUInt> s = B_mat.getSize();
    
    static std::vector<FESystemUInt> derivatives_x(1);
    derivatives_x[0] = 1;
    
    static FESystem::Numerics::LocalVector<FESystemDouble> Nvec;
    Nvec.resize(n);
    B_mat.zero();

    FESystemAssert4(((s.first == 2) && (s.second == 4*n)), FESystem::Numerics::MatrixSizeMismatch, 2, 4*n, s.first, s.second);
    
    Nvec.zero();
    this->finite_element->getShapeFunction(pt, Nvec);
    B_mat.setRowVals(0, 2*n,  3*n-1, Nvec); // gamma-xz:  thetay
    Nvec.scale(-1.0);
    B_mat.setRowVals(1, 3*n,  4*n-1, Nvec); // gamma-xy : thetaz
    
    Nvec.zero();
    this->finite_element->getShapeFunctionDerivativeForPhysicalCoordinates(derivatives_x, pt, Nvec);
    B_mat.setRowVals(0,   n,  2*n-1, Nvec); // gamma-xz:  w
    B_mat.setRowVals(1,   0,    n-1, Nvec); // gamma-xy:  v
}



void
FESystem::Structures::TimoshenkoBeam::getMaterialComplianceMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& bend_mat, FESystem::Numerics::MatrixBase<FESystemDouble>& shear_mat)
{
    const std::pair<FESystemUInt, FESystemUInt> s_b = bend_mat.getSize();
    const std::pair<FESystemUInt, FESystemUInt> s_s = shear_mat.getSize();
    
    FESystemAssert4(((s_b.first == 2) && (s_b.second== 2)), FESystem::Numerics::MatrixSizeMismatch, 2, 2, s_b.first, s_b.second);
    FESystemAssert4(((s_s.first == 2) && (s_s.second== 2)), FESystem::Numerics::MatrixSizeMismatch, 2, 2, s_s.first, s_s.second);
    
    shear_mat.setVal(0, 0, this->G_val);
    shear_mat.setVal(1, 1, this->G_val);

    bend_mat.setVal(0, 0, this->E_val * this->I_ch_val);
    bend_mat.setVal(1, 1, this->E_val * this->I_tr_val);
    
    shear_mat.scale(this->kappa*this->area_val);
}


