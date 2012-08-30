//
//  EulerBernoulliBeam.cpp
//  FESystem
//
//  Created by Manav Bhatia on 8/13/12.
//
//


// FESystem includes
#include "Disciplines/Structure/EulerBernoulliBeam.h"
#include "Base/FESystemExceptions.h"
#include "Mesh/ElemBase.h"
#include "Numerics/DenseMatrix.h"
#include "Numerics/LocalVector.h"
#include "FiniteElems/FiniteElementBase.h"
#include "Quadrature/QuadratureBase.h"
#include "Geom/Point.h"


FESystem::Structures::EulerBernoulliBeam::EulerBernoulliBeam():
FESystem::Structures::LinearBeamElementBase()
{
    
}


FESystem::Structures::EulerBernoulliBeam::~EulerBernoulliBeam()
{
    
}



void
FESystem::Structures::EulerBernoulliBeam::clear()
{
    FESystem::Structures::LinearBeamElementBase::clear();
}





void
FESystem::Structures::EulerBernoulliBeam::getStressTensor(const FESystem::Numerics::VectorBase<FESystemDouble>& pt, const FESystem::Numerics::VectorBase<FESystemDouble>& sol,
                                                      FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}




void
FESystem::Structures::EulerBernoulliBeam::initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, const FESystem::Quadrature::QuadratureBase& q_rule,
                                                 FESystemDouble E, FESystemDouble nu, FESystemDouble rho, FESystemDouble I_tr, FESystemDouble I_ch, FESystemDouble A, FESystemBoolean if_lateral)
{
    FESystemAssert0(elem.getElementType() == FESystem::Mesh::EDGE2, FESystem::Exception::InvalidValue);

    FESystem::Structures::LinearBeamElementBase::initialize(elem, fe, q_rule, E, nu, rho, I_tr, I_ch, A, if_lateral);
}



void
FESystem::Structures::EulerBernoulliBeam::calculateStiffnessMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    const FESystemUInt n = this->geometric_elem->getNNodes();
    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    FESystemDouble length=this->geometric_elem->getElementSize(*(this->finite_element), *(this->quadrature));
    
    static FESystem::Numerics::DenseMatrix<FESystemDouble> B_mat, B_mat_shear, C_mat_bend, tmp_mat1, tmp_mat2;
    
    if (this->if_include_lateral_stiffness)
    {
        FESystemAssert4(((s.first == 4*n) && (s.second == 4*n)), FESystem::Numerics::MatrixSizeMismatch, 4*n, 4*n, s.first, s.second);
        C_mat_bend.resize(1,1); B_mat.resize(1, 4*n); B_mat_shear.resize(2, 4*n); tmp_mat1.resize(1, 4*n), tmp_mat2.resize(4*n, 4*n);
    }
    else
    {
        FESystemAssert4(((s.first == 2*n) && (s.second == 2*n)), FESystem::Numerics::MatrixSizeMismatch, 2*n, 2*n, s.first, s.second);
        C_mat_bend.resize(1,1); B_mat.resize(1, 2*n); B_mat_shear.resize(1, 2*n); tmp_mat1.resize(1, 2*n), tmp_mat2.resize(2*n, 2*n);
    }
    
    C_mat_bend.zero(); B_mat.zero(); tmp_mat1.zero(); tmp_mat2.zero();
    
    
    const std::vector<FESystem::Geometry::Point*>& q_pts_bend = this->quadrature->getQuadraturePoints();
    const std::vector<FESystemDouble>& q_weight_bend = this->quadrature->getQuadraturePointWeights();
    
    FESystemDouble jac=0.0;
    mat.zero();
    this->getMaterialComplianceMatrix(C_mat_bend);
    
    // bending contribution
    for (FESystemUInt i=0; i<q_pts_bend.size(); i++)
    {
        jac = this->finite_element->getJacobianValue(*(q_pts_bend[i]));
        this->calculateBendingOperatorMatrix(*(q_pts_bend[i]), length, B_mat);
        
        C_mat_bend.matrixRightMultiply(1.0, B_mat, tmp_mat1);
        B_mat.matrixTransposeRightMultiply(1.0, tmp_mat1, tmp_mat2);
        
        mat.add(q_weight_bend[i]*jac, tmp_mat2);
    }
}




void
FESystem::Structures::EulerBernoulliBeam::calculateBendingOperatorMatrix(const FESystem::Geometry::Point& pt, FESystemDouble length, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat)
{
    const FESystemUInt n = this->geometric_elem->getNNodes();
    const std::pair<FESystemUInt, FESystemUInt> s = B_mat.getSize();
    
    static FESystemDouble N1, N2, N3, N4;
    FESystemDouble xi = pt.getVal(0);
    
    // shape function values
    // N1 = (length/8.0) * (4.0/length -  6.0/length*xi + 0.0 +  2.0/length*pow(xi,3));
    // N2 = (length/8.0) * (4.0/length +  6.0/length*xi + 0.0 -  2.0/length*pow(xi,3));
    // N3 = (length/8.0) * ( 1.0 - xi - pow(xi,2) + pow(xi,3));  // needs a -1.0 factor for theta_y
    // N4 = (length/8.0) * (-1.0 - xi + pow(xi,2) + pow(xi,3));  // needs a -1.0 factor for theta_y
    
    // shape function derivative
    // N1 = (1.0/4.0) * (0.0 -  6.0/length + 0.0     +  6.0/length*pow(xi,2));
    // N2 = (1.0/4.0) * (0.0 +  6.0/length + 0.0     -  6.0/length*pow(xi,2));
    // N3 = (1.0/4.0) * (0.0 -  1.0        - 2.0*xi  +         3.0*pow(xi,2));  // needs a -1.0 factor for theta_y
    // N4 = (1.0/4.0) * (0.0 -  1.0        + 2.0*xi  +         3.0*pow(xi,2));  // needs a -1.0 factor for theta_y

    // second order shape function derivative
    N1 = (0.5/length) * (  0.0     +  12.0/length*xi);
    N2 = (0.5/length) * (  0.0     -  12.0/length*xi);
    N3 = (0.5/length) * (- 2.0     +          6.0*xi);  // needs a -1.0 factor for theta_y
    N4 = (0.5/length) * (  2.0     +          6.0*xi);  // needs a -1.0 factor for theta_y
    
    
    B_mat.zero();
    if (this->if_include_lateral_stiffness)
    {
        FESystemAssert4(((s.first == 1) && (s.second == 4*n)), FESystem::Numerics::MatrixSizeMismatch, 1, 4*n, s.first, s.second);
        B_mat.setVal(0, 0, N1*(this->I_ch_val/this->I_tr_val)); B_mat.setVal(0, 1, N2*(this->I_ch_val/this->I_tr_val)); // v-disp
        B_mat.setVal(0, 2, N1); B_mat.setVal(0, 3, N2); // w-disp

        B_mat.setVal(0, 4,-N3); B_mat.setVal(0, 5,-N4); // theta-y
        B_mat.setVal(0, 6, N3*(this->I_ch_val/this->I_tr_val)); B_mat.setVal(0, 7, N4*(this->I_ch_val/this->I_tr_val)); // theta-z
    }
    else
    {
        FESystemAssert4(((s.first == 1) && (s.second == 2*n)), FESystem::Numerics::MatrixSizeMismatch, 1, 2*n, s.first, s.second);
        B_mat.setVal(0, 0, N1); B_mat.setVal(0, 1, N2); // w-disp
        B_mat.setVal(0, 2,-N3); B_mat.setVal(0, 3,-N4); // theta-y
    }

}



void
FESystem::Structures::EulerBernoulliBeam::getMaterialComplianceMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& bend_mat)
{
    const std::pair<FESystemUInt, FESystemUInt> s_b = bend_mat.getSize();
    
    FESystemAssert4(((s_b.first == 1) && (s_b.second== 1)), FESystem::Numerics::MatrixSizeMismatch, 1, 1, s_b.first, s_b.second);
    bend_mat.setVal(0, 0, this->E_val * this->I_tr_val);
}



