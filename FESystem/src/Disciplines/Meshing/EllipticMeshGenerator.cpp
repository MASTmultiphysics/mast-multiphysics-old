//
//  EllipticMeshGenerator.cpp
//  FESystem
//
//  Created by Manav Bhatia on 12/23/12.
//
//


// FESystem includes
#include "Disciplines/Meshing/EllipticMeshGenerator.h"
#include "Base/FESystemExceptions.h"
#include "Numerics/DenseMatrix.h"
#include "Numerics/LocalVector.h"
#include "Mesh/ElemBase.h"
#include "FiniteElems/FiniteElementBase.h"
#include "Quadrature/QuadratureBase.h"
#include "Geom/Point.h"



FESystem::Meshing::EllipticMeshGenerator::EllipticMeshGenerator():
if_initialized(false),
geometric_elem(NULL),
quadrature(NULL),
finite_element(NULL),
solution(NULL),
P(0.0),
Q(0.0),
coupling_coeff(0.0)
{
    
}




FESystem::Meshing::EllipticMeshGenerator::~EllipticMeshGenerator()
{
    
}




void
FESystem::Meshing::EllipticMeshGenerator::clear()
{
    this->if_initialized = false;
    this->geometric_elem = NULL;
    this->quadrature = NULL;
    this->finite_element = NULL;
    this->solution = NULL;
    this->P = 0.0;
    this->Q = 0.0;
    this->coupling_coeff = 0.0;
}



void
FESystem::Meshing::EllipticMeshGenerator::initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, const FESystem::Quadrature::QuadratureBase& q_rule,
                                                     const FESystem::Numerics::VectorBase<FESystemDouble>& sol, const FESystemDouble p_val, const FESystemDouble q_val, const FESystemDouble coeff)
{
    FESystemAssert0(!this->if_initialized, FESystem::Exception::InvalidState);
    
    this->geometric_elem = &elem;
    this->quadrature = &q_rule;
    this->finite_element = &fe;
    this->solution = &sol;
    this->P = p_val;
    this->Q = q_val;
    this->coupling_coeff = coeff;
    
    this->if_initialized = true;
}




void
FESystem::Meshing::EllipticMeshGenerator::calculateResidual(FESystem::Numerics::VectorBase<FESystemDouble>& res)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    FESystemUInt dim = this->geometric_elem->getDimension(), n=this->geometric_elem->getNNodes(), n1 = dim*n;
    
    FESystemAssert2(res.getSize() == n1, FESystem::Exception::DimensionsDoNotMatch, res.getSize(), n1);
    
    FESystem::Numerics::DenseMatrix<FESystemDouble> Bmat, Bmat_dxi, Bmat_deta;
    FESystem::Numerics::LocalVector<FESystemDouble> Nvec, dNdxi, dNdeta, xy, dxy_xi, dxy_eta, tmp_vec1, tmp_vec2;
    Nvec.resize(n); dNdxi.resize(n); dNdeta.resize(n); xy.resize(dim); dxy_xi.resize(dim); dxy_eta.resize(dim);
    Bmat.resize(dim, n1); Bmat_dxi.resize(dim, n1); Bmat_deta.resize(dim, n1); tmp_vec1.resize(dim); tmp_vec2.resize(n1);
    
    const std::vector<FESystem::Geometry::Point*>& q_pts = this->quadrature->getQuadraturePoints();
    const std::vector<FESystemDouble>& q_weight = this->quadrature->getQuadraturePointWeights();
    
    FESystemDouble jac, g11, g22, g12, gval;
    std::vector<FESystemUInt> dxi_index(dim);
    
    res.zero();
    
    for (FESystemUInt i=0; i<q_pts.size(); i++)
    {
        jac = this->finite_element->getJacobianValue(*(q_pts[i]));

        this->finite_element->getShapeFunction(*(q_pts[i]), Nvec);
        std::fill(dxi_index.begin(), dxi_index.end(), 0); dxi_index[0] = 1;
        this->finite_element->getShapeFunctionDerivativeForPhysicalCoordinates(dxi_index, *(q_pts[i]), dNdxi);
        std::fill(dxi_index.begin(), dxi_index.end(), 0); dxi_index[1] = 1;
        this->finite_element->getShapeFunctionDerivativeForPhysicalCoordinates(dxi_index, *(q_pts[i]), dNdeta);
        
        Bmat.setRowVals(0, 0, n-1, Nvec); Bmat.setRowVals(1, n, n1-1, Nvec);
        Bmat_dxi.setRowVals(0, 0, n-1, dNdxi); Bmat_dxi.setRowVals(1, n, n1-1, dNdxi);
        Bmat_deta.setRowVals(0, 0, n-1, dNdeta); Bmat_deta.setRowVals(1, n, n1-1, dNdeta);
        
        Bmat.rightVectorMultiply(*(this->solution), xy);
        Bmat_dxi.rightVectorMultiply(*(this->solution), dxy_xi);
        Bmat_deta.rightVectorMultiply(*(this->solution), dxy_eta);
        
        g11 = pow(dxy_xi.getL2Norm(),2);  g22 = pow(dxy_eta.getL2Norm(),2);  g12 = dxy_xi.getVal(0)*dxy_eta.getVal(0) + dxy_xi.getVal(1)*dxy_eta.getVal(1);
        gval = g11*g22-g12*g12;
        
        Bmat_dxi.rightVectorMultiply(*(this->solution), tmp_vec1);  // dX/dxi
        Bmat_dxi.leftVectorMultiply(tmp_vec1, tmp_vec2); // dB/dxi^T  dX/dxi
        res.add(q_weight[i]*jac*g22, tmp_vec2); // dB/dxi^T g22  dX/dxi
        
        Bmat_deta.leftVectorMultiply(tmp_vec1, tmp_vec2); // dB/deta^T  dX/dxi
        res.add(-q_weight[i]*jac*g12*coupling_coeff, tmp_vec2); // -dB/deta^T g12  dX/dxi
        
        Bmat.leftVectorMultiply(tmp_vec1, tmp_vec2); // B^T  dX/dxi
        res.add(-q_weight[i]*jac*gval*this->P, tmp_vec2); // -B^T g dX/dxi

        
        Bmat_deta.rightVectorMultiply(*(this->solution), tmp_vec1);  // dX/deta
        Bmat_dxi.leftVectorMultiply(tmp_vec1, tmp_vec2); // dB/dxi^T  dX/deta
        res.add(-q_weight[i]*jac*g12*coupling_coeff, tmp_vec2); // -dB/dxi^T g12  dX/deta
        
        Bmat_deta.leftVectorMultiply(tmp_vec1, tmp_vec2); // dB/deta^T  dX/deta
        res.add(q_weight[i]*jac*g11, tmp_vec2); // dB/deta^T g11  dX/deta
        
        Bmat.leftVectorMultiply(tmp_vec1, tmp_vec2); // B^T  dX/deta
        res.add(-q_weight[i]*jac*gval*this->Q, tmp_vec2); // -B^T g  dX/deta
    }
}




void
FESystem::Meshing::EllipticMeshGenerator::calculateTangentMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    FESystemUInt dim = this->geometric_elem->getDimension(), n=this->geometric_elem->getNNodes(), n1 = dim*n;
    
    std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    
    FESystemAssert4((s.first == n1) && (s.second = n1), FESystem::Numerics::MatrixSizeMismatch, s.first, s.second, n1, n1);
    
    FESystem::Numerics::DenseMatrix<FESystemDouble> Bmat, Bmat_dxi, Bmat_deta, tmp_mat;
    FESystem::Numerics::LocalVector<FESystemDouble> Nvec, dNdxi, dNdeta, xy, dxy_xi, dxy_eta;
    Nvec.resize(n); dNdxi.resize(n); dNdeta.resize(n); xy.resize(dim); dxy_xi.resize(dim); dxy_eta.resize(dim);
    Bmat.resize(dim, n1); Bmat_dxi.resize(dim, n1); Bmat_deta.resize(dim, n1); tmp_mat.resize(n1, n1);
    
    const std::vector<FESystem::Geometry::Point*>& q_pts = this->quadrature->getQuadraturePoints();
    const std::vector<FESystemDouble>& q_weight = this->quadrature->getQuadraturePointWeights();
    
    FESystemDouble jac, g11, g22, g12, gval;
    std::vector<FESystemUInt> dxi_index(dim);
    
    mat.zero();
    
    for (FESystemUInt i=0; i<q_pts.size(); i++)
    {
        jac = this->finite_element->getJacobianValue(*(q_pts[i]));
        
        this->finite_element->getShapeFunction(*(q_pts[i]), Nvec);
        std::fill(dxi_index.begin(), dxi_index.end(), 0); dxi_index[0] = 1;
        this->finite_element->getShapeFunctionDerivativeForPhysicalCoordinates(dxi_index, *(q_pts[i]), dNdxi);
        std::fill(dxi_index.begin(), dxi_index.end(), 0); dxi_index[1] = 1;
        this->finite_element->getShapeFunctionDerivativeForPhysicalCoordinates(dxi_index, *(q_pts[i]), dNdeta);
        
        Bmat.setRowVals(0, 0, n-1, Nvec); Bmat.setRowVals(1, n, n1-1, Nvec);
        Bmat_dxi.setRowVals(0, 0, n-1, dNdxi); Bmat_dxi.setRowVals(1, n, n1-1, dNdxi);
        Bmat_deta.setRowVals(0, 0, n-1, dNdeta); Bmat_deta.setRowVals(1, n, n1-1, dNdeta);
        
        Bmat.rightVectorMultiply(*(this->solution), xy);
        Bmat_dxi.rightVectorMultiply(*(this->solution), dxy_xi);
        Bmat_deta.rightVectorMultiply(*(this->solution), dxy_eta);
        
        g11 = pow(dxy_xi.getL2Norm(),2);  g22 = pow(dxy_eta.getL2Norm(),2);  g12 = dxy_xi.getVal(0)*dxy_eta.getVal(0) + dxy_xi.getVal(1)*dxy_eta.getVal(1);
        gval = g11*g22-g12*g12;
        
        Bmat_dxi.matrixTransposeRightMultiply(1.0, Bmat_dxi, tmp_mat); // dB/dxi^T  dB/dxi
        mat.add(q_weight[i]*jac*g22, tmp_mat); // dB/dxi^T g22  dB/dxi
        
        Bmat_deta.matrixTransposeRightMultiply(1.0, Bmat_dxi, tmp_mat); // dB/deta^T  dB/dxi
        mat.add(-q_weight[i]*jac*g12*coupling_coeff, tmp_mat); // -dB/deta^T g12  dX/dxi
        
        Bmat.matrixTransposeRightMultiply(1.0, Bmat_dxi, tmp_mat); // B^T  dX/dxi
        mat.add(-q_weight[i]*jac*gval*this->P, tmp_mat); // -B^T g dB/dxi
        
        
        Bmat_dxi.matrixTransposeRightMultiply(1.0, Bmat_deta, tmp_mat); // dB/dxi^T  dB/deta
        mat.add(-q_weight[i]*jac*g12*coupling_coeff, tmp_mat); // -dB/dxi^T g12  dB/deta
        
        Bmat_deta.matrixTransposeRightMultiply(1.0, Bmat_deta, tmp_mat); // dB/deta^T  dB/deta
        mat.add(q_weight[i]*jac*g11, tmp_mat); // dB/deta^T g11  dB/deta
        
        Bmat.matrixTransposeRightMultiply(1.0, Bmat_deta, tmp_mat); // B^T  dB/deta
        mat.add(-q_weight[i]*jac*gval*this->Q, tmp_mat); // -B^T g  dB/deta
    }
}


