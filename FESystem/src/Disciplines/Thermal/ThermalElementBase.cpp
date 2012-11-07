//
//  ThermalElementBase.cpp
//  FESystem
//
//  Created by Manav Bhatia on 9/6/12.
//
//


// C++ includes
#include <algorithm>


// FESystem includes
#include "Disciplines/Thermal/ThermalElementBase.h"
#include "Base/FESystemExceptions.h"
#include "Mesh/ElemBase.h"
#include "FiniteElems/FiniteElementBase.h"
#include "Quadrature/QuadratureBase.h"
#include "Geom/Point.h"
#include "Functions/FunctionMappingBase.h"
#include "Numerics/DenseMatrix.h"
#include "Numerics/LocalVector.h"



FESystem::Thermal::ThermalElementBase::ThermalElementBase():
if_initialized(false),
geometric_elem(NULL),
quadrature(NULL),
finite_element(NULL),
k_val(0.0),
cp_val(0.0),
rho_val(0.0),
cross_section_val(0.0)
{
    
}


FESystem::Thermal::ThermalElementBase::~ThermalElementBase()
{
    
}


std::string
FESystem::Thermal::ThermalElementBase::getVariableName(FESystem::Thermal::ThermalVariable var)
{
    std::string name;
    
    switch (var)
    {
        case FESystem::Thermal::TEMP:
            name = "TEMP";
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
    }
    
    return name;
}


FESystem::Thermal::ThermalVariable
FESystem::Thermal::ThermalElementBase::getVariableEnum(const std::string &var)
{
    FESystem::Thermal::ThermalVariable name;
    
    if (var == "TEMP")
        name = FESystem::Thermal::TEMP;
    else
        FESystemAssert0(false, FESystem::Exception::InvalidValue);
    
    return name;
}


void
FESystem::Thermal::ThermalElementBase::clear()
{
    this->geometric_elem = NULL;
    this->quadrature = NULL;
    this->finite_element = NULL;
    this->k_val  = 0.0;
    this->cp_val = 0.0;
    this->rho_val = 0.0;
    this->cross_section_val = 0.0;
    
    this->if_initialized = false;
}



void
FESystem::Thermal::ThermalElementBase::initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe,
                                                  const FESystem::Quadrature::QuadratureBase& q_rule, FESystemDouble k, FESystemDouble cp, FESystemDouble rho, FESystemDouble cs)
{
    FESystemAssert0(!this->if_initialized, FESystem::Exception::InvalidFunctionCall);
    
    this->if_initialized = true;
    this->geometric_elem = &elem;
    this->quadrature = &q_rule;
    this->finite_element = &fe;
    
    this->k_val  = k;
    this->cp_val = cp;
    this->rho_val = rho;
    if (this->geometric_elem->getDimension() < 3)
        this->cross_section_val = cs;
    else
        this->cross_section_val = 1.0;  // set this to one for 3D cases since that calculates the volume through quadrature
}




void
FESystem::Thermal::ThermalElementBase::getFluxVector(const FESystem::Numerics::VectorBase<FESystemDouble>& pt, const FESystem::Numerics::VectorBase<FESystemDouble>& sol,
                                                     FESystem::Numerics::VectorBase<FESystemDouble>& vec)
{
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}




void
FESystem::Thermal::ThermalElementBase::calculateConsistentCapacitanceMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    const FESystemUInt n = this->geometric_elem->getNNodes();
    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    
    FESystemAssert4(((s.first == n) && (s.second== n)), FESystem::Numerics::MatrixSizeMismatch, n, n, s.first, s.second);
    
    FESystem::Numerics::DenseMatrix<FESystemDouble> B_mat, C_mat, tmp_mat1, tmp_mat2;
    C_mat.resize(1,1); B_mat.resize(1, n); tmp_mat1.resize(1, n), tmp_mat2.resize(n, n);
    C_mat.zero(); B_mat.zero(); tmp_mat1.zero(); tmp_mat2.zero();
    
    const std::vector<FESystem::Geometry::Point*>& q_pts = this->quadrature->getQuadraturePoints();
    const std::vector<FESystemDouble>& q_weight = this->quadrature->getQuadraturePointWeights();
    
    FESystemDouble jac=0.0;
    mat.zero();
    this->getMaterialCapacitanceMatrix(C_mat);
    
    for (FESystemUInt i=0; i<q_pts.size(); i++)
    {
        jac = this->finite_element->getJacobianValue(*(q_pts[i]));
        this->calculateCapacitanceOperatorMatrix(*(q_pts[i]), B_mat);
        
        C_mat.matrixRightMultiply(1.0, B_mat, tmp_mat1);
        B_mat.matrixTransposeRightMultiply(1.0, tmp_mat1, tmp_mat2);
        
        mat.add(q_weight[i]*jac, tmp_mat2);
    }
}



void
FESystem::Thermal::ThermalElementBase::calculateDiagonalCapaciatanceMatrix(FESystem::Numerics::VectorBase<FESystemDouble>& vec)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    const FESystemUInt n = this->geometric_elem->getNNodes();
    
    FESystemAssert2(vec.getSize() == n, FESystem::Exception::DimensionsDoNotMatch, n, vec.getSize());
    
    vec.zero();
    FESystemDouble val = this->geometric_elem->getElementSize(*(this->finite_element), *(this->quadrature)) * this->rho_val * this->cp_val * this->cross_section_val;
    
    val /= (1.0*n);
    vec.setAllVals(val);
}




void
FESystem::Thermal::ThermalElementBase::calculateConductanceMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);

    const FESystemUInt n = this->geometric_elem->getNNodes(), dim = this->geometric_elem->getDimension();
    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    
    FESystemAssert4(((s.first == n) && (s.second== n)), FESystem::Numerics::MatrixSizeMismatch, n, n, s.first, s.second);
    
    FESystem::Numerics::DenseMatrix<FESystemDouble> B_mat, C_mat, tmp_mat1, tmp_mat2;
    C_mat.resize(dim,dim); B_mat.resize(dim, n); tmp_mat1.resize(dim, n), tmp_mat2.resize(n, n);
    C_mat.zero(); B_mat.zero(); tmp_mat1.zero(); tmp_mat2.zero();
    
    const std::vector<FESystem::Geometry::Point*>& q_pts = this->quadrature->getQuadraturePoints();
    const std::vector<FESystemDouble>& q_weight = this->quadrature->getQuadraturePointWeights();
    
    FESystemDouble jac=0.0;
    mat.zero();
    this->getMaterialConductanceMatrix(C_mat);
    
    for (FESystemUInt i=0; i<q_pts.size(); i++)
    {
        jac = this->finite_element->getJacobianValue(*(q_pts[i]));
        this->calculateConductanceOperatorMatrix(*(q_pts[i]), B_mat);
        
        C_mat.matrixRightMultiply(1.0, B_mat, tmp_mat1);
        B_mat.matrixTransposeRightMultiply(1.0, tmp_mat1, tmp_mat2);
        
        mat.add(q_weight[i]*jac, tmp_mat2);
    }
}



void
FESystem::Thermal::ThermalElementBase::calculateCapacitanceOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);

    const FESystemUInt n = this->geometric_elem->getNNodes(), dim = this->geometric_elem->getDimension();
    const std::pair<FESystemUInt, FESystemUInt> s = B_mat.getSize();
    
    FESystemAssert4(((s.first == dim) && (s.second== n)), FESystem::Numerics::MatrixSizeMismatch, dim, n, s.first, s.second);
    
    FESystem::Numerics::LocalVector<FESystemDouble> Nvec;
    Nvec.resize(n); Nvec.zero();
    B_mat.zero();
    
    this->finite_element->getShapeFunction(pt, Nvec);
    
    for (FESystemUInt i_dim=0; i_dim < dim; i_dim++)
        B_mat.setRowVals(i_dim, 0, n-1, Nvec);
}




void
FESystem::Thermal::ThermalElementBase::calculateConductanceOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    const FESystemUInt n = this->geometric_elem->getNNodes(), dim = this->geometric_elem->getDimension();
    const std::pair<FESystemUInt, FESystemUInt> s = B_mat.getSize();
    
    FESystemAssert4(((s.first == dim) && (s.second== n)), FESystem::Numerics::MatrixSizeMismatch, dim, n, s.first, s.second);
    
    FESystem::Numerics::LocalVector<FESystemDouble> Nvec;
    Nvec.resize(n); Nvec.zero();
    B_mat.zero();
    
    std::vector<FESystemUInt> derivatives;
    derivatives.resize(dim);
    
    for (FESystemUInt i_dim=0; i_dim<dim; i_dim++)
    {
        std::fill(derivatives.begin(), derivatives.end(), 0);
        derivatives[i_dim-1] = 1; // first order derivative of the i_dim dimension
        
        this->finite_element->getShapeFunctionDerivativeForPhysicalCoordinates(derivatives, pt, Nvec);
        B_mat.setRowVals(i_dim, 0, n-1, Nvec);
    }
}




void
FESystem::Thermal::ThermalElementBase::getMaterialCapacitanceMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    const FESystemUInt dim = this->geometric_elem->getDimension();
    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    
    FESystemAssert4(((s.first == dim) && (s.second== dim)), FESystem::Numerics::MatrixSizeMismatch, dim, dim, s.first, s.second);

    mat.zero();
    for (FESystemUInt i_dim=0; i_dim < dim; i_dim++)
        mat.setVal(i_dim, i_dim, this->rho_val * this->cp_val * this->cross_section_val);
}



void
FESystem::Thermal::ThermalElementBase::getMaterialConductanceMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    const FESystemUInt dim = this->geometric_elem->getDimension();
    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    
    FESystemAssert4(((s.first == dim) && (s.second== dim)), FESystem::Numerics::MatrixSizeMismatch, dim, dim, s.first, s.second);
        
    mat.zero();
    for (FESystemUInt i_dim=0; i_dim < dim; i_dim++)
        mat.setVal(i_dim, i_dim, this->k_val * this->cross_section_val);
}


