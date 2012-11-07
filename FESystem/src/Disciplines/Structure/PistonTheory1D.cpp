//
//  PistonTheory.cpp
//  FESystem
//
//  Created by Manav Bhatia on 9/13/12.
//
//

// FESystem includes
#include "Disciplines/Structure/PistonTheory1D.h"
#include "Base/FESystemExceptions.h"
#include "Mesh/ElemBase.h"
#include "Numerics/DenseMatrix.h"
#include "Numerics/LocalVector.h"
#include "FiniteElems/FiniteElementBase.h"
#include "Quadrature/QuadratureBase.h"
#include "Geom/Point.h"


FESystem::Structures::PistonTheory1D::PistonTheory1D():
FESystem::Structures::Structural1DElementBase(),
theory_order(0),
mach(0.0),
a_inf(0.0),
gamma(0.0),
u_inf(0.0)
{
    
}


FESystem::Structures::PistonTheory1D::~PistonTheory1D()
{
    
}


void
FESystem::Structures::PistonTheory1D::clear()
{
    FESystem::Structures::Structural1DElementBase::clear();
    this->theory_order = 0;
    this->mach = 0.0;
    this->a_inf = 0.0;
    this->gamma = 0.0;
    this->u_inf = 0.0;
}



FESystemUInt
FESystem::Structures::PistonTheory1D::getNElemDofs() const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    return this->geometric_elem->getNNodes();
}


void
FESystem::Structures::PistonTheory1D::getActiveElementMatrixIndices(std::vector<FESystemUInt>& vec)
{
    FESystemUInt n = this->geometric_elem->getNNodes();
    vec.resize(n);
    
    for (FESystemUInt i=0; i<n; i++) vec[i] = 2*n+i; // w-displacement
}


void
FESystem::Structures::PistonTheory1D::initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, const FESystem::Quadrature::QuadratureBase& q_rule,
                                                 FESystemUInt order, FESystemDouble m, FESystemDouble a, FESystemDouble u, FESystemDouble g)
{
    FESystemAssert0((order >= 1) && (order <= 3), FESystem::Exception::InvalidFunctionCall);
    
    FESystem::Structures::Structural1DElementBase::initialize(elem, fe, q_rule, 0.0, 0.0, 0.0, 0.0);
    this->theory_order = order;
    this->mach = m;
    this->a_inf = a;
    this->gamma = g;
    this->u_inf = u;
}


void
FESystem::Structures::PistonTheory1D::getStressTensor(const FESystem::Geometry::Point& pt, const FESystem::Numerics::VectorBase<FESystemDouble>& sol,
                                                    FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}



void
FESystem::Structures::PistonTheory1D::calculateForceVector(const FESystem::Numerics::VectorBase<FESystemDouble>& sol, const FESystem::Numerics::VectorBase<FESystemDouble>& vel, FESystem::Numerics::VectorBase<FESystemDouble>& force)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    const FESystemUInt n = this->geometric_elem->getNNodes();
    
    FESystemAssert2(sol.getSize() == n, FESystem::Exception::DimensionsDoNotMatch, n, sol.getSize());
    FESystemAssert2(vel.getSize() == n, FESystem::Exception::DimensionsDoNotMatch, n, vel.getSize());
    FESystemAssert2(force.getSize() == n, FESystem::Exception::DimensionsDoNotMatch, n, force.getSize());
    
    FESystem::Numerics::LocalVector<FESystemDouble> tmp_vec1, Nvec;
    std::vector<FESystemUInt> derivatives(1);
    tmp_vec1.resize(n); Nvec.resize(n);
    tmp_vec1.zero(); Nvec.zero();
    derivatives[0] = 1;
    
    const std::vector<FESystem::Geometry::Point*>& q_pts = this->quadrature->getQuadraturePoints();
    const std::vector<FESystemDouble>& q_weight = this->quadrature->getQuadraturePointWeights();
    
    FESystemDouble jac=0.0, wdot, cp;
    force.zero();
    
    for (FESystemUInt i=0; i<q_pts.size(); i++)
    {
        jac = this->finite_element->getJacobianValue(*(q_pts[i]));
        
        // calculate the coefficient of the force vector
        this->finite_element->getShapeFunctionDerivativeForPhysicalCoordinates(derivatives, *(q_pts[i]), Nvec);
        wdot = this->u_inf*Nvec.dotProduct(sol);
        
        this->finite_element->getShapeFunction(*(q_pts[i]), Nvec);
        wdot += Nvec.dotProduct(vel);
        wdot /= this->a_inf;
        
        // calculate the cp value
        cp = 0.0;
        switch (this->theory_order)
        {
            case 3:
                cp += (gamma+1.0)/12.0*pow(wdot,3);
            case 2:
                cp += (gamma+1.0)/4.0*pow(wdot,2);
            case 1:
                cp += wdot;
                cp *= -2.0/pow(this->mach,2); // multiply by -1 since +ve Cp implies -ve force
                break;

            default:
                FESystemAssert0(false, FESystem::Exception::InvalidValue);
                break;
        }
        force.add(q_weight[i]*jac*cp, Nvec);
    }
}



void
FESystem::Structures::PistonTheory1D::calculateTangentMatrix(const FESystem::Numerics::VectorBase<FESystemDouble>& sol, const FESystem::Numerics::VectorBase<FESystemDouble>& vel,
                                                             FESystem::Numerics::MatrixBase<FESystemDouble>& dfdw, FESystem::Numerics::MatrixBase<FESystemDouble>& dfdwdot)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    const FESystemUInt n = this->geometric_elem->getNNodes();
    const std::pair<FESystemUInt, FESystemUInt> s1 = dfdw.getSize(), s2 = dfdwdot.getSize();
    
    FESystemAssert2(sol.getSize() == n, FESystem::Exception::DimensionsDoNotMatch, n, sol.getSize());
    FESystemAssert2(vel.getSize() == n, FESystem::Exception::DimensionsDoNotMatch, n, vel.getSize());
    FESystemAssert4((s1.first == n) && (s1.second == n), FESystem::Numerics::MatrixSizeMismatch, n, n, s1.first, s1.second);
    FESystemAssert4((s2.first == n) && (s2.second == n), FESystem::Numerics::MatrixSizeMismatch, n, n, s2.first, s2.second);
    
    
    FESystem::Numerics::DenseMatrix<FESystemDouble> Bmat, Bmatdx, tmp_mat1;
    FESystem::Numerics::LocalVector<FESystemDouble> tmp_vec1, Nvec;
    std::vector<FESystemUInt> derivatives(1);
    Bmat.resize(1, n); Bmatdx.resize(1, n); tmp_mat1.resize(n, n);
    tmp_vec1.resize(n); Nvec.resize(n);
    tmp_vec1.zero(); Nvec.zero();
    derivatives[0] = 1;
    
    const std::vector<FESystem::Geometry::Point*>& q_pts = this->quadrature->getQuadraturePoints();
    const std::vector<FESystemDouble>& q_weight = this->quadrature->getQuadraturePointWeights();
    
    FESystemDouble jac=0.0, wdot, cp;
    dfdw.zero(); dfdwdot.zero();
    
    for (FESystemUInt i=0; i<q_pts.size(); i++)
    {
        jac = this->finite_element->getJacobianValue(*(q_pts[i]));
        
        // calculate the coefficient of the force vector
        this->finite_element->getShapeFunctionDerivativeForPhysicalCoordinates(derivatives, *(q_pts[i]), Nvec);
        Bmatdx.setRowVals(0, 0, n-1, Nvec);
        wdot = this->u_inf*Nvec.dotProduct(sol);
        
        this->finite_element->getShapeFunction(*(q_pts[i]), Nvec);
        Bmat.setRowVals(0, 0, n-1, Nvec);
        wdot += Nvec.dotProduct(vel);
        wdot /= this->a_inf;
        
        // calculate the cp value
        cp = 0.0;
        switch (this->theory_order)
        {
            case 3:
                cp += (gamma+1.0)/4.0*pow(wdot,2);
            case 2:
                cp += (gamma+1.0)/2.0*wdot;
            case 1:
                cp += 1.0;
                cp *= -2.0/pow(this->mach,2); // multiply by -1 since +ve Cp implies -ve force
                break;
                
            default:
                FESystemAssert0(false, FESystem::Exception::InvalidValue);
                break;
        }
        // calculate the mass term
        tmp_mat1.zero();
        Bmat.matrixTransposeRightMultiply(q_weight[i]*jac*cp/this->a_inf, Bmat, tmp_mat1);
        dfdwdot.add(1.0, tmp_mat1);
        
        // calculate the stiffness term
        tmp_mat1.zero();
        Bmat.matrixTransposeRightMultiply(q_weight[i]*jac*cp*this->u_inf/this->a_inf, Bmatdx, tmp_mat1);
        dfdw.add(1.0, tmp_mat1);
    }
}



