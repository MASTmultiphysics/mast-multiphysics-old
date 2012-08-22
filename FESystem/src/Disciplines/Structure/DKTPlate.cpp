//
//  DKTPlate.cpp
//  FESystem
//
//  Created by Manav Bhatia on 8/13/12.
//
//



// FESystem includes
#include "Disciplines/Structure/DKTPlate.h"
#include "Base/FESystemExceptions.h"
#include "Mesh/ElemBase.h"
#include "Numerics/DenseMatrix.h"
#include "Numerics/LocalVector.h"
#include "FiniteElems/FiniteElementBase.h"
#include "Quadrature/QuadratureBase.h"
#include "Geom/Point.h"
#include "Mesh/Node.h"
#include "Mesh/Tri6.h"




namespace FESystem
{
    namespace FEBatoz
    {
        /*!
         *   returns the length of the side defined by vector from node i to j
         */
        double getEdgeLength(const FESystem::Mesh::ElemBase& elem, const FESystemUInt i, const FESystemUInt j)
        {
            const FESystem::Mesh::Node& point_i = elem.getNode(i);
            const FESystem::Mesh::Node& point_j = elem.getNode(j);
            
            FESystem::Numerics::LocalVector<FESystemDouble> vec_ji;
            vec_ji.resize(3);
            vec_ji.copyVector(point_j);
            vec_ji.add(-1.0, point_i);
            
            return vec_ji.getL2Norm();
        }
        
        /*!
         *   returns the cos of normal to the side defined by vector from node i to j
         */
        void getEdgeNormalSineCosine(const FESystem::Mesh::ElemBase& elem, const FESystemUInt i, const FESystemUInt j, FESystemDouble& sine, FESystemDouble& cosine)
        {
            static FESystem::Numerics::LocalVector<FESystemDouble> vec0, vec1, vec2, vec3;
            vec0.resize(3); vec1.resize(3); vec2.resize(3); vec3.resize(3);
            
            // calculate the normal to the element
            elem.getNodeLocationInLocalPhysicalCoordinate(0, vec0);
            elem.getNodeLocationInLocalPhysicalCoordinate(1, vec1);
            elem.getNodeLocationInLocalPhysicalCoordinate(2, vec2);
            
            vec1.add(-1.0, vec0);
            vec2.add(-1.0, vec0);
            vec1.crossProduct(vec2, vec3);
            // this is the unit vector normal to the plane of the triangle
            vec1.copyVector(vec3);
            vec1.scaleToUnitLength();
            
            // cross product of the length vector with the surface normal will
            // give the needed vector
            elem.getNodeLocationInLocalPhysicalCoordinate(i, vec3);
            elem.getNodeLocationInLocalPhysicalCoordinate(j, vec2);
            
            vec2.add(-1.0, vec3);
            vec2.crossProduct(vec1, vec3);
            vec3.scaleToUnitLength();        // this is the unit vector needed
            
            // cos of angle between this and the x-axis is simply the
            // 0th component of this vector
            cosine = vec3.getVal(0);
            sine = vec3.getVal(1);
        }
        
        
        void calculateDKTShapeFunction(FESystem::Mesh::ElemBase& elem,
                                       FESystem::Numerics::LocalVector<FESystemDouble>& tri6_shape_funcs,
                                       FESystem::Numerics::LocalVector<FESystemDouble>& betax,
                                       FESystem::Numerics::LocalVector<FESystemDouble>& betay)
        {
            // -- keep in mind that the index numbers for the elems start at 0.
            // -- also, the mid side node numbers in the Batoz's paper are different from
            //   the ones used in this library. Hence, use the following association
            //                   BATOZ TRI6 node #              FESystem TRI6 node #
            //                          1                               0
            //                          2                               1
            //                          3                               2
            //                          4 (on edge 2,3)                 4 (on edge 1,2)
            //                          5 (on edge 3,1)                 5 (on edge 2,0)
            //                          6 (on edge 1,2)                 3 (on edge 0,1)
            // -- all shape functions for this element are in a variable major format. Hence,
            //  they follow Hx for w1, w2, w3, thetax1, thetax2, thetax3, thetay1, thetay2, thetay3.
            // And then the same this is followed for Hy (from dofs 9-17)
            
            // local static variables for shape functions
            double N1, N2, N3, N4, N5, N6;
            
            // local static variables for edge lengths and sine/cosines
            double l12, l23, l31, cos4, cos5, cos6, sin4, sin5, sin6;
            
            N1 = tri6_shape_funcs.getVal(0);
            N2 = tri6_shape_funcs.getVal(1);
            N3 = tri6_shape_funcs.getVal(2);
            N4 = tri6_shape_funcs.getVal(4);
            N5 = tri6_shape_funcs.getVal(5);
            N6 = tri6_shape_funcs.getVal(3);
            
            FEBatoz::getEdgeNormalSineCosine(elem, 1, 2, sin4, cos4);
            FEBatoz::getEdgeNormalSineCosine(elem, 2, 0, sin5, cos5);
            FEBatoz::getEdgeNormalSineCosine(elem, 0, 1, sin6, cos6);
            
            l12 = FEBatoz::getEdgeLength(elem, 0, 1);
            l23 = FEBatoz::getEdgeLength(elem, 1, 2);
            l31 = FEBatoz::getEdgeLength(elem, 2, 0);
            
            betax.setVal( 0, 1.5 * (N5 * sin5 / l31 - N6 * sin6 / l12)); // Hx, w1
            betax.setVal( 1, 1.5 * (N6 * sin6 / l12 - N4 * sin4 / l23)); // Hx, w2
            betax.setVal( 2, 1.5 * (N4 * sin4 / l23 - N5 * sin5 / l31)); // Hx, w3
            betax.setVal( 3, (-.75) * (N5 * sin5 * cos5 + N6 * sin6 * cos6)); // Hx, thetax1
            betax.setVal( 4, (-.75) * (N4 * sin4 * cos4 + N6 * sin6 * cos6)); // Hx, thetax2
            betax.setVal( 5, (-.75) * (N5 * sin5 * cos5 + N4 * sin4 * cos4)); // Hx, thetax3
            betax.setVal( 6, (N1 + N5 * (0.5 * cos5 * cos5 - 0.25 * sin5 * sin5) + N6 * (0.5 * cos6 * cos6 - 0.25 * sin6 * sin6))); // Hx, thetay1
            betax.setVal( 7, (N2 + N4 * (0.5 * cos4 * cos4 - 0.25 * sin4 * sin4) + N6 * (0.5 * cos6 * cos6 - 0.25 * sin6 * sin6))); // Hx, thetay2
            betax.setVal( 8, (N3 + N5 * (0.5 * cos5 * cos5 - 0.25 * sin5 * sin5) + N4 * (0.5 * cos4 * cos4 - 0.25 * sin4 * sin4))); // Hx, thetay3
            
            betay.setVal( 0, 1.5 * (-N5 * cos5 / l31 + N6 * cos6 / l12)); // Hy, w1
            betay.setVal( 1, 1.5 * (-N6 * cos6 / l12 + N4 * cos4 / l23)); // Hy, w2
            betay.setVal( 2, 1.5 * (-N4 * cos4 / l23 + N5 * cos5 / l31)); // Hy, w3
            betay.setVal( 3, (-N1 + N5 * (0.25 * cos5 * cos5 - 0.5 * sin5 * sin5) + N6 * (0.25 * cos6 * cos6 - 0.5 * sin6 * sin6))); // Hy, thetax1
            betay.setVal( 4, (-N2 + N4 * (0.25 * cos4 * cos4 - 0.5 * sin4 * sin4) + N6 * (0.25 * cos6 * cos6 - 0.5 * sin6 * sin6))); // Hy, thetax2
            betay.setVal( 5, (-N3 + N5 * (0.25 * cos5 * cos5 - 0.5 * sin5 * sin5) + N4 * (0.25 * cos4 * cos4 - 0.5 * sin4 * sin4))); // Hy, thetax3
            betay.setVal( 6, (-1.0) * betax.getVal(3)); // Hy, thetay1
            betay.setVal( 7, (-1.0) * betax.getVal(4)); // Hy, thetay2
            betay.setVal( 8, (-1.0) * betax.getVal(5)); // Hy, thetay3
        }
    }
}



FESystem::Structures::DKTPlate::DKTPlate():
FESystem::Structures::Structural2DElementBase(),
finite_element_tri6(NULL),
quadrature_tri6(NULL)
{
    this->origin = new FESystem::Geometry::Point(2);
    this->tri6_elem = new FESystem::Mesh::Tri6;
    this->node0 = new FESystem::Mesh::Node(this->origin->getCoordinateSystem());
    this->node1 = new FESystem::Mesh::Node(this->origin->getCoordinateSystem());
    this->node2 = new FESystem::Mesh::Node(this->origin->getCoordinateSystem());
    this->node3 = new FESystem::Mesh::Node(this->origin->getCoordinateSystem());
    this->node4 = new FESystem::Mesh::Node(this->origin->getCoordinateSystem());
    this->node5 = new FESystem::Mesh::Node(this->origin->getCoordinateSystem());
    
}


FESystem::Structures::DKTPlate::~DKTPlate()
{
    // delete the element
    delete this->tri6_elem;
    delete this->node0;
    delete this->node1;
    delete this->node2;
    delete this->node3;
    delete this->node4;
    delete this->node5;
    delete this->origin;
}




void
FESystem::Structures::DKTPlate::clear()
{
    FESystem::Structures::Structural2DElementBase::clear();
    this->finite_element_tri6 = NULL;
    this->quadrature_tri6 = NULL;
}




void
FESystem::Structures::DKTPlate::transformMatrixToGlobalCoordinate(const std::vector<FESystem::Structures::StructuralVariable>& vars,
                                                                  const FESystem::Numerics::MatrixBase<FESystemDouble>& elem_cs_mat,
                                                                  FESystem::Numerics::MatrixBase<FESystemDouble>& global_cs_mat)
{
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}


void
FESystem::Structures::DKTPlate::getStressTensor(const FESystem::Numerics::VectorBase<FESystemDouble>& pt, const FESystem::Numerics::VectorBase<FESystemDouble>& sol,
                                                FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}



void
FESystem::Structures::DKTPlate::initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe_tri3, FESystem::FiniteElement::FiniteElementBase& fe_tri6,
                                           const FESystem::Quadrature::QuadratureBase& q_tri3, const FESystem::Quadrature::QuadratureBase& q_tri6,
                                           FESystemDouble E, FESystemDouble nu, FESystemDouble rho, FESystemDouble th)
{
    FESystem::Structures::Structural2DElementBase::initialize(elem, fe_tri3, q_tri3, E, nu, rho, th);
    
    this->finite_element_tri6 = &fe_tri6;
    this->quadrature_tri6 = &q_tri6;
    
    // create the TRI6 element
    this->node0->copyVector(elem.getNode(0));
    this->node1->copyVector(elem.getNode(1));
    this->node2->copyVector(elem.getNode(2));
    this->node3->copyVector(elem.getNode(0)); this->node3->add(1.0, elem.getNode(1)); this->node3->scale(0.5);
    this->node4->copyVector(elem.getNode(1)); this->node4->add(1.0, elem.getNode(2)); this->node4->scale(0.5);
    this->node5->copyVector(elem.getNode(2)); this->node5->add(1.0, elem.getNode(0)); this->node5->scale(0.5);
    
    this->tri6_elem->setNode(0, *(this->node0));
    this->tri6_elem->setNode(1, *(this->node1));
    this->tri6_elem->setNode(2, *(this->node2));
    this->tri6_elem->setNode(3, *(this->node3));
    this->tri6_elem->setNode(4, *(this->node4));
    this->tri6_elem->setNode(5, *(this->node5));
}


void
FESystem::Structures::DKTPlate::calculateConsistentMassMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    const FESystemUInt n = this->finite_element->getNShapeFunctions();
    
    static FESystem::Numerics::DenseMatrix<FESystemDouble> B_mat, C_mat, tmp_mat1, tmp_mat2;
    C_mat.resize(3,3); B_mat.resize(3, 3*n); tmp_mat1.resize(3, 3*n), tmp_mat2.resize(3*n, 3*n);
    C_mat.zero(); B_mat.zero(); tmp_mat1.zero(); tmp_mat2.zero();
    
    const std::vector<FESystem::Geometry::Point*>& q_pts = this->quadrature->getQuadraturePoints();
    const std::vector<FESystemDouble>& q_weight = this->quadrature->getQuadraturePointWeights();
    
    FESystemDouble jac=0.0;
    mat.zero();
    this->getMaterialMassMatrix(C_mat);
    
    for (FESystemUInt i=0; i<q_pts.size(); i++)
    {
        jac = this->finite_element->getJacobianValue(*(q_pts[i]));
        this->calculateInertiaOperatorMatrix(*(q_pts[i]), B_mat);
        
        C_mat.matrixRightMultiply(1.0, B_mat, tmp_mat1);
        B_mat.matrixTransposeRightMultiply(1.0, tmp_mat1, tmp_mat2);
        
        mat.add(q_weight[i]*jac, tmp_mat2);
    }
}



void
FESystem::Structures::DKTPlate::calculateStiffnessMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    const FESystemUInt n = this->finite_element->getNShapeFunctions();
    
    static FESystem::Numerics::DenseMatrix<FESystemDouble> B_mat, C_mat, tmp_mat1, tmp_mat2;
    C_mat.resize(3,3); B_mat.resize(3, 3*n); tmp_mat1.resize(3, 3*n), tmp_mat2.resize(3*n, 3*n);
    C_mat.zero(); B_mat.zero(); tmp_mat1.zero(); tmp_mat2.zero();
    
    const std::vector<FESystem::Geometry::Point*>& q_pts = this->quadrature->getQuadraturePoints();
    const std::vector<FESystemDouble>& q_weight = this->quadrature->getQuadraturePointWeights();
    
    FESystemDouble jac=0.0;
    mat.zero();
    this->getMaterialComplianceMatrix(C_mat);
    
    for (FESystemUInt i=0; i<q_pts.size(); i++)
    {
        jac = this->finite_element->getJacobianValue(*(q_pts[i]));
        this->calculateBendingOperatorMatrix(*(q_pts[i]), B_mat);
        
        C_mat.matrixRightMultiply(1.0, B_mat, tmp_mat1);
        B_mat.matrixTransposeRightMultiply(1.0, tmp_mat1, tmp_mat2);
        
        mat.add(q_weight[i]*jac, tmp_mat2);
    }
}




void
FESystem::Structures::DKTPlate::calculateInertiaOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat)
{
    const FESystemUInt n = this->finite_element->getNShapeFunctions();
    static FESystem::Numerics::LocalVector<FESystemDouble> Nvec;
    Nvec.resize(n);
    B_mat.zero();
    
    Nvec.zero();
    this->finite_element->getShapeFunction(pt, Nvec);
    B_mat.setRowVals(0,   0,   n-1, Nvec); // w
    B_mat.setRowVals(1,   n, 2*n-1, Nvec); // theta_x
    B_mat.setRowVals(2, 2*n, 3*n-1, Nvec); // theta_y
}



void
FESystem::Structures::DKTPlate::calculateBendingOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat)
{
    const FESystemUInt n = this->finite_element->getNShapeFunctions();
    
    static std::vector<FESystemUInt> derivatives_x(2),derivatives_y(2);
    derivatives_x[0] = 1; derivatives_x[1] = 0;
    derivatives_y[0] = 0; derivatives_y[1] = 1;
    static FESystem::Numerics::LocalVector<FESystemDouble> Nvec, dbetaxdx, dbetaxdy, dbetaydx, dbetaydy;
    Nvec.resize(n); dbetaxdx.resize(9); dbetaxdy.resize(9); dbetaydx.resize(9); dbetaydy.resize(9);
    B_mat.zero(); dbetaxdx.zero(); dbetaxdy.zero(); dbetaydx.zero(); dbetaydy.zero();
    
    Nvec.zero();
    this->finite_element_tri6->getShapeFunctionDerivativeForPhysicalCoordinates(derivatives_x, pt, Nvec);
    FESystem::FEBatoz::calculateDKTShapeFunction(*tri6_elem, Nvec, dbetaxdx, dbetaydx);
    
    Nvec.zero();
    this->finite_element_tri6->getShapeFunctionDerivativeForPhysicalCoordinates(derivatives_y, pt, Nvec);
    FESystem::FEBatoz::calculateDKTShapeFunction(*tri6_elem, Nvec, dbetaxdy, dbetaydy);
    
    B_mat.setRowVals(0, 0,  8, dbetaxdx); // epsilon-x: phix
    B_mat.setRowVals(1, 0,  8, dbetaydy); // epsilon-y: phiy
    dbetaxdy.add(1.0, dbetaydx);
    B_mat.setRowVals(2, 0,  8, dbetaxdy); // gamma-xy : phix
    
}


void
FESystem::Structures::DKTPlate::getMaterialMassMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    
    FESystemAssert4(((s.first == 3) && (s.second== 3)), FESystem::Numerics::MatrixSizeMismatch, 3, 3, s.first, s.second);
    
    mat.setVal(0, 0, this->rho_val * this->th_val);
    mat.setVal(1, 1, 1.0e-12 * this->rho_val * this->th_val);
    mat.setVal(2, 2, 1.0e-12 * this->rho_val * this->th_val);
}



void
FESystem::Structures::DKTPlate::getMaterialComplianceMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    
    FESystemAssert4(((s.first == 3) && (s.second== 3)), FESystem::Numerics::MatrixSizeMismatch, 3, 3, s.first, s.second);
    
    FESystemDouble val = this->E_val/(1.0-this->nu_val*this->nu_val);
    
    mat.setVal(0, 0, val);
    mat.setVal(0, 1, this->nu_val*val);
    mat.setVal(1, 0, this->nu_val*val);
    mat.setVal(1, 1, val);
    mat.setVal(2, 2, this->G_val);
    mat.scale(pow(this->th_val,3)/12.0);
}




