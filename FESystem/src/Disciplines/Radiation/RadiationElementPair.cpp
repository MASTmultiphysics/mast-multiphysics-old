////
////  RadiationElementPair.cpp
////  FESystem
////
////  Created by Manav Bhatia on 1/25/13.
////
////
//
//// C++ includes
//
//
//// FESystem includes
//#include "Disciplines/Radiation/RadiationElementPair.h"
//#include "Geom/Point.h"
//#include "Numerics/LocalVector.h"
//#include "Numerics/DenseMatrix.h"
//#include "Mesh/Node.h"
//#include "Mesh/FaceElemBase.h"
//#include "FiniteElems/FELagrange.h"
//#include "Quadrature/TrapezoidQuadrature.h"
//#include "Base/FESystemExceptions.h"
//
//
//
//FESystem::Radiation::RadiationElementPair::RadiationElementPair():
//initialized(false),
//rad_elem1(NULL),
//rad_elem2(NULL),
//rad_elem1_orientation(0),
//rad_elem2_orientation(0)
//{
//}
//
//
//
//
//
//FESystem::Radiation::RadiationElementPair::~RadiationElementPair()
//{
//    
//}
//
//
//
//void
//FESystem::Radiation::RadiationElementPair::reinit(const FESystem::Mesh::FaceElemBase& el1, const FESystemInt orient1, const FESystem::Mesh::FaceElemBase& el2, const FESystemInt orient2)
//{
//    // this element should be initialized
//    FESystemAssert0(!this->initialized, FESystem::Exception::InvalidState);
//    
//    this->rad_elem1 = &el1;
//    this->rad_elem2 = &el2;
//    this->rad_elem1_orientation = orient1;
//    this->rad_elem2_orientation = orient2;
//    
//    this->initialized = true;
//}
//
//
//
//
//void
//FESystem::Radiation::RadiationElementPair::clear()
//{
//    // this element should be initialized
//    this->rad_elem1 = NULL;
//    this->rad_elem2 = NULL;
//    
//    this->initialized = false;
//}
//
//
//
//
//
//void
//FESystem::Radiation::RadiationElementPair::getShapeFactor(FESystem::Radiation::ShapeFactorMethod method, FESystemDouble& f1, FESystemDouble& f2)
//{
//    // return zero if the surface elements have 1st order geometry and don't see each other
//    if ((this->rad_elem1->getGeometryOrder() == 1) && (this->rad_elem2->getGeometryOrder() == 1))
//    {
//        FESystem::Numerics::LocalVector<FESystemDouble> normal1, normal2;
//        normal1.resize(3); normal2.resize(3);
//        this->rad_elem1->calculateSurfaceNormal(normal1);
//        this->rad_elem2->calculateSurfaceNormal(normal2);
//        
//        // since the cavity is assumed to be convex, if the normals are parallel,
//        // and if the dot product of the normals is not negative, then they lie in
//        // the same plane, and their shape factor should be zero.
//        if (fabs(normal1.dotProduct(normal2) - 1.0) <= FESystem::Base::getMachineEpsilon<FESystemDouble>())
//        {
//            f1 = 0.0; f2 = 0.0;
//            return;
//        }
//    }
//    
//    // otherwise, if the geometry order is > 1, or if the elems see each other, then calculate the factors
//    switch (method)
//    {
//        case DOUBLE_AREA_INT:
//            return this->doubleAreaIntShapeFactor(f1, f2);
//            break;
//            
//            
//        case CONTOUR_INT:
//            return this->contourIntShapeFactor(f1, f2);
//            break;
//            
//        case MITALAS_CONTOUR_INT:
//            return this->mitalasContourIntShapeFactor(f1, f2);
//            break;
//            
//        default:
//            FESystemAssert1(false, FESystem::Exception::EnumerationNotHandled, method);
//    }
//    
//}
//
//
//
//
//
//void
//FESystem::Radiation::RadiationElementPair::doubleAreaIntShapeFactor(FESystemDouble& f1, FESystemDouble& f2)
//{
//    FESystem::FiniteElement::FELagrange fe_elem1, fe_elem2;
//    FESystem::Quadrature::TrapezoidQuadrature q_rule;
//    q_rule.init(2, 10);
//    
//    FESystemDouble factor = 0.0, area1 = 0.0, area2 = 0.0, length = 0.0, cos_theta1 = 0.0, cos_theta2 = 0.0;
//    FESystemUInt n_nodes1 = this->rad_elem1->getNNodes(), n_nodes2 = this->rad_elem2->getNNodes();
//    
//    area1 = this->rad_elem1->getElementSize(fe_elem1, q_rule), area2 = this->rad_elem2->getElementSize(fe_elem2, q_rule);
//    
//    fe_elem1.clear();
//    fe_elem2.clear();
//    
//    fe_elem1.reinit(*this->rad_elem1);
//    fe_elem2.reinit(*this->rad_elem2);
//    
//    FESystem::Numerics::LocalVector<FESystemDouble> xy1, xy2, xy1_nodal, xy2_nodal, Nvec1, Nvec2, normal1, normal2;
//    FESystem::Numerics::DenseMatrix<FESystemDouble> B1, B2;
//    
//    xy1.resize(3); xy2.resize(3); Nvec1.resize(n_nodes1); Nvec2.resize(n_nodes2);
//    xy1_nodal.resize(3*n_nodes1); xy2_nodal.resize(3*n_nodes2); normal1.resize(3); normal2.resize(3);
//    B1.resize(3, 3*n_nodes1); B2.resize(3, 3*n_nodes2);
//    
//    const std::vector<FESystem::Geometry::Point*>& q_pts = q_rule.getQuadraturePoints();
//    const std::vector<FESystemDouble>& q_weight = q_rule.getQuadraturePointWeights();
//    
//    // get the location of the nodes of the two elements
//    // this is available from the local_elem_node pointers
//    
//    for (FESystemUInt n1=0; n1<n_nodes1; n1++)
//        for (FESystemUInt i=0; i<3; i++)
//            xy1_nodal.setVal(i*n_nodes1+n1, this->rad_elem1->getNode(n1).getVal(i));
//    
//    for (FESystemUInt n1=0; n1<n_nodes2; n1++)
//        for (FESystemUInt i=0; i<3; i++)
//            xy2_nodal.setVal(i*n_nodes2+n1, this->rad_elem2->getNode(n1).getVal(i));
//    
//    // iterate  over the quadrature points of both the elements
//    for (FESystemUInt qp1=0; qp1<q_pts.size(); qp1++)
//    {
//        // calculate the xy location on the first element
//        fe_elem1.getShapeFunction(*(q_pts[qp1]), Nvec1);
//        B1.zero();
//        for (FESystemUInt n1=0; n1<3; n1++)
//            B1.setRowVals(n1, n1*n_nodes1, (n1+1)*n_nodes1-1, Nvec1);
//        
//        B1.rightVectorMultiply(xy1_nodal, xy1);
//		
//        for (FESystemUInt qp2=0; qp2<q_pts.size(); qp2++)
//        {
//            // calculate the xy location on second element
//            B2.zero();
//            fe_elem2.getShapeFunction(*(q_pts[qp2]), Nvec2);
//            for (FESystemUInt n1=0; n1<3; n1++)
//                B2.setRowVals(n1, n1*n_nodes2, (n1+1)*n_nodes2-1, Nvec2);
//            
//            B2.rightVectorMultiply(xy2_nodal, xy2);
//            
//            // now calculate the factor
//            length = 0.0; cos_theta1 = 0.0; cos_theta2= 0.0;
//            
//            // calculate the vector from point 1 to 2
//            xy2.add(-1.0, xy1);
//            length = xy2.getL2Norm();
//            
//            // calculate the value of cos(theta) at point 1 and 2
//            if (length >  1.0e-12)
//            {
//                cos_theta1 = fabs(normal1.dotProduct(xy2)) / length;
//                cos_theta2 = fabs(normal2.dotProduct(xy2)) / length;
//                
//                factor += (cos_theta1 * cos_theta2)/(PI_VAL* length * length) * fe_elem1.getJacobianValue(*q_pts[qp1]) * fe_elem1.getJacobianValue(*q_pts[qp1]) * q_weight[qp1] * q_weight[qp2];
//            }
//        }
//    }
//    
//    if (factor < 0.0)
//        factor = 0.0;
//    
//    f1 = factor/area1;
//    f2 = factor/area2;
//}
//
//
//
//
//FESystemDouble mitalasEdgeFactor(const FESystem::Numerics::VectorBase<FESystemDouble>& line1_p1,
//                                 const FESystem::Numerics::VectorBase<FESystemDouble>& line1_q1,
//                                 const FESystem::Numerics::VectorBase<FESystemDouble>& line2_p2,
//                                 const FESystem::Numerics::VectorBase<FESystemDouble>& line2_q2)
//{
//    // calculate the lengths for this pair
//    FESystem::Numerics::LocalVector<FESystemDouble> unit_vec1, unit_vec2, u_point, du_point, temp_vec1,
//    p1, p2, q1, q2;
//    FESystemDouble length1, length2, integral, d1, d2, f1, f2, T, J, V, cos_lambda,
//    cos_theta, omega, du, factor = 0.0;
//    FESystemUInt sum  = 0;
//    const FESystemUInt n_divs = 10;
//    FESystemDouble cos_vec_in[2*n_divs], cos_vec_out[2*n_divs],
//    log_vec_in[2*n_divs], log_vec_out[2*n_divs], V_vec[n_divs];
//    
//    FESystemBoolean orientation_switch = false;
//    
//    unit_vec1.zero(); unit_vec2.zero(); u_point.zero(); du_point.zero(); temp_vec1.zero();
//    length1=0.0; length2=0.0; integral=0.0;
//    d1=0.0; d2=0.0; f1=0.0; f2=0.0;
//    T=0.0; J=0.0; V=0.0; cos_lambda=0.0; cos_theta=0.0; omega=0.0; du=0.0;
//    
//    // calculate the lengths
//    unit_vec1.copyVector(line1_q1); unit_vec1.add(-1.0, line1_p1); length1 = unit_vec1.getL2Norm();
//    unit_vec2.copyVector(line2_q2); unit_vec2.add(-1.0, line2_p2); length2 = unit_vec2.getL2Norm();
//    
//    // calculate the unit vectors
//    factor = 1.0/length1; unit_vec1.scale(factor);
//    factor = 1.0/length2; unit_vec2.scale(factor);
//    
//    // if the two vectors have opposite orientation, switch the point of
//    // one of the vectors and then proceede with the calculations. Finally, multiply the
//    // integral with -1.
//    p1.copyVector(line1_p1);
//    q1.copyVector(line1_q1);
//    if (unit_vec1.dotProduct(unit_vec2) <= 0.0)
//    {
//        p2.copyVector(line2_q2);
//        q2.copyVector(line2_p2);
//        unit_vec2.scale(-1.0);
//        orientation_switch = true;
//    }
//    else
//    {
//        p2.copyVector(line2_p2);
//        q2.copyVector(line2_q2);
//    }
//    
//    // calculate the unit vecs to find if the points p1,q1 lie on the other line segment
//    temp_vec1.copyVector(p1); temp_vec1.add(-1.0, p2); d1 = temp_vec1.getL2Norm();
//    if (d1 > FESystem::Base::getMachineEpsilon<FESystemDouble>())
//    {
//        factor = 1.0/d1; temp_vec1.scale(factor);
//        if (fabs(fabs(temp_vec1.dotProduct(unit_vec2)) -1) <= FESystem::Base::getMachineEpsilon<FESystemDouble>())
//            sum += 1;
//    }
//    else
//        sum += 1; // since the point coincides
//    
//    temp_vec1.copyVector(q1); temp_vec1.add(-1.0, p2); f1 = temp_vec1.getL2Norm();
//    if (f1 > FESystem::Base::getMachineEpsilon<FESystemDouble>())
//    {
//        factor = 1.0/f1; temp_vec1.scale(factor);
//        if (fabs(fabs(temp_vec1.dotProduct(unit_vec2)) -1) <= FESystem::Base::getMachineEpsilon<FESystemDouble>())
//            sum += 2;
//    }
//    else
//        sum += 2; // since the point coincides
//    
//    
//    temp_vec1.copyVector(p1); temp_vec1.add(-1.0, q2); d2 = temp_vec1.getL2Norm();
//    temp_vec1.copyVector(q1); temp_vec1.add(-1.0, q2); f2 = temp_vec1.getL2Norm();
//    
//    integral = 0.0;
//    du = length2 / (1.0 * n_divs);
//    // this is the starting point for intgration
//    u_point.copyVector(p2);
//    
//    // this is the increment in the integration point
//    // this is scaled to half the value in order to increment the initial
//    // point. Later is it increased to the actual value through multiplication by 2.0
//    du_point.copyVector(unit_vec2);
//    du_point.scale(du/2.0);
//    
//    // now increment the point to the middle of the first interval
//    u_point.add(1.0,du_point);
//    du_point.scale(2.0);
//    
//    switch (sum)
//    {
//        case 0:
//        {
//            // calculate all terms by numerical integration
//            for (FESystemUInt i=0; i < n_divs; i++)
//            {
//                // calculate T, J, V, cos lambda, cos theta, omega
//                temp_vec1.copyVector(p1); temp_vec1.add(-1.0, u_point);
//                J = temp_vec1.getL2Norm();
//                factor = 1.0/J; temp_vec1.scale(factor);
//                cos_theta = (-unit_vec1.dotProduct(temp_vec1));
//                if (fabs(cos_theta) > 1.0)
//                {
//                    if (cos_theta < 0.0)
//                        cos_theta = -1.0;
//                    else
//                        cos_theta = 1.0;
//                }
//                
//                // calculate the point for calculation of V
//                temp_vec1.copyVector(p1);
//                temp_vec1.add(cos_theta * J, unit_vec1);
//                temp_vec1.add(-1.0, u_point);
//                V = temp_vec1.getL2Norm();
//                
//                // repeat as above for T and cos_lambda
//                temp_vec1.copyVector(q1); temp_vec1.add(-1.0, u_point);
//                T = temp_vec1.getL2Norm();
//                factor = 1.0/T; temp_vec1.scale(factor);
//                cos_lambda = (unit_vec1.dotProduct(temp_vec1));
//                if (fabs(cos_lambda) > 1.0)
//                {
//                    if (cos_lambda < 0.0)
//                        cos_lambda = -1.0;
//                    else
//                        cos_lambda = 1.0;
//                }
//                
//                // cos theta starts at 0, and cos lambda starts at n_divs
//                cos_vec_in[i] = cos_lambda;
//                cos_vec_in[i+n_divs] = cos_theta;
//                
//                // T starts at 0, and J starts at n_divs
//                log_vec_in[i] = T;
//                log_vec_in[i+n_divs] = J;
//                V_vec[i] = V;
//                // increment the point
//                u_point.add(1.0, du_point);
//            }
//            
//            
//#ifdef ENABLE_ACCELERATE_FRAMEWORK
//            // calculate the cos values
//            int_val = 2*n_divs;
//            vvacos(cos_vec_out, cos_vec_in, &int_val);
//            
//            // calculate the log values
//            vvlog(log_vec_out, log_vec_in, &int_val);
//#else
//            for (FESystemUInt i=0; i<2*n_divs; i++)
//            {
//                cos_vec_out[i] = acos(cos_vec_in[i]);
//                log_vec_out[i] = log(log_vec_in[i]);
//            }
//#endif  // ENABLE_ACCELERATE_FRAMEWORK
//            
//            // calculate the integral values
//            for (FESystemUInt i=0; i<n_divs; i++)
//            {
//                //
//                // the original statement which has been vectorized is
//                // omega = FESystemNumbers::Pi - (acos(cos_lambda) + acos(cos_theta));
//                //
//                omega = PI_VAL - (cos_vec_out[i] + cos_vec_out[i+n_divs]);
//                
//                //
//                // the original statement which has been vectorized is
//                // integral += (T * cos_lambda * log(T) + J * cos_theta * log(J) + V * omega) * du;
//                //
//                integral += (log_vec_in[i] * cos_vec_in[i] * log_vec_out[i] +
//                             log_vec_in[i+n_divs] * cos_vec_in[i+n_divs] * log_vec_out[i+n_divs] +
//                             V_vec[i] * omega) * du;
//            }
//        }
//            break;
//            
//        case 1:
//        {
//            // direct calculation of cos theta term
//            cos_theta = unit_vec1.dotProduct(unit_vec2);
//            if (fabs(cos_theta) > 1.0)
//            {
//                if (cos_theta < 0.0)
//                    cos_theta = -1.0;
//                else
//                    cos_theta = 1.0;
//            }
//            if (d2 > 0.0)
//                integral += cos_theta * (d2*d2 * (0.5 * log(d2) - 0.25));
//            if (d1 > 0.0)
//                integral -= cos_theta * (d1*d1 * (0.5 * log(d1) - 0.25));
//            
//            // if the point p1 is on q2 or beyond it on unit_vec2, then
//            // cos theta needs to be scaled by -1.0
//            if ((d2 == 0.0 || d1 > d2) &&
//                ((d1 + d2 - length2) >= 0.0))
//                cos_theta *= -1.0;
//            
//            // if p1 does not lie between the two points p2 and q2, the perform the
//            // integration, else break the line segment
//            if ((d1 + d2 - length2) >= 0.0)
//            {
//                // calculate all other terms by numerical integration
//                for (FESystemUInt i=0; i < n_divs; i++)
//                {
//                    // calculate T, V, cos lambda, omega
//                    temp_vec1.copyVector(q1); temp_vec1.add(-1.0, u_point);
//                    T = temp_vec1.getL2Norm();
//                    factor = 1.0/T; temp_vec1.scale(factor);
//                    cos_lambda = (unit_vec1.dotProduct(temp_vec1));
//                    if (fabs(cos_lambda) > 1.0)
//                    {
//                        if (cos_lambda < 0.0)
//                            cos_lambda = -1.0;
//                        else
//                            cos_lambda = 1.0;
//                    }
//                    
//                    // calculate the point for calculation of V
//                    temp_vec1.copyVector(q1);
//                    temp_vec1.add(-cos_lambda * T, unit_vec1);
//                    temp_vec1.add(-1.0, u_point);
//                    V = temp_vec1.getL2Norm();
//                    
//                    // cos theta starts at 0, and cos lambda starts at n_divs
//                    cos_vec_in[i] = cos_lambda;
//                    
//                    // T starts at 0, and J starts at n_divs
//                    log_vec_in[i] = T;
//                    V_vec[i] = V;
//                    
//                    // increment the point
//                    u_point.add(1.0, du_point);
//                }
//                
//                // cos_theta is constant, hence, needs to be evaluated only once
//                cos_vec_in[n_divs] = cos_theta;
//                
//#ifdef ENABLE_ACCELERATE_FRAMEWORK
//                // calculate the cos values
//                int_val = 1+n_divs;
//                vvacos(cos_vec_out, cos_vec_in, &int_val);
//                
//                // calculate the log values
//                int_val = n_divs;
//                vvlog(log_vec_out, log_vec_in, &int_val);
//#else
//                for (FESystemUInt i=0; i<n_divs; i++)
//                {
//                    cos_vec_out[i] = acos(cos_vec_in[i]);
//                    log_vec_out[i] = log(log_vec_in[i]);
//                }
//                cos_vec_out[n_divs] = acos(cos_vec_in[n_divs]);
//#endif // ENABLE_ACCELERATE_FRAMEWORK
//                
//                // calculate the integral values
//                for (FESystemUInt i=0; i<n_divs; i++)
//                {
//                    //
//                    // the original statement which has been vectorized is
//                    // omega = FESystemNumbers::Pi - (acos(cos_lambda) + acos(cos_theta));
//                    //
//                    omega = PI_VAL - (cos_vec_out[i] + cos_vec_out[n_divs]);
//                    
//                    //
//                    // the original statement which has been vectorized is
//                    // integral += (T * cos_lambda * log(T) + V * omega) * du;
//                    //
//                    integral += (log_vec_in[i] * cos_vec_in[i] * log_vec_out[i] +
//                                 V_vec[i] * omega) * du;
//                }
//            }
//            else
//            {
//                // needs to be implemented
//                FESystemAssert0(false, FESystem::Exception::InvalidValue);
//            }
//        }
//            break;
//            
//        case 2:
//        {
//            // direct calculation of cos lambda term
//            cos_lambda = unit_vec1.dotProduct(unit_vec2);
//            if (fabs(cos_lambda) > 1.0)
//            {
//                if (cos_lambda < 0.0)
//                    cos_lambda = -1.0;
//                else
//                    cos_lambda = 1.0;
//            }
//            
//            if (f2 > 0.0)
//                integral -= cos_lambda * (f2*f2 * (0.5 * log(f2) - 0.25));
//            if (f1 > 0.0)
//                integral += cos_lambda * (f1*f1 * (0.5 * log(f1) - 0.25));
//            
//            // if the point p1 is on q2 or beyond it on unit_vec2, then
//            // cos theta needs to be scaled by -1.0
//            if ((f1 == 0.0 || f2 > f1) &&
//                ((f1 + f2 - length2) >= 0.0))
//                cos_lambda *= -1.0;
//            
//            // if p1 does not lie between the two points p2 and q2, the perform the
//            // integration, else break the line segment
//            if ((f1 + f2 - length2) >= 0.0)
//            {
//                // calculate all terms by numerical integration
//                for (FESystemUInt i=0; i < n_divs; i++)
//                {
//                    // calculate J, V, cos theta, omega
//                    temp_vec1.copyVector(p1); temp_vec1.add(-1.0, u_point);
//                    J = temp_vec1.getL2Norm();
//                    factor = 1.0/J; temp_vec1.scale(factor);
//                    cos_theta = (-unit_vec1.dotProduct(temp_vec1));
//                    if (fabs(cos_theta) > 1.0)
//                    {
//                        if (cos_theta < 0.0)
//                            cos_theta = -1.0;
//                        else
//                            cos_theta = 1.0;
//                    }
//                    
//                    // calculate the point for calculation of V
//                    temp_vec1.copyVector(p1);
//                    temp_vec1.add(cos_theta * J, unit_vec1);
//                    temp_vec1.add(-1.0, u_point);
//                    V = temp_vec1.getL2Norm();
//                    
//                    // cos theta starts at 0, and cos lambda starts at n_divs
//                    cos_vec_in[1+i] = cos_theta;
//                    
//                    // T starts at 0, and J starts at n_divs
//                    log_vec_in[i] = J;
//                    V_vec[i] = V;
//                    // increment the point
//                    u_point.add(1.0, du_point);
//                }
//                
//                // cos_lambda is constant, hence is evaluated only once
//                cos_vec_in[0] = cos_lambda;
//                
//#ifdef ENABLE_ACCELERATE_FRAMEWORK
//                // calculate the cos values
//                int_val = 1+n_divs;
//                vvacos(cos_vec_out, cos_vec_in, &int_val);
//                
//                // calculate the log values
//                int_val = n_divs;
//                vvlog(log_vec_out, log_vec_in, &int_val);
//#else
//                for (FESystemUInt i=0; i<n_divs; i++)
//                {
//                    cos_vec_out[i] = acos(cos_vec_in[i]);
//                    log_vec_out[i] = log(log_vec_in[i]);
//                }
//                cos_vec_out[n_divs] = acos(cos_vec_in[n_divs]);
//#endif // ENABLE_ACCELERATE_FRAMEWORK
//                
//                // calculate the integral values
//                for (FESystemUInt i=0; i<n_divs; i++)
//                {
//                    //
//                    // the original statement which has been vectorized is
//                    // omega = FESystemNumbers::Pi - (acos(cos_lambda) + acos(cos_theta));
//                    //
//                    omega = PI_VAL - (cos_vec_out[0] + cos_vec_out[1+i]);
//                    
//                    //
//                    // the original statement which has been vectorized is
//                    // integral += (J * cos_theta * log(J) + V * omega) * du;
//                    //
//                    integral += (log_vec_in[i] * cos_vec_in[1+i] * log_vec_out[i] +
//                                 V_vec[i] * omega) * du;
//                }
//            }
//            else
//            {
//                FESystemAssert0(false, FESystem::Exception::InvalidValue);
//            }
//        }
//            break;
//            
//        case 3:
//        {
//            // direct calculation of all terms
//            if (d2 > 0.0)
//                integral += (d2*d2 * (0.5 * log(d2) - 0.25));
//            if (d1 > 0.0)
//                integral -= (d1*d1 * (0.5 * log(d1) - 0.25));
//            if (f2 > 0.0)
//                integral -= (f2*f2 * (0.5 * log(f2) - 0.25));
//            if (f1 > 0.0)
//                integral += (f1*f1 * (0.5 * log(f1) - 0.25));
//            
//        }
//            break;
//            
//        default:
//            FESystemAssert0(false, FESystem::Exception::InvalidValue);
//    }
//    
//    integral -= length1 * length2;
//    integral *= (unit_vec1.dotProduct(unit_vec2));
//    if (orientation_switch)
//        integral *= -1.0;
//    
//    return integral;
//}
//
//
//FESystemDouble NumIntEdgeFactor(const FESystem::Numerics::VectorBase<FESystemDouble>& line1_p1,
//                                const FESystem::Numerics::VectorBase<FESystemDouble>& line1_q1,
//                                const FESystem::Numerics::VectorBase<FESystemDouble>& line2_p2,
//                                const FESystem::Numerics::VectorBase<FESystemDouble>& line2_q2)
//{
//    // calculate the lengths for this pair
//    FESystem::Numerics::LocalVector<FESystemDouble> unit_vec1, unit_vec2, u_point, du_point, e_point, de_point_b2,
//    de_point, temp_vec1,
//    p1, p2, q1, q2;
//    FESystemDouble length1, length2, integral, d1, d2, f1, f2, T, J, V, cos_lambda,
//    cos_theta, omega, du;
//    const FESystemUInt n_divs = 50;
//    FESystemBoolean orientation_switch = false;
//    
//    unit_vec1.zero(); unit_vec2.zero(); u_point.zero(); du_point.zero();
//    e_point.zero(); temp_vec1.zero(); de_point.zero();
//    length1=0.0; length2=0.0; integral=0.0;
//    d1=0.0; d2=0.0; f1=0.0; f2=0.0;
//    T=0.0; J=0.0; V=0.0; cos_lambda=0.0; cos_theta=0.0; omega=0.0; du=0.0;
//    
//    // calculate the lengths
//    unit_vec1.copyVector(line1_q1); unit_vec1.add(-1.0, line1_p1); length1 = unit_vec1.getL2Norm();
//    unit_vec2.copyVector(line2_q2); unit_vec2.add(-1.0, line2_p2); length2 = unit_vec2.getL2Norm();
//    
//    // calculate the unit vectors
//    unit_vec1.scale(1.0/length1);
//    unit_vec2.scale(1.0/length2);
//    
//    // if the two vectors have opposite orientation, switch the point of
//    // one of the vectors and then proceede with the calculations. Finally, multiply the
//    // integral with -1.
//    p1.copyVector(line1_p1);
//    q1.copyVector(line1_q1);
//    if (unit_vec1.dotProduct(unit_vec2) <= 0.0)
//    {
//        p2.copyVector(line2_q2);
//        q2.copyVector(line2_p2);
//        unit_vec2.scale(-1.0);
//        orientation_switch = true;
//    }
//    else
//    {
//        p2.copyVector(line2_p2);
//        q2.copyVector(line2_q2);
//    }
//    
//    integral = 0.0;
//    du = length2 / (1.0 * n_divs);
//    // this is the starting point for intgration
//    u_point.copyVector(p2);
//    
//    // this is the increment in the integration point
//    // this is scaled to half the value in order to increment the initial
//    // point. Later is it increased to the actual value through multiplication by 2.0
//    du_point.copyVector(unit_vec2);
//    du_point.scale(du/2.0);
//    
//    de_point.copyVector(unit_vec1);
//    de_point.scale(du);
//    de_point_b2.copyVector(unit_vec1);
//    de_point_b2.scale(du/2);
//    
//    // now increment the point to the middle of the first interval
//    u_point.add(1.0, du_point);
//    du_point.scale(2.0);
//    
//    // calculate all terms by numerical integration
//    for (FESystemUInt i=0; i < n_divs; i++)
//    {
//        
//        e_point.copyVector(p1);
//        e_point.add(1.0, de_point_b2);
//        for (FESystemUInt i=0; i < n_divs; i++)
//        {
//            temp_vec1.copyVector(e_point);
//            temp_vec1.add(-1.0, u_point);
//            
//            integral += log(temp_vec1.getL2Norm()) * du * du;
//            
//            e_point.add(1.0, de_point);
//        }
//        // once all the integration is done, increment the point
//        u_point.add(1.0, du_point);
//    }
//    
//    integral *= (unit_vec1.dotProduct(unit_vec2));
//    if (orientation_switch)
//        integral *= -1.0;
//    
//    return integral;
//}
//
//
//FESystemDouble	F1(const FESystemDouble &s0,
//                   const FESystemDouble &s1,
//                   const FESystemDouble &s2,
//                   const FESystemDouble &x)
//{
//    FESystemDouble b,f;
//    b = sqrt(4.0 * s0 * s2 - s1 * s1);
//    f = -4.0 * s2 * x;
//    f += 2.0 * b * atan2((s1 + 2.0 * s2 * x), b);
//    f += (s1 + 2.0 * s2 * x) * log(s0 + x * (s1 + s2 * x));
//    f /= (2.0 * s2);
//    return f;
//}
//
//
//FESystemDouble	F2(const FESystemDouble &s0,
//                   const FESystemDouble &s1,
//                   const FESystemDouble &s2,
//                   const FESystemDouble &x)
//{
//    FESystemDouble b,f;
//    b = sqrt(4.0 * s0 * s2 - s1 * s1);
//    f = 2.0 * (-s1) * b * b * atan2((s1 + 2.0 * s2 * x), b);
//    f += b * (2.0 * s2 * x * (s1 - s2 * x) +
//              (2.0 * s2 * (s2 * x * x + s0) - s1 * s1) * log(s0 + x * (s1 + s2 * x)));
//    f /= (4.0 * s2 * s2 * b);
//    return f;
//}
//
//
//FESystemDouble	F3(const FESystemDouble *b,
//                   const FESystemDouble &length1,
//                   const FESystemDouble &length2)
//{
//    FESystemDouble e[5],f,tmp1,tmp2,x;
//    
//    e[0] = 4.0 * b[3] * b[0] - b[1] * b[1];
//    e[1] = 4.0 * b[3] * b[2] - 2.0 * b[1] * b[5];
//    e[2] = 4.0 * b[3] * b[4] - b[5];
//    e[3] = b[1] + 2.0 * b[3] * length1;
//    e[4] = b[5];
//    
//    // use trapezoidal rule to integrate
//    FESystemUInt ndivs = 10;
//    FESystemDouble dx = length2/ndivs;
//    
//    f = 0.0;
//    // this evaluation takes care of both the upper and the lower limits of the
//    // integral. Hence, the atan2 function appears twice
//    for (FESystemUInt i=0; i < ndivs; i++)
//    {
//        x = dx*(i+0.5);
//        tmp1 = sqrtf(e[0] + e[1] * x + e[2] * x * x);
//        tmp2 = e[3] + e[4] * x;
//        f += dx * (tmp1 / b[3]) * (atan2(tmp2, tmp1) - atan2((b[1] + b[5] * x), tmp1));
//    }
//    
//    return f;
//}
//
//
//
//FESystemDouble contourIntegrationEdgeFactor(const FESystem::Numerics::LocalVector<FESystemDouble>& p1,
//                                            const FESystem::Numerics::LocalVector<FESystemDouble>& q1,
//                                            const FESystem::Numerics::LocalVector<FESystemDouble>& p2,
//                                            const FESystem::Numerics::LocalVector<FESystemDouble>& q2)
//{
//    // calculate the lengths for this pair
//    FESystem::Numerics::LocalVector<FESystemDouble> unit_vec1, unit_vec2, vec12;
//    FESystemDouble length1, length2, integral, upper_integral, lower_integral, b[6], d[5];
//    
//    // calculate the lengths
//    unit_vec1.copyVector(q1); unit_vec1.add(-1.0, p1); length1 = unit_vec1.getL2Norm();
//    unit_vec2.copyVector(q2); unit_vec2.add(-1.0, p2); length2 = unit_vec2.getL2Norm();
//    vec12.copyVector(p2); vec12.add(-1.0, p1);
//    
//    // calculate the unit vectors
//    unit_vec1.scale(1.0/length1);
//    unit_vec2.scale(1.0/length2);
//    
//    // zero the entries
//    for (FESystemUInt i=0; i<6; i++)
//        b[i] = 0.0;
//    
//    for (FESystemUInt i=0; i<3; i++)
//    {
//        b[0] += vec12.getVal(i) * vec12.getVal(i);
//        b[1] -= 2.0 * vec12.getVal(i) * unit_vec1.getVal(i);
//        b[2] += 2.0 * vec12.getVal(i) * unit_vec2.getVal(i);
//        b[3] += unit_vec1.getVal(i) * unit_vec1.getVal(i);
//        b[4] += unit_vec2.getVal(i) * unit_vec2.getVal(i);
//        b[5] -= 2.0 * unit_vec1.getVal(i) * unit_vec2.getVal(i);
//    }
//    
//    // now the d values
//    d[0] = b[0] + length1 * (b[1] + length1 * b[3]);
//    d[1] = b[2] + length1 * b[5];
//    d[2] = b[4];
//    d[3] = b[1] + 2.0 * b[3] * length1;
//    d[4] = b[5];
//    
//    
//    upper_integral =
//    -2.0 * length1 * length2 + // term1
//    (0.5 / b[3]) * (d[3] * F1(d[0], d[1], d[2], length2) +
//                    d[4] * F2(d[0], d[1], d[2], length2)) ;
//    
//    lower_integral =
//    (0.5 / b[3]) * (b[1] * F1(b[0], b[2], b[4], length2) +
//                    b[5] * F2(b[0], b[2], b[4], length2));
//    
//    integral = (0.25 / PI_VAL) * (upper_integral - lower_integral + F3(b, length1, length2));
//    return integral;
//}
//
//
//
//
//
//
//void
//FESystem::Radiation::RadiationElementPair::mitalasContourIntShapeFactor(FESystemDouble& f1, FESystemDouble& f2)
//{
//    FESystemDouble factor;
//    FESystemUInt n_sides1, n_sides2, j1, j2;
//    factor = 0.0;
//    
//    n_sides1 = this->rad_elem1->getNBoundaries();
//    n_sides2 = this->rad_elem2->getNBoundaries();
//    
//    for (FESystemUInt i=0; i < n_sides1; i++)
//    {
//        // get the nodes for this side
//        if (i == (n_sides1-1))
//            j1 = 0;
//        else
//            j1 = i + 1;
//        
//        const FESystem::Numerics::VectorBase<FESystemDouble>& p1 = this->rad_elem1->getNode(i);
//        const FESystem::Numerics::VectorBase<FESystemDouble>& q1 = this->rad_elem1->getNode(j1);
//        
//        for (FESystemUInt j=0; j < n_sides2; j++)
//        {
//            if (j == (n_sides2-1))
//                j2 = 0;
//            else
//                j2 = j + 1;
//            
//            const FESystem::Numerics::VectorBase<FESystemDouble>& p2 = this->rad_elem2->getNode(j);
//            const FESystem::Numerics::VectorBase<FESystemDouble>& q2 = this->rad_elem2->getNode(j2);
//            
//            factor += mitalasEdgeFactor(p1, q1, p2, q2);
//        }
//    }
//    
//    if (this->rad_elem1_orientation < 0)
//        factor *= -1.0;
//    if (this->rad_elem2_orientation < 0)
//        factor *= -1.0;
//    
//    factor /= (2.0 * PI_VAL);
//    
//
//    FESystem::FiniteElement::FELagrange fe;
//    FESystem::Quadrature::TrapezoidQuadrature q_rule;
//    q_rule.init(2, 2);
//    
//    FESystemDouble area1 = this->rad_elem1->getElementSize(fe, q_rule);
//    FESystemDouble area2 = this->rad_elem2->getElementSize(fe, q_rule);
//    
//    if (factor < 0.0)
//        factor = 0.0;
//    
//    f1 = factor/area1;
//    f2 = factor/area2;
//}
//
//
//
