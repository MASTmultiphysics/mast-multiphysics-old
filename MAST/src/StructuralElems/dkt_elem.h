//
//  dkt_elem.h
//  MAST
//
//  Created by Manav Bhatia on 10/24/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_dkt_elem_h__
#define __MAST_dkt_elem_h__

// libMesh includes
#include "libmesh/elem.h"
#include "libmesh/face_tri3.h"
#include "libmesh/face_tri6.h"
#include "libmesh/node.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature.h"


namespace MAST
{
    class DKTElem {
    public:
        DKTElem(const Elem& elem, QBase& qrule):
        _global_elem(elem),
        _qrule(qrule),
        _fe(NULL),
        _tri6(NULL),
        _nodes(3)
        {
            libmesh_assert(elem.type() == TRI3);
            
            _tri6 = new Tri6;
            
            // next three nodes are created using the corner nodes
            for (unsigned int i=0; i<3; i++)
                _nodes[i] = new Node;
            
            // first edge
            Point p;
            p = *elem.get_node(0); p += *elem.get_node(1); p*= 0.5;
            (*_nodes[0]) = p;
            _nodes[0]->set_id(3);
            
            // second edge
            p = *elem.get_node(1); p += *elem.get_node(2); p*= 0.5;
            (*_nodes[1]) = p;
            _nodes[1]->set_id(4);
            
            // third edge
            p = *elem.get_node(2); p += *elem.get_node(0); p*= 0.5;
            (*_nodes[2]) = p;
            _nodes[2]->set_id(5);

            // first three nodes are same as that of the original element
            for (unsigned int i=0; i<3; i++) {
                _tri6->set_node(  i) = elem.get_node(i);
                _tri6->set_node(i+3) = _nodes[i];
            }

            // now setup the shape functions
            _fe = FEBase::build(2, FEType(SECOND, LAGRANGE)).release();
            _fe->attach_quadrature_rule(&qrule);
            _fe->get_phi();
            _fe->get_dphi();
            _fe->reinit(_tri6);
        }
        
        
        /*!
         *   destructor
         */
        ~DKTElem() {
            delete _fe;
            delete _tri6;
            for (unsigned int i=0; i<3; i++)
                delete _nodes[i];
        }
        
        /*!
         *   initialze the bending strain operator for DKT element
         */
        void initialize_dkt_bending_strain_operator (const unsigned int qp,
                                                     FEMOperatorMatrix& Bmat);
        
    protected:
        
        /*!
         *   returns the length of the side defined by vector from node i to j
         */
        Real _get_edge_length(unsigned int i, unsigned int j);
        
        /*!
         *   returns the cos of normal to the side defined by vector from node i to j
         */
        void _get_edge_normal_sine_cosine(unsigned int i, unsigned int j,
                                         Real& sine, Real& cosine);
        
        
        void _calculate_dkt_shape_functions(const DenseVector<Real>& phi,
                                            DenseVector<Real>& betax,
                                            DenseVector<Real>& betay);
        
        /*!
         *    element for which DKT is created
         */
        const Elem& _global_elem;
        
        /*!
         *   quadrature rule to be used. This should already be of the correct order
         */
        QBase& _qrule;
        
        /*!
         *   FE object to get shape functions for this element
         */
        FEBase* _fe;
        
        /*!
         *   6 noded triangle that is used to calculate the shape functions
         */
        Elem* _tri6;
        
        /*!
         *    vector to store node pointers for the mid-sized nodes
         */
        std::vector<Node*> _nodes;
        
        
    };
}


inline
void
MAST::DKTElem::initialize_dkt_bending_strain_operator (const unsigned int qp,
                                                       FEMOperatorMatrix& Bmat) {
    
    const std::vector<std::vector<RealVectorValue> >& dphi = _fe->get_dphi();
    const unsigned int n_phi = (unsigned int)dphi.size();

    DenseVector<Real> phi, dbetaxdx, dbetaxdy, dbetaydx, dbetaydy,
    w, thetax, thetay;
    phi.resize(n_phi); dbetaxdx.resize(9); dbetaxdy.resize(9);
    dbetaydx.resize(9); dbetaydy.resize(9); w.resize(3); thetax.resize(3);
    thetay.resize(3);
    
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
        phi(i_nd) = dphi[i_nd][qp](0);  // dphi/dx
    _calculate_dkt_shape_functions(phi, dbetaxdx, dbetaydx);

    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
        phi(i_nd) = dphi[i_nd][qp](1);  // dphi/dy
    _calculate_dkt_shape_functions(phi, dbetaxdy, dbetaydy);
    
    // add the shear components together
    dbetaxdy.add(1.0, dbetaydx);


    // now set values for individual dofs
    for (unsigned int i=0; i<3; i++) {
        w(i)      = dbetaxdx(  i);
        thetax(i) = dbetaxdx(3+i);
        thetay(i) = dbetaxdx(6+i);
    }
    Bmat.set_shape_function(0, 2,      w); // epsilon-x: w
    Bmat.set_shape_function(0, 3, thetax); // epsilon-x: tx
    Bmat.set_shape_function(0, 4, thetay); // epsilon-x: ty

    
    for (unsigned int i=0; i<3; i++) {
        w(i)      = dbetaydy(  i);
        thetax(i) = dbetaydy(3+i);
        thetay(i) = dbetaydy(6+i);
    }
    Bmat.set_shape_function(1, 2,      w); // epsilon-x: w
    Bmat.set_shape_function(1, 3, thetax); // epsilon-x: tx
    Bmat.set_shape_function(1, 4, thetay); // epsilon-x: ty

    for (unsigned int i=0; i<3; i++) {
        w(i)      = dbetaxdy(  i);
        thetax(i) = dbetaxdy(3+i);
        thetay(i) = dbetaxdy(6+i);
    }
    Bmat.set_shape_function(2, 2,      w); // epsilon-x: w
    Bmat.set_shape_function(2, 3, thetax); // epsilon-x: tx
    Bmat.set_shape_function(2, 4, thetay); // epsilon-x: ty
}



inline
Real
MAST::DKTElem::_get_edge_length(unsigned int i, unsigned int j)
{
    Point l = *_tri6->get_node(j);
    l -= *_tri6->get_node(i);
    
    return l.size();
}




inline
void
MAST::DKTElem::_get_edge_normal_sine_cosine(unsigned int i, unsigned int j,
                                           Real& sine, Real& cosine)
{
    Point vec0, vec1, vec2, vec3;
    
    // calculate the normal to the element
    vec0 = *_global_elem.get_node(0);
    vec1 = *_global_elem.get_node(1);
    vec2 = *_global_elem.get_node(2);
    
    vec1 -= vec0;
    vec2 -= vec0;
    vec3 = vec1.cross(vec2);
    // this is the unit vector normal to the plane of the triangle
    vec1 = vec3;
    vec1 /= vec1.size();
    
    // cross product of the length vector with the surface normal will
    // give the needed vector
    vec3 = *_global_elem.get_node(i);
    vec2 = *_global_elem.get_node(j);
    
    vec2 -= vec3;
    vec3 = vec2.cross(vec1);
    vec3 /= vec3.size();        // this is the unit vector needed
    
    // cos of angle between this and the x-axis is simply the
    // 0th component of this vector
    cosine = vec3(0);
    sine   = vec3(1);
}




inline
void
MAST::DKTElem::_calculate_dkt_shape_functions(const DenseVector<Real>& phi,
                                              DenseVector<Real>& betax,
                                              DenseVector<Real>& betay)
{
    // -- keep in mind that the index numbers for the elems start at 0.
    // -- also, the mid side node numbers in the Batoz's paper are different from
    //   the ones used in this library. Hence, use the following association
    //                   BATOZ TRI6 node #              MAST TRI6 node #
    //                          1                               0
    //                          2                               1
    //                          3                               2
    //                          4 (on edge 2,3)                 4 (on edge 1,2)
    //                          5 (on edge 3,1)                 5 (on edge 2,0)
    //                          6 (on edge 1,2)                 3 (on edge 0,1)
    // -- all shape functions for this element are in a variable major format. Hence,
    //  they follow Hx for w1, w2, w3, thetax1, thetax2, thetax3, thetay1, thetay2, thetay3.
    // And then the same this is followed for Hy (from dofs 9-17)
    
    // local variables for shape functions
    Real N1, N2, N3, N4, N5, N6;
    
    // local variables for edge lengths and sine/cosines
    Real l12, l23, l31, cos4, cos5, cos6, sin4, sin5, sin6;
    
    N1 = phi(0);
    N2 = phi(1);
    N3 = phi(2);
    N4 = phi(4);
    N5 = phi(5);
    N6 = phi(3);
    
    _get_edge_normal_sine_cosine(1, 2, sin4, cos4);
    _get_edge_normal_sine_cosine(2, 0, sin5, cos5);
    _get_edge_normal_sine_cosine(0, 1, sin6, cos6);
    
    l12 = _get_edge_length(0, 1);
    l23 = _get_edge_length(1, 2);
    l31 = _get_edge_length(2, 0);
    
    betax(0) =  1.5 * (N5 * sin5 / l31 - N6 * sin6 / l12); // Hx, w1
    betax(1) =  1.5 * (N6 * sin6 / l12 - N4 * sin4 / l23); // Hx, w2
    betax(2) =  1.5 * (N4 * sin4 / l23 - N5 * sin5 / l31); // Hx, w3
    betax(3) =  (-.75) * (N5 * sin5 * cos5 + N6 * sin6 * cos6); // Hx, thetax1
    betax(4) =  (-.75) * (N4 * sin4 * cos4 + N6 * sin6 * cos6); // Hx, thetax2
    betax(5) =  (-.75) * (N5 * sin5 * cos5 + N4 * sin4 * cos4); // Hx, thetax3
    betax(6) =  (N1 + N5 * (0.5 * cos5 * cos5 - 0.25 * sin5 * sin5) + N6 * (0.5 * cos6 * cos6 - 0.25 * sin6 * sin6)); // Hx, thetay1
    betax(7) =  (N2 + N4 * (0.5 * cos4 * cos4 - 0.25 * sin4 * sin4) + N6 * (0.5 * cos6 * cos6 - 0.25 * sin6 * sin6)); // Hx, thetay2
    betax(8) =  (N3 + N5 * (0.5 * cos5 * cos5 - 0.25 * sin5 * sin5) + N4 * (0.5 * cos4 * cos4 - 0.25 * sin4 * sin4)); // Hx, thetay3
    
    betay(0) =  1.5 * (-N5 * cos5 / l31 + N6 * cos6 / l12); // Hy, w1
    betay(1) =  1.5 * (-N6 * cos6 / l12 + N4 * cos4 / l23); // Hy, w2
    betay(2) =  1.5 * (-N4 * cos4 / l23 + N5 * cos5 / l31); // Hy, w3
    betay(3) =  (-N1 + N5 * (0.25 * cos5 * cos5 - 0.5 * sin5 * sin5) + N6 * (0.25 * cos6 * cos6 - 0.5 * sin6 * sin6)); // Hy, thetax1
    betay(4) =  (-N2 + N4 * (0.25 * cos4 * cos4 - 0.5 * sin4 * sin4) + N6 * (0.25 * cos6 * cos6 - 0.5 * sin6 * sin6)); // Hy, thetax2
    betay(5) =  (-N3 + N5 * (0.25 * cos5 * cos5 - 0.5 * sin5 * sin5) + N4 * (0.25 * cos4 * cos4 - 0.5 * sin4 * sin4)); // Hy, thetax3
    betay(6) =  (-1.0) * betax(3); // Hy, thetay1
    betay(7) =  (-1.0) * betax(4); // Hy, thetay2
    betay(8) =  (-1.0) * betax(5); // Hy, thetay3
}



#endif // __MAST_dkt_elem_h__
