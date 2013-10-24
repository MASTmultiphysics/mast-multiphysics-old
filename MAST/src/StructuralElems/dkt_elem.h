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
            // first three nodes are same as that of the original element
            for (unsigned int i=0; i<3; i++)
                _tri6->set_node(i) = elem.get_node(i);
            
            // next three nodes are created using the corner nodes
            for (unsigned int i=0; i<3; i++)
                _nodes[i] = new Node;

            // first edge
            (*_nodes[0])  = *elem.get_node(0);
            (*_nodes[0]) += *elem.get_node(1);
            (*_nodes[0]) *= 0.5;
            
            // second edge
            (*_nodes[1])  = *elem.get_node(1);
            (*_nodes[1]) += *elem.get_node(2);
            (*_nodes[1]) *= 0.5;

            // third edge
            (*_nodes[2])  = *elem.get_node(2);
            (*_nodes[2]) += *elem.get_node(0);
            (*_nodes[2]) *= 0.5;

            // now setup the shape functions
            _fe = FEBase::build(2, FEType(SECOND, LAGRANGE)).release();
            _fe->attach_quadrature_rule(&qrule);
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
    libmesh_error();
}




#endif // __MAST_dkt_elem_h__
