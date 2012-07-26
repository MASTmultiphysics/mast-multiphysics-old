
/*
 *  Hex8.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/27/09.
 *  Copyright 2009 . All rights reserved.
 *
 */


#ifndef __fesystem_hex_8_h__
#define __fesystem_hex_8_h__

// FESystem includes
#include "Mesh/VolumeElemBase.h"


namespace FESystem
{
	namespace Mesh
	{
        /*!
         *   Defines an 8-noded hexahedral element. The node connectivity is defined as shown below. 
         *   The element has a local 3 dimensional coordinate system as shown below, and 
         *   the node and coordinate association is as follows:
         *
         *   Node0 -> (xi, eta, phi) = (-1, -1, -1)
         *   Node1 -> (xi, eta, phi) = ( 1, -1, -1)
         *   Node2 -> (xi, eta, phi) = ( 1,  1, -1)
         *   Node3 -> (xi, eta, phi) = (-1,  1, -1)
         *   Node4 -> (xi, eta, phi) = (-1, -1,  1)
         *   Node5 -> (xi, eta, phi) = ( 1, -1,  1)
         *   Node6 -> (xi, eta, phi) = ( 1,  1,  1)
         *   Node7 -> (xi, eta, phi) = (-1,  1,  1)
         *
         *   The following six sides (boundaries) are defined for this element, and for each side, the corresponding surface normal 
         *   is defined as that of a QUAD4 element. Each boundary is defined such that the surface normal is always outwards. 
         *   Side0 -> 0, 1, 2, 3 -> -phi
         *   Side1 -> 3, 2, 6, 7 -> eta
         *   Side2 -> 7, 6, 5, 4 -> phi
         *   Side3 -> 4, 5, 1, 0 -> -eta
         *   Side4 -> 1, 5, 6, 2 -> xi
         *   Side5 -> 0, 3, 7, 4 -> -xi
         *
         *   \verbatim
         *
         *             eta 
         *             ^     ^ phi
         *        Node7|    /    Node 6   
         *       o-----|---/---o 
         *      /|     |  /   /|
         * Node3 |       /   / |
         *    o-------------o Node 2 
         *    |  |          | -----------> xi
         *    |  | Node4    |  | Node 5
         *    |  o----------|--o 
         *    | /           | /
         *    |/            |/
         *    o-------------o
         *   Node0        Node1
         *
         *   \endverbatim
         */
		class Hex8: public FESystem::Mesh::VolumeElemBase
		{
		public:
            /*!
             *   Constructor takes the mesh as input. 
             */
			Hex8(FESystem::Mesh::MeshBase& m);
			
			virtual ~Hex8();
			
		protected:
			
		};
	}
}


#endif // __fesystem_hex_8_h__

