
/*
 *  Tet4.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/27/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

#ifndef __fesystem_tet_4_h__
#define __fesystem_tet_4_h__

// FESystem includes
#include "Mesh/VolumeElemBase.h"


namespace FESystem
{
	namespace Mesh
	{
        /*!
         *   Defines an 4-noded tetrahedral element. The node connectivity is defined as shown below. The element has a 
         *   local 3 dimensional coordinate system as shown below, and the node and coordinate association is as follows:
         *
         *   Node0 -> (xi, eta, phi) = (-1, -1, -1)
         *   Node1 -> (xi, eta, phi) = ( 1, -1, -1)
         *   Node2 -> (xi, eta, phi) = ( 1,  1, -1)
         *   Node3 -> (xi, eta, phi) = (-1,  1, -1)
         *   Node4 -> (xi, eta, phi) = ( 0,  0,  1)
         *
         *   The following six sides (boundaries) are defined for this element, and for each side, the corresponding surface normal 
         *   is defined as that of a QUAD4 element. Each boundary is defined such that the surface normal is always outwards. 
         *   Side0 -> 0, 2, 1  -> -phi
         *   Side1 -> 0, 3, 2
         *   Side2 -> 0, 2, 3
         *   Side3 -> 0, 1, 3
         *
         *               ^ phi 
         *               |
         *               |
         *               o  Node3
         *              .|.
         *            .  | 
         *          .    | .
         *         .     o    Node2
         *       .    .   \ .
         *     .  .        \   ------------> xi
         *    o-------------o
         *   Node0        Node1
         *
         *
         */
		class Tet4: public FESystem::Mesh::VolumeElemBase
		{
		public:
            /*!
             *   Constructor takes the mesh as input. 
             */
			Tet4(FESystem::Mesh::MeshBase& m);
			
			virtual ~Tet4();
			
		protected:
			
		};
	}
}

#endif // __fesystem_tet_4_h__
