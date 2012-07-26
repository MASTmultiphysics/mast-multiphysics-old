
/*
 *  Pyramid5.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/27/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

#ifndef __fesystem_pyramid_5_h__
#define __fesystem_pyramid_5_h__

// FESystem includes
#include "Mesh/VolumeElemBase.h"

namespace FESystem
{
	namespace Mesh
	{
        /*!
         *   Defines an 5-noded prism element. The node connectivity is defined as shown below. The element has a 
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
         *   Side0 -> 0, 3, 2, 1 -> -phi
         *   Side1 -> 0, 4, 3    
         *   Side2 -> 1, 2, 4
         *   Side3 -> 3, 2, 4
         *
         *              ^ phi 
         *              |
         *              |
         *             Node4
         *              o  
         *             /.\   ^ eta
         *            /. .\ / 
         *           /.   ./
         *          /.    / \
         *         /.    / . \
         *        /.    /   . \
         * Node3 o-------------o Node2 
         *      /.          . /
         *     /.            /------------> xi
         *    o-------------o
         *   Node0        Node1
         *
         *
         */
		class Prism5: public FESystem::Mesh::VolumeElemBase
		{
		public:
            /*!
             *   Constructor takes the mesh as input. 
             */
			Prism5(FESystem::Mesh::MeshBase& m);
			
			virtual ~Prism5();
			
		protected:
			
		};
	}
}

#endif // __fesystem_pyramid_5_h__
