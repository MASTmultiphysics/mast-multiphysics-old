
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
         *
         *
         *   \verbatim
         *
         *             ^ eta
         *             |   ^ -phi
         *        Node3|  /      Node 2
         *       o-----|-/-----o
         *       |     |    .  |
         *       |             |
         *       |        .    | -------> xi
         *       | .           |
         *       | Node0       | Node 1
         *       o-----.-------o
         *        \ . .
         *          o   Node 4
         *
         *
         *
         *   \endverbatim
         *
         *
         */
		class Pyramid5: public FESystem::Mesh::PyramidElemBase
		{
		public:
            /*!
             *   Constructor takes the mesh as input.
             */
			Pyramid5(FESystemBoolean local_cs_same_as_global);
			
			virtual ~Pyramid5();
			
            /*!
             *   Returns the geometry order for this element.
             */
            virtual FESystemUInt getGeometryOrder() const;
            
            /*!
             *   Returns the matrix that maps the shape functions from the parent nondegenrate to this element
             */
            virtual const FESystem::Numerics::MatrixBase<FESystemDouble>& getParentToDegenerateElemMappingMatrix() const;
            
		protected:
			
            /*!
             *   Matrix that stores the mapping from the nondegenerate to degenerate element. This is a unit matrix for the QUAD elements
             */
            static std::auto_ptr<FESystem::Numerics::MatrixBase<FESystemDouble> > pyramid5_nondegenerate_to_degenerate_element_mapping;
            
            /*!
             *   Initialize parent nondegenerate element. Is defined for each inherited element
             */
            virtual void initializeParentNondegenerateElement();
            
            /*!
             *    clears the parent nondegenerate element before it can be updated
             */
            virtual void clearParentNondegenerateElement();
		};
	}
}

#endif // __fesystem_pyramid_5_h__
