//
//  Prism6.h
//  FESystem
//
//  Created by Manav Bhatia on 1/24/13.
//
//

#ifndef __FESystem__Prism6__
#define __FESystem__Prism6__

// FESystem includes
#include "Mesh/VolumeElemBase.h"


namespace FESystem
{
	namespace Mesh
	{
        /*!
         *   Defines an 6-noded prism element. The node connectivity is defined as shown below.
         *   The element has a local 3 dimensional coordinate system as shown below, and
         *   the node and coordinate association is as follows:
         *
         *   Node0 -> (xi, eta, phi) = (-1, -1, -1)
         *   Node1 -> (xi, eta, phi) = ( 1, -1, -1)
         *   Node2 -> (xi, eta, phi) = ( 1,  1, -1)
         *   Node3 -> (xi, eta, phi) = (-1,  1, -1)
         *   Node4 -> (xi, eta, phi) = (-1, -1,  1)
         *   Node5 -> (xi, eta, phi) = ( 1, -1,  1)
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
         *             ^     ^ -phi
         *        Node3|    /    Node 2
         *       o-----|---/---o
         *       |  .  |  /   /|
         *       |       .   / |
         *       |          o  Node 5
         *       |          | -----------> xi
         *       | Node0    |  | Node 1
         *       o----------|--o
         *          .       | /
         *              .   |/
         *                  o
         *                    Node 4
         *
         *   \endverbatim
         */
		class Prism6: public FESystem::Mesh::PrismElemBase
		{
		public:
            /*!
             *   Constructor takes the mesh as input.
             */
			Prism6(FESystemBoolean local_cs_same_as_global);
			
			virtual ~Prism6();
			
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
            static std::auto_ptr<FESystem::Numerics::MatrixBase<FESystemDouble> > prism6_nondegenerate_to_degenerate_element_mapping;
            
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

#endif /* defined(__FESystem__Prism6__) */
