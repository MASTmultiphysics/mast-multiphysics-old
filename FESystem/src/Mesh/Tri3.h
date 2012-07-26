
/*
 *  Tri3.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/27/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

#ifndef __fesystem_tri_3_h__
#define __fesystem_tri_3_h__

// FESystem includes
#include "Mesh/FaceElemBase.h"


namespace FESystem
{
	namespace Mesh
	{
        /*!
         *   Defines a 3-noded triangular element lying in a plane surface. The node connectivity is defined as shown below. 
         *   The surface normal for the face is defined based on the right-hand-rule going from Node0->Node1->Node2. So,
         *   for the element shown here, the surface normal will be a vector coming out of the screen. 
         *   This is a degenerate element created using a QUAD element, where the last two nodes of the quad are coincident 
         *   with Node2. 
         *   Node 0 -> (xi, eta) = (-1,-1) 
         *   Node 1 -> (xi, eta) = ( 1,-1) 
         *   Node 2 -> (xi, eta) = ( 1, 1) = ( -1, 1)
         *
         *   \verbatim
         *
         *         ^ eta
         *         |
         *         |
         *           Node2
         *         o
         *        . .
         *       .   . 
         *      .     . --------> xi
         *     .       .
         *    .         .
         *    o----------o    
         *   Node0        Node1
         *
         *  \endverbatim
         *
         */
		class Tri3: public FESystem::Mesh::TriElemBase
		{
		public:
            /*!
             *   Constructor takes the mesh as input. 
             */
			Tri3();
			
			virtual ~Tri3();
			
            /*!
             *   Returns the geometry order for this element. 
             */
            virtual FESystemUInt getGeometryOrder() const; 

            /*!
             *   Returns the matrix that maps the shape functions from the parent nondegenrate to this element
             */
            virtual const FESystem::Numerics::MatrixBase<FESystemDouble>& getParentToDegenerateElemMappingMatrix() const;
                        
            /*!
             *   Matrix that stores the mapping from the nondegenerate to degenerate element. This is a unit matrix for the QUAD elements
             */
            static std::auto_ptr<FESystem::Numerics::MatrixBase<FESystemDouble> > tri3_nondegenerate_to_degenerate_element_mapping;

		protected:
			
            /*!
             *   Initialize parent nondegenerate element. Is defined for each inherited element
             */
            virtual void initializeParentNondegenerateElement();

		};
	}
}


#endif // __fesystem_tri_3_h__
