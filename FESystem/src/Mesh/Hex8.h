
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
         *             ^     ^ -phi
         *        Node3|    /    Node 2
         *       o-----|---/---o 
         *      /|     |  /   /|
         * Node7 |       /   / |
         *    o-------------o Node 6
         *    |  |          | -----------> xi
         *    |  | Node0    |  | Node 1
         *    |  o----------|--o 
         *    | /           | /
         *    |/            |/
         *    o-------------o
         *   Node4        Node5
         *
         *   \endverbatim
         */
		class Hex8: public FESystem::Mesh::HexElemBase
		{
		public:
            /*!
             *   Constructor takes the mesh as input.
             */
			Hex8(FESystemBoolean local_cs_same_as_global);
			
			virtual ~Hex8();
			
            /*!
             *   Returns the geometry order for this element.
             */
            virtual FESystemUInt getGeometryOrder() const;
            
            /*!
             *   Returns a constant reference to the table that identifies the ID of the node along the local dimension
             *   This is useful for Lagrange function mapping for this element. The table is 2-D, with number of rows
             *   as the number of nodes, and the number of columns as the number of local dimensions.
             */
            virtual const FESystem::Utility::Table<FESystemUInt>& getPointDimensionIDTable() const;
            
            /*!
             *   Returns the local coordinate for \p id^th node in the element.
             */
            void getLocalComputationalCoordinateForLocalNodeID(FESystemUInt id, FESystem::Geometry::Point& p) const;
            
            /*!
             *   This defines the number of nodes along the \p n^th dimension of the element. This
             *   is useful for tasks like tensorial products of unidimensional shape functions to create
             *   shape functions in multiple dimensions.
             */
            virtual FESystemUInt getNNodesAlongDimension(FESystemUInt n) const;
            
            /*!
             *   This method supports the tensor prouct for creation of shape functions. It returns the id of the node along the
             *   dimension \p d
             */
            virtual FESystemUInt getLocalNodeIDALongDim(FESystemUInt n_id, FESystemUInt d) const;
            
            /*!
             *   Returns true if the element is a degerate element from Quad or Hex element. This is used for
             *   elements like Triangle, Tetrahedron, Prisms, Pyramids, etc.
             */
            virtual bool ifDegerateElement() const;
            
            /*!
             *   Returns the parent nondegerate element from which this element is derived.
             */
            virtual const FESystem::Mesh::ElemBase& getParentNondegenerateElem() const;
            
            /*!
             *   Returns the matrix that maps the shape functions from the parent nondegenrate to this element
             */
            virtual const FESystem::Numerics::MatrixBase<FESystemDouble>& getParentToDegenerateElemMappingMatrix() const;
            
            
		protected:
            
            /*!
             *   Table that stores the local node ID along the direction of the computational coordinates
             */
            static std::auto_ptr<FESystem::Utility::Table<FESystemUInt> > hex8_node_dim_ID_table;
            
            /*!
             *   Matrix that stores the mapping from the nondegenerate to degenerate element. This is a unit matrix for the QUAD elements
             */
            static std::auto_ptr<FESystem::Numerics::MatrixBase<FESystemDouble> > hex8_nondegenerate_to_degenerate_element_mapping;
            
            /*!
             *   Dummy function call for Quad. This throws an error, since it does not exist for Quad elements
             */
            virtual void initializeParentNondegenerateElement();
            
            /*!
             *    clears the parent nondegenerate element before it can be updated
             */
            virtual void clearParentNondegenerateElement();
            
		};
	}
}


#endif // __fesystem_hex_8_h__

