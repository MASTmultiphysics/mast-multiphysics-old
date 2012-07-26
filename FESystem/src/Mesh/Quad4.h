
/*
 *  Quad4.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/27/09.
 *  Copyright 2009 . All rights reserved.
 *
 */


#ifndef __fesystem_quad_4_h__
#define __fesystem_quad_4_h__

// FESystem includes
#include "Mesh/FaceElemBase.h"
#include "Utils/Table.h"

namespace FESystem
{
	namespace Mesh
	{
        /*!
         *   Defines a 4-noded quadrilateral element lying in a plane surface. The node connectivity is defined as shown below. 
         *   The surface normal for the face is defined based on the right-hand-rule going from Node0->Node1->Node2->Node3. So,
         *   for the element shown here, the surface normal will be a vector coming out of the screen. The element ensures that 
         *   all four nodes lie in a single plane. The element has a coordinate system as shown in this sketch, where the coordinate 
         *   and node association is
         *
         *   Node 0 -> (xi, eta) = (-1, -1)
         *   Node 1 -> (xi, eta) = ( 1, -1)
         *   Node 2 -> (xi, eta) = ( 1,  1)
         *   Node 3 -> (xi, eta) = (-1,  1)
         *
         *   The following four sides (boundaries) are defined for this element, and for each side, the corresponding normal 
         *   is defined as that of a EDGE2 element with xi from first to second node, and eta obtained by cross-product of QUAD4 surface normal 
         *   with the EDGE2 xi.  
         *   Side0 -> 0, 1 -> eta
         *   Side1 -> 1, 2 -> -xi
         *   Side2 -> 2, 3 -> -eta
         *   Side3 -> 3, 0 -> xi
         *
         *   \verbatim
         *
         *
         *           ^  eta
         *           |
         *           |
         *           |
         *  Node3          Node2
         *    o-------------o
         *    |             |
         *    |             |
         *    |             | ------> xi
         *    |             |
         *    |             |
         *    o-------------o
         *   Node0        Node1
         *
         *   \endverbatim
         *
         */
		class Quad4: public FESystem::Mesh::QuadElemBase
		{
		public:
            /*!
             *   Constructor takes the mesh as input. 
             */
			Quad4();
			
			virtual ~Quad4();
			
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
            static std::auto_ptr<FESystem::Utility::Table<FESystemUInt> > quad4_node_dim_ID_table;
            
            /*!
             *   Matrix that stores the mapping from the nondegenerate to degenerate element. This is a unit matrix for the QUAD elements
             */
            static std::auto_ptr<FESystem::Numerics::MatrixBase<FESystemDouble> > quad4_nondegenerate_to_degenerate_element_mapping;

            /*!
             *   Dummy function call for Quad. This throws an error, since it does not exist for Quad elements
             */
            virtual void initializeParentNondegenerateElement();

			
		};
        
	}
}

#endif // __fesystem_quad_4_h__
