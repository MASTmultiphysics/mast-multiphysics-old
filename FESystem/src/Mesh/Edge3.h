//
//  Edge3.h
//  FESystem
//
//  Created by Manav Bhatia on 5/17/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem_edge3_h__
#define __fesystem_edge3_h__


// FESystem include
#include "Mesh/EdgeElemBase.h"
#include "Utils/Table.h"


namespace FESystem
{
    // Forward declerations
    namespace Numerics {template <typename ValType> class VectorBase;}
    
	namespace Mesh
	{
        /*!
         *   Defines a 2-noded edge element lying on a straight line. The node connectivity is defined as shown below. 
         *   For this node, xi is the local coordinate and the node to local coordinate mapping is
         *   Node 0 -> xi = -1
         *   Node 1 -> xi =  1
         *   Node 2 -> xi =  0
         *
         *          Node2
         *    o------o------o  -----> xi
         *   Node0        Node1
         *
         *
         */
		class Edge3: public FESystem::Mesh::EdgeElemBase
		{
		public:
            /*!
             *   Constructor takes the mesh as input. 
             */
			Edge3();
			
			virtual ~Edge3();
			
            /*!
             *   Returns the geometry order for this element. 
             */
            virtual FESystemUInt getGeometryOrder() const; 

            /*!
             *   This assigns a unit vector to lie in the x-y plane. This needs to be defined for this element 
             *   as the transformation matrix to map the local axis quantities to global coordinate depends on the orientation
             *   of the local x-y-z axis. This vector, in combination with the direction of the element defines the x-y plane. 
             */
            void setVectorForXYPlane(const FESystem::Numerics::VectorBase<FESystemDouble>& vec);
            
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
             *   
             */
            virtual FESystemUInt getLocalNodeIDALongDim(FESystemUInt n_id, FESystemUInt d) const;
            
            /*!
             *   Returns the matrix that maps the shape functions from the parent nondegenrate to this 
             *   element
             */
            virtual const FESystem::Numerics::MatrixBase<FESystemDouble>& getParentToDegenerateElemMappingMatrix() const;
            
            /*!
             *   Table that stores the local node ID along the direction of the computational coordinates
             */
            static std::auto_ptr<FESystem::Utility::Table<FESystemUInt> > edge3_node_dim_ID_table;
            
            /*!
             *   Matrix that stores the mapping from the nondegenerate to degenerate element. This is a unit matrix for the EDGE2 element
             */
            static std::auto_ptr<FESystem::Numerics::MatrixBase<FESystemDouble> > edge3_nondegenerate_to_degenerate_element_mapping;
            
        protected:
			
		};
        
	}
}



#endif  // __fesystem_edge3_h__
