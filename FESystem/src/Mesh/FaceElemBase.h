
/*
 *  FaceElemBase.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/27/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

#ifndef __fesystem_face_elem_base_h__
#define __fesystem_face_elem_base_h__

// FESystem includes
#include "Mesh/ElemBase.h"


namespace FESystem
{
	namespace Mesh
	{
        /*!
         *   This element is the base class for elements that are faces defined in a two dimensional space
         */
		class FaceElemBase: public ElemBase
		{
		public:
            /*!
             *   The constructor takes the number of nodes in this element
             */
			FaceElemBase(FESystemUInt nnodes, FESystem::Mesh::ElementType type, FESystemBoolean local_cs_same_as_global);
			
			virtual ~FaceElemBase();
			
            /*!
             *   Returns the number of dimensions of this element
             */
            virtual FESystemUInt getDimension() const;
            
		protected:
			
            /*!
             *   Initializes the local coordinate system: uses the nodes 0, 1 and 2 of the element.
             */
            virtual void initializeLocalPhysicalCoordinateSystem();
            
            /*!
             *   clears the data structure for the local coordinate system
             */
            virtual void clearLocalPhysicalCoordinateSystem();
            
        };
        
        
        /*!
         *   This defines the base element for quads.
         */
        class QuadElemBase: public FaceElemBase
        {
        public:
            /*!
             *   The constructor takes the number of nodes in this element
             */
            QuadElemBase(FESystemUInt nnodes, FESystem::Mesh::ElementType type, FESystemBoolean local_cs_same_as_global);
            
            virtual ~QuadElemBase();
            
            /*!
             *   Each bounadry in FE is defined by one of the computational coordinates being constant. This method returns the
             *   number of that computational coordinate and its value for the boundary, in the provided \p coord_id and \p coord_val variables.
             *   This is useful for doing computations on the boundaries, for example, application of boundary conditions
             */
            virtual void getConstantCoordinateIDAndValueForBoundary(const FESystemUInt b_id, FESystemUInt& coord_id, FESystemDouble& coord_val) const;
            
            /*!
             *   Returns the map of boundary id and nodes defining the boundary
             */
            virtual const std::map<FESystemUInt, std::vector<FESystemUInt> >& getBoundaryIDAndBoundaryNodeMap() const;
            
        protected:
            
            /*!
             *   Map of boundary to set of nodes on the boundary
             */
            static std::auto_ptr<std::map<FESystemUInt, std::vector<FESystemUInt> > > quad_boundary_node_set;
            
        };
        
        
        /*!
         *   This defines the base element for triangles.
         */
        class TriElemBase: public FaceElemBase
        {
        public:
            /*!
             *   The constructor takes the number of nodes in this element
             */
            TriElemBase(FESystemUInt nnodes, FESystem::Mesh::ElementType type, FESystemBoolean local_cs_same_as_global);
            
            virtual ~TriElemBase();
            
            /*!
             *   Dummy function call for Tri. This throws an error, since it does not exist for Tri elements
             */
            void getLocalComputationalCoordinateForLocalNodeID(FESystemUInt id, FESystem::Geometry::Point& p) const;
            
            /*!
             *   Dummy function call for Tri. This throws an error, since it does not exist for Tri elements
             */
            virtual FESystemUInt getNNodesAlongDimension(FESystemUInt n) const;
            
            /*!
             *   Dummy function call for Tri. This throws an error, since it does not exist for Tri elements
             */
            virtual FESystemUInt getLocalNodeIDALongDim(FESystemUInt n_id, FESystemUInt d) const;
            
            /*!
             *   Dummy function call for Tri. This throws an error, since it does not exist for Tri elements
             */
            const FESystem::Utility::Table<FESystemUInt>& getPointDimensionIDTable() const;
            
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
             *   Each bounadry in FE is defined by one of the computational coordinates being constant. This method returns the
             *   number of that computational coordinate and its value for the boundary, in the provided \p coord_id and \p coord_val variables.
             *   This is useful for doing computations on the boundaries, for example, application of boundary conditions
             */
            virtual void getConstantCoordinateIDAndValueForBoundary(const FESystemUInt b_id, FESystemUInt& coord_id, FESystemDouble& coord_val) const;
            
            /*!
             *   Returns the map of boundary id and nodes defining the boundary
             */
            virtual const std::map<FESystemUInt, std::vector<FESystemUInt> >& getBoundaryIDAndBoundaryNodeMap() const;
            
        protected:
            
            /*!
             *   Map of boundary to set of nodes on the boundary
             */
            static std::auto_ptr<std::map<FESystemUInt, std::vector<FESystemUInt> > > tri_boundary_node_set;
            
        };
        
    }
}


#endif // __fesystem_face_elem_base_h__
