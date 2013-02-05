
/*
 *  VolumeElemBase.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/27/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

#ifndef __fesystem_volume_elem_base_h__
#define __fesystem_volume_elem_base_h__

// FESystem includes
#include "Mesh/ElemBase.h"


namespace FESystem
{
	namespace Mesh
	{
        /*!
         *   This element is the base class for elements that are faces defined in a two dimensional space
         */
		class VolumeElemBase: public ElemBase
		{
		public:
            /*!
             *   The constructor takes the number of nodes in this element
             */
			VolumeElemBase(FESystemUInt nnodes, FESystem::Mesh::ElementType type, FESystemBoolean local_cs_same_as_global);
			
			virtual ~VolumeElemBase();
			
            /*!
             *   Returns the number of dimensions of this element
             */
            virtual FESystemUInt getDimension() const;
            
            /*!
             *    This is implementation of the pure function in the class ElemBase, which throws an exception for a volumen element.
             */
            virtual void calculateSurfaceNormal(FESystem::Numerics::VectorBase<FESystemDouble>& n_vec) const;

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
        class HexElemBase: public VolumeElemBase
        {
        public:
            /*!
             *   The constructor takes the number of nodes in this element
             */
            HexElemBase(FESystemUInt nnodes, FESystem::Mesh::ElementType type, FESystemBoolean local_cs_same_as_global);
            
            virtual ~HexElemBase();
            
            /*!
             *   Returns the number of boundaries for this element. A boundary is defined as a side of dimension dim-1 (eg. face for volume and line for face, etc.)
             */
            virtual FESystemUInt getNBoundaries() const;

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
            static std::auto_ptr<std::map<FESystemUInt, std::vector<FESystemUInt> > > hex_boundary_node_set;
            
        };
        
        
        /*!
         *   This defines the base element for triangles.
         */
        class PrismElemBase: public VolumeElemBase
        {
        public:
            /*!
             *   The constructor takes the number of nodes in this element
             */
            PrismElemBase(FESystemUInt nnodes, FESystem::Mesh::ElementType type, FESystemBoolean local_cs_same_as_global);
            
            virtual ~PrismElemBase();
            
            /*!
             *   Returns the number of boundaries for this element. A boundary is defined as a side of dimension dim-1 (eg. face for volume and line for face, etc.)
             */
            virtual FESystemUInt getNBoundaries() const;

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
             *   elements like Triangle, Tetrahedron, Prisms, Tets, etc.
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
            static std::auto_ptr<std::map<FESystemUInt, std::vector<FESystemUInt> > > prism_boundary_node_set;
            
        };

        
        /*!
         *   This defines the base element for triangles.
         */
        class TetElemBase: public VolumeElemBase
        {
        public:
            /*!
             *   The constructor takes the number of nodes in this element
             */
            TetElemBase(FESystemUInt nnodes, FESystem::Mesh::ElementType type, FESystemBoolean local_cs_same_as_global);
            
            virtual ~TetElemBase();
            
            /*!
             *   Returns the number of boundaries for this element. A boundary is defined as a side of dimension dim-1 (eg. face for volume and line for face, etc.)
             */
            virtual FESystemUInt getNBoundaries() const;

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
             *   elements like Triangle, Tetrahedron, Prisms, Tets, etc.
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
            static std::auto_ptr<std::map<FESystemUInt, std::vector<FESystemUInt> > > tet_boundary_node_set;
            
        };

        
        
        /*!
         *   This defines the base element for triangles.
         */
        class PyramidElemBase: public VolumeElemBase
        {
        public:
            /*!
             *   The constructor takes the number of nodes in this element
             */
            PyramidElemBase(FESystemUInt nnodes, FESystem::Mesh::ElementType type, FESystemBoolean local_cs_same_as_global);
            
            virtual ~PyramidElemBase();
            
            /*!
             *   Returns the number of boundaries for this element. A boundary is defined as a side of dimension dim-1 (eg. face for volume and line for face, etc.)
             */
            virtual FESystemUInt getNBoundaries() const;

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
             *   elements like Triangle, Tetrahedron, Prisms, Tets, etc.
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
            static std::auto_ptr<std::map<FESystemUInt, std::vector<FESystemUInt> > > pyramid_boundary_node_set;
            
        };

    }
}

#endif // __fesystem_volume_elem_base_h__
