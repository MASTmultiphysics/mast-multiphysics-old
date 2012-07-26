
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
         *   This element is the base class for elements that are volumes defined in a three dimensional space
         */         
		class VolumeElemBase: public ElemBase
		{
		public:
            /*!
             *   The constructor takes the mesh to which this element belongs and the number of nodes in this element
             */
			VolumeElemBase(FESystem::Mesh::MeshBase& m, FESystemUInt nnodes);
			
			virtual ~VolumeElemBase();
			
            /*!
             *   Returns the number of dimensions of this element
             */
            virtual FESystemUInt getDimension() const; 

		protected:
			
            /*!
             *   Initializes the local coordinate system: uses the nodes 0, 1 and 2 of the element. 
             */
            virtual void initializeLocalPhysicalCoordinateSystem();
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
			HexElemBase(FESystemUInt nnodes, FESystem::Mesh::ElementType type);
			
			virtual ~HexElemBase();
            
            /*!
             *   Returns true if the element is a degerate element from Quad or Hex element. This is used for 
             *   elements like Triangle, Tetrahedron, Prisms, Pyramids, etc. 
             */
            virtual bool ifDegerateElement() const;
            
            /*!
             *   Returns the parent nondegerate element from which this element is derived. 
             */
            virtual const FESystem::Mesh::ElemBase& getParentNondegenerateElem();
            
		protected:
			
		};
        
        
        /*!
         *   This defines the base element for triangles.
         */         
		class VolumeDegenerateElemBase: public VolumeElemBase
		{
		public:
            /*!
             *   The constructor takes the number of nodes in this element
             */
			VolumeDegenerateElemBase(FESystemUInt nnodes, FESystem::Mesh::ElementType type);
			
			virtual ~VolumeDegenerateElemBase();
			
            /*!
             *   Dummy function call for Tri3. This throws an error, since it does not exist for Tri3 elements
             */
            void getLocalComputationalCoordinateForLocalNodeID(FESystemUInt id, FESystem::Geometry::Point& p) const;
            
            /*! 
             *   Dummy function call for Tri3. This throws an error, since it does not exist for Tri3 elements
             */
            virtual FESystemUInt getNNodesAlongDimension(FESystemUInt n) const;            
            
            /*!
             *   Dummy function call for Tri3. This throws an error, since it does not exist for Tri3 elements
             */
            virtual FESystemUInt getLocalNodeIDALongDim(FESystemUInt n_id, FESystemUInt d) const;
            
            /*!
             *   Dummy function call for Tri3. This throws an error, since it does not exist for Tri3 elements
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
            virtual const FESystem::Mesh::ElemBase& getParentNondegenerateElem();
            
		protected:
			
            /*!
             *   Initialize parent nondegenerate element. Is defined for each inherited element
             */
            virtual void initializeParentNondegenerateElement()=0;
            
            
            /*!
             *  parent nondegenrate element
             */
            FESystem::Mesh::ElemBase* parent_nondegenerate_elem;
		};
    
    }
}

#endif // __fesystem_volume_elem_base_h__
