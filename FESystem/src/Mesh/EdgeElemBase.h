
/*
 *  EdgeElemBase.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/27/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

#ifndef __fesystem_edge_elem_base_h__
#define __fesystem_edge_elem_base_h__

// FESystem includes
#include "Mesh/ElemBase.h"


namespace FESystem
{
	namespace Mesh
	{
        /*!
         *   Defines a 1-dimensional element
         */
		class EdgeElemBase : public ElemBase
		{
		public:
            /*!
             *   Constructor takes the mesh and number of nodes as input. 
             */
			EdgeElemBase(FESystemUInt n_nodes, FESystem::Mesh::ElementType type, FESystemBoolean local_cs_same_as_global);
			
			virtual ~EdgeElemBase();

            /*!
             *   This assigns a unit vector to lie in the x-y plane. This needs to be defined for this element 
             *   as the transformation matrix to map the local axis quantities to global coordinate depends on the orientation
             *   of the local x-y-z axis. This vector, in combination with the direction of the element defines the x-y plane. 
             */
            void setVectorForXYPlane(const FESystem::Numerics::VectorBase<FESystemDouble>& vec);

            /*!
             *   Returns the number of dimensions of this element
             */
            virtual FESystemUInt getDimension() const; 

            /*!
             *   Returns the number of boundaries for this element. A boundary is defined as a side of dimension dim-1 (eg. face for volume and line for face, etc.)
             */
            virtual FESystemUInt getNBoundaries() const;

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
             *    Calculates the normal for this element.
             */
            virtual void calculateSurfaceNormal(FESystem::Numerics::VectorBase<FESystemDouble>& n_vec) const;

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
            static std::auto_ptr<std::map<FESystemUInt, std::vector<FESystemUInt> > > edge_boundary_node_set;

            /*!
             *   Initializes the local coordinate system
             */
            virtual void initializeLocalPhysicalCoordinateSystem();
            
            /*!
             *   clears the local coordinate system
             */
            virtual void clearLocalPhysicalCoordinateSystem();
            
            /*!
             *   Dummy call, nothing to be done here
             */
            virtual void initializeParentNondegenerateElement();

            /*!
             *    clears the parent nondegenerate element before it can be updated
             */
            virtual void clearParentNondegenerateElement();

            /*!
             *   Initializes the local coordinate system
             */
            FESystem::Numerics::VectorBase<FESystemDouble>* y_unit_vec;
};
	}
}


#endif // __fesystem_edge_elem_base_h__
