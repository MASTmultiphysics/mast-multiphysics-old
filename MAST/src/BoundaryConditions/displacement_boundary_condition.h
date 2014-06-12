//
//  displacement_boundary_condition.h
//  MAST
//
//  Created by Manav Bhatia on 12/16/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_displacement_boundary_condition_h__
#define __MAST_displacement_boundary_condition_h__

// C++ includes
#include <memory>

// MAST includes
#include "BoundaryConditions/boundary_condition.h"

// libmesh includes
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/zero_function.h"


namespace MAST {
    
    class DisplacementDirichletBoundaryCondition: public MAST::BoundaryCondition {
      
    public:
        DisplacementDirichletBoundaryCondition():
        MAST::BoundaryCondition(MAST::DISPLACEMENT_DIRICHLET)
        { }
        
        virtual ~DisplacementDirichletBoundaryCondition() { }

        /*!
         *   initializes the object for the specified domain id (either boundary, 
         *   or subdomain), for the displacement components initialized using 
         *   a bitwise operator. This method initializes the components to zero
         *   value on the domain.
         */
        virtual void init(const libMesh::boundary_id_type bid,
                          const std::vector<unsigned int>& constrained_comp);
        
        
        /*!
         *    Returns a reference to the Dirichlet boundary condition object
         */
        DirichletBoundary& dirichlet_boundary() {
            return *_dirichlet_boundary;
        }
        
    protected:
        
        /*!
         *    Dirichlet boundary function for this boundary
         */
        std::auto_ptr<DirichletBoundary> _dirichlet_boundary;
    };
}



inline
void
MAST::DisplacementDirichletBoundaryCondition::init(const libMesh::boundary_id_type bid,
                                                   const std::vector<unsigned int>& constrained_comp) {
    // should not have been initialized if this is called
    libmesh_assert(_dirichlet_boundary.get() == NULL);

    ZeroFunction<libMesh::Real> zero_function;
    std::set<libMesh::boundary_id_type> bid_set; bid_set.insert(bid);
    
    _dirichlet_boundary.reset(new DirichletBoundary(bid_set,
                                                    constrained_comp,
                                                    &zero_function));
}


#endif  // __MAST_displacement_boundary_condition_h__
