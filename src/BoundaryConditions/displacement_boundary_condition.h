/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2014  Manav Bhatia
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

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
        libMesh::DirichletBoundary& dirichlet_boundary() {
            return *_dirichlet_boundary;
        }
        
    protected:
        
        /*!
         *    Dirichlet boundary function for this boundary
         */
        std::auto_ptr<libMesh::DirichletBoundary> _dirichlet_boundary;
    };
}



inline
void
MAST::DisplacementDirichletBoundaryCondition::init(const libMesh::boundary_id_type bid,
                                                   const std::vector<unsigned int>& constrained_comp) {
    // should not have been initialized if this is called
    libmesh_assert(_dirichlet_boundary.get() == NULL);

    libMesh::ZeroFunction<Real> zero_function;
    std::set<libMesh::boundary_id_type> bid_set; bid_set.insert(bid);
    
    _dirichlet_boundary.reset(new libMesh::DirichletBoundary(bid_set,
                                                    constrained_comp,
                                                    &zero_function));
}


#endif  // __MAST_displacement_boundary_condition_h__
