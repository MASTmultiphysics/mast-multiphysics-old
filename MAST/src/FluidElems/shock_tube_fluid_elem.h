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

#ifndef __MAST_shock_tube_fluid_elem_h__
#define __MAST_shock_tube_fluid_elem_h__

// libMesh includes
#include "libmesh/libmesh_config.h"



// MAST includes
#include "FluidElems/fluid_system.h"

namespace MAST {
    
    class ShockTubeFluidElem: public FluidSystem {
    public:
        ShockTubeFluidElem(libMesh::EquationSystems& es,
                           const std::string& name_in,
                           const unsigned int number_in);
    
        virtual ~ShockTubeFluidElem();
        
        virtual bool element_time_derivative (bool request_jacobian,
                                              libMesh::DiffContext &context);
        
        
    protected:
        
        void calculate_source_flux(const PrimitiveSolution& sol,
                                   libMesh::DenseVector<Real>& vec);

        void calculate_source_flux_jacobian(const PrimitiveSolution& sol,
                                            libMesh::DenseMatrix<Real>& mat);
    };
}



#endif // __MAST_shock_tube_fluid_elem_h__
