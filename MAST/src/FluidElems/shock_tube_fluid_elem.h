//
//  shock_tube_fluid_elem.h
//  MAST
//
//  Created by Manav Bhatia on 5/16/14.
//  Copyright (c) 2014 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_shock_tube_fluid_elem_h__
#define __MAST_shock_tube_fluid_elem_h__

// libMesh includes
#include "libmesh/libmesh_config.h"

#ifndef LIBMESH_USE_COMPLEX_NUMBERS

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


#endif

#endif // __MAST_shock_tube_fluid_elem_h__
