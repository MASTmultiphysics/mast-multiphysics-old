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

#ifndef __MAST__aerodynamic_qoi__
#define __MAST__aerodynamic_qoi__

// libMesh includes
#include "libmesh/getpot.h"
#include "libmesh/diff_qoi.h"

// MAST includes
#include "FluidElems/fluid_elem_base.h"





class AerodynamicQoI : public FluidElemBase, public libMesh::DifferentiableQoI
{
public:
    AerodynamicQoI(GetPot& infile):
    FluidElemBase(infile),
    DifferentiableQoI()
    {
        // set the dimension
        this->dim = infile("dimension", 0);
        
        // initialize the variable vector
        _fluid_vars.resize(dim+2);
        for (unsigned int i=0; i<dim+2; i++)
            _fluid_vars[i] = i;
            
        // this object calculates the lift, which are calculated on the side
        assemble_qoi_sides = true;
        assemble_qoi_internal_sides = false;
        assemble_qoi_elements = true;
        
        FluidElemBase::init_data();
    }
    
    virtual ~AerodynamicQoI(){}
    
    virtual void init_qoi( std::vector<Real>& sys_qoi);

    virtual void element_qoi_derivative (libMesh::DiffContext&, const libMesh::QoISet& qois);

    virtual void element_qoi (libMesh::DiffContext&, const libMesh::QoISet& qois);

    virtual void side_qoi_derivative(libMesh::DiffContext &context, const libMesh::QoISet & qois);
    
    virtual void side_qoi(libMesh::DiffContext &context, const libMesh::QoISet & qois);
    
    virtual libMesh::AutoPtr<libMesh::DifferentiableQoI> clone( )
    {
        AerodynamicQoI* new_qoi = new AerodynamicQoI(_infile);
        new_qoi->flight_condition = this->flight_condition;
        
        libMesh::AutoPtr<libMesh::DifferentiableQoI> my_clone(new_qoi);
        *my_clone = *this;
        
        return my_clone;
    }

private:
    
    std::vector<unsigned int> _fluid_vars;
};


#endif /* defined(__MAST__aerodynamic_qoi__) */
