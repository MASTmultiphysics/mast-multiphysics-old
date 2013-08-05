//
//  aerodynamic_qoi.h
//  MAST
//
//  Created by Manav Bhatia on 6/11/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST__aerodynamic_qoi__
#define __MAST__aerodynamic_qoi__

// libMesh includes
#include "libmesh/diff_qoi.h"

// FESystem includes
#include "FluidElems/fluid_elem_base.h"


// Bring in everything from the libMesh namespace
using namespace libMesh;

class AerodynamicQoI : public EulerElemBase, public DifferentiableQoI
{
public:
    AerodynamicQoI(const unsigned int dimension):
    DifferentiableQoI()
    {
        // set the dimension
        this->dim = dimension;
        
        // initialize the variable vector
        _fluid_vars.resize(dim+2);
        for (unsigned int i=0; i<dim+2; i++)
            _fluid_vars[i] = i;
            
        // this object calculates the lift, which are calculated on the side
        assemble_qoi_sides = true;
        assemble_qoi_internal_sides = false;
        assemble_qoi_elements = false;
        
        EulerElemBase::init_data();
    }
    
    virtual ~AerodynamicQoI(){}
    
    virtual void init_qoi( std::vector<Number>& sys_qoi);

    virtual void element_qoi_derivative (DiffContext&, const QoISet& qois);

    virtual void element_qoi (DiffContext&, const QoISet& qois);

    virtual void side_qoi_derivative(DiffContext &context, const QoISet & qois);
    
    virtual void side_qoi(DiffContext &context, const QoISet & qois);
    
    virtual AutoPtr<DifferentiableQoI> clone( )
    {
        AerodynamicQoI* new_qoi = new AerodynamicQoI(dim);
        new_qoi->flight_condition = this->flight_condition;
        
        AutoPtr<DifferentiableQoI> my_clone(new_qoi);
        *my_clone = *this;
        
        return my_clone;
    }

private:
    
    std::vector<unsigned int> _fluid_vars;
};


#endif /* defined(__MAST__aerodynamic_qoi__) */
