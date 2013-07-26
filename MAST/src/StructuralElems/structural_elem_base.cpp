
// MAST includes
#include "StructuralElems/structural_elem_base.h"

#ifndef LIBMESH_USE_COMPLEX_NUMBERS

StructuralElementBase::StructuralElementBase():
_if_initialized(false),
_elem(NULL),
_system(NULL)
{
    
}


StructuralElementBase::~StructuralElementBase()
{
    
}



void
StructuralElementBase::clear()
{
    _if_initialized = false;
    _elem = NULL;
    _system = NULL;
}



void
StructuralElementBase::initialize(const Elem& elem, const StructuralSystemBase& system)
{
    libmesh_assert(!_if_initialized);
    
    _elem = &elem;
    _system = &system;
}

#endif // LIBMESH_USE_COMPLEX_NUMBERS

