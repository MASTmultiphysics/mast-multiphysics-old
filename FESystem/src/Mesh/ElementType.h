//
//  ElementTypes.h
//  FESystem
//
//  Created by Manav Bhatia on 3/19/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem_element_type_h__
#define __fesystem_element_type_h__


namespace FESystem
{
    namespace Mesh
    {
        enum ElementType
        {
            // first order elementsO
            EDGE2,
            TRI3,
            QUAD4,
            TET4,
            HEX8,
            PRISM6,
            PYRAMID5,
            // second order elements
            EDGE3,
            EDGE5,
            TRI6,
            TRI7,
            QUAD8,
            QUAD9,
            TET18,
            HEX27,
            PRISM21
        };
    }
}


#endif //__fesystem_element_type_h__
