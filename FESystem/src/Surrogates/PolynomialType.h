/*
 *  PolynomialType.h
 *  FESystem
 *
 *  Created by Manav  Bhatia on 2/16/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __fesystem_polynomial_type_h__
#define __fesystem_polynomial_type_h__

namespace FESystem
{
  namespace Surrogates
  {
    enum PolynomialType
    {
      ZEROTH_ORDER,
      FIRST_ORDER,
      SECOND_ORDER_COMPLETE,
      SECOND_ORDER_NO_COUPLING,
      THIRD_ORDER_COMPLETE,
      THIRD_ORDER_NO_COUPLING
    };
  }
}

#endif // --fesystem_polynomial_type_h__
