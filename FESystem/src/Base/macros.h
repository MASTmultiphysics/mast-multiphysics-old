//
//  macros.h
//  FESystem
//
//  Created by Manav Bhatia on 2/13/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem_macros_h__
#define __fesystem_macros_h__


#define INSTANTIATE_CLASS_FOR_ALL_DATA_TYPES(class_name) \
template class class_name<FESystemFloat> ;\
template class class_name<FESystemDouble> ;\
template class class_name<FESystemComplexFloat> ;\
template class class_name<FESystemComplexDouble> 

#define INSTANTIATE_CLASS_FOR_ONLY_REAL_DATA_TYPES(class_name) \
template class class_name<FESystemFloat> ;\
template class class_name<FESystemDouble> 

#define INSTANTIATE_CLASS_FOR_ONLY_COMPLEX_DATA_TYPES(class_name) \
template class class_name<FESystemComplexFloat> ;\
template class class_name<FESystemComplexDouble> 



#endif // __fesystem_macros_h__
