
/*
 *  ParameterBase.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/27/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

#ifndef __fesystem_parameter_base_h__
#define __fesystem_parameter_base_h__


namespace FESystem
{
	namespace Base
	{
		template <typename DataType>
		class ParameterBound
		{
		public:
			ParameterBound();
			
			virtual ~ParameterBound();
			
		protected:
			
			DataType lower_bound;
			
			DataType upper_bound;			
		};
			
		
		
		template <typename DataType>
		class ParameterBase
		{
		public:
			ParameterBase();
			
			virtual ~ParameterBase();
			
		protected:
						
			FESystem::Base::ParameterBound<DataType> bounds;
			
		};
		
		
//		template < >
//		class ParameterBase<FESystemComplexDouble>
//		{
//		public:
//			
//			virtual ~ParameterBase();
//			
//		protected:
//			ParameterBase();
//			
//			FESystem::Base::ParameterBound<FESystemComplexDouble> bounds;
//			
//		};
//
//		
//		template < >
//		class ParameterBase<FESystemComplexFloat>
//		{
//		public:
//			
//			virtual ~ParameterBase();
//			
//		protected:
//			ParameterBase();
//			
//			FESystem::Base::ParameterBound<FESystemComplexFloat> bounds;
//			
//		};
		
		
	}
}


#endif // __fesystem_parameter_base_h__
