
/*
 *  DataSourceBase.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 12/3/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

#ifndef __fesystem_data_source_base_h__
#define __fesystem_data_source_base_h__


namespace FESystem
{
	// Forward declerations
	namespace Base { class ParameterBase; }
	
	namespace Data
	{
		template <typename DataType> 
		class DataSourceBase
		{
		public:
			DataSourceBase();
			
			virtual ~DataSourceBase();
			
			virtual DataType& getData() = 0;

			virtual DataType& getDataSensitivity(const FESystem::Base::ParameterBase& param) = 0;

		protected:
			
			
			
		};
	}
}


#endif // __fesystem_data_source_base_h__
