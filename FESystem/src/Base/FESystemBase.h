
/*
 *  FESystemBase.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/4/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

#ifndef __fesystem_fesystem_base_h__
#define __fesystem_fesystem_base_h__


// C++ includes
#include <map>
#include <vector>

// FESystem includes
#include "Base/FESystemExceptions.h"

// Forward declerations
class LibMeshInit;


namespace FESystem
{
	// Forward declerations
	namespace Parallel { class GlobalCommunicator;}
	
	namespace Utility { class CommandLineArguments;}
    
    namespace Geometry  {class CoordinateSystemBase;}
		
	namespace Base 
	{
        /*! 
         *   Enumerations for all library names that this code may on.
         */
		enum LibraryEnumerations 
		{
			FESystem = 0,
			MPI,
			libMesh,
			Petsc,
			Slepc,
			HDF5,
			Trilinos
		};
		
        /*!
         *   Base class for initializaiton of libraries that this code depends on.
         */
		class InitializerBase
		{
		public:
			InitializerBase();
			
			virtual ~InitializerBase();
			
			virtual void initialize() = 0;
			
			virtual void close() = 0;
			
			void setCommandArguments(const FESystem::Utility::CommandLineArguments& args);
			
			const FESystem::Utility::CommandLineArguments& getCommandArguments() const;
			
			FESystemBoolean isInitialized() const;
			
			FESystemBoolean dependencyLibraryCreated(const FESystem::Base::LibraryEnumerations& lib) const;
			
			FESystemBoolean dependencyLibraryIsIntialized(const FESystem::Base::LibraryEnumerations& lib) const;
			
			FESystem::Base::InitializerBase& getDependencyLibrary
			(const FESystem::Base::LibraryEnumerations& lib) const;

		protected:
			
			FESystemBoolean addAndInitialize(const FESystem::Base::LibraryEnumerations& lib);
			
			FESystemBoolean initialized;
			
			const FESystem::Utility::CommandLineArguments* command_line_arguments;

			typedef std::map<FESystem::Base::LibraryEnumerations, FESystem::Base::InitializerBase*>
			InitializerBaseMap;
			
			FESystem::Base::InitializerBase::InitializerBaseMap dependent_initializers_map;
			
			typedef std::vector<FESystem::Base::LibraryEnumerations> InitializerBaseVector;

			FESystem::Base::InitializerBase::InitializerBaseVector  initialization_sequence;
		};
		
		        
        /*!
         *   This is the class that initializes the FESystem library, including the data for
         */
		class FESystemLibrary: public FESystem::Base::InitializerBase
		{
		public:
			FESystemLibrary();
			
			virtual ~FESystemLibrary(); 
			
			virtual void initialize(); 
			
			virtual void close(); 
            			
		protected: 
			
		};


		
		class MPILibrary: public FESystem::Base::InitializerBase
		{
		public: 
			MPILibrary();
			
			virtual ~MPILibrary();
			
			virtual void initialize();
			
			virtual void close();
			
		protected:
			
			//FESystem::Parallel::GlobalCommunicator comm_world;
		};

	
		class libMeshLibrary: public FESystem::Base::InitializerBase
		{
		public:
			libMeshLibrary();
			
			virtual ~libMeshLibrary(); 
			
			virtual void initialize(); 
			
			virtual void close(); 
			
		protected: 
			
			std::auto_ptr<LibMeshInit> libmesh_init;
		};


		class PetscLibrary: public FESystem::Base::InitializerBase
		{
		public:
			PetscLibrary();
			
			virtual ~PetscLibrary(); 
			
			virtual void initialize(); 
			
			virtual void close(); 
			
		protected: 
			
		};

		
		class SlepcLibrary: public FESystem::Base::InitializerBase
		{
		public:
			SlepcLibrary();
			
			virtual ~SlepcLibrary(); 
			
			virtual void initialize(); 
			
			virtual void close(); 
			
		protected: 
			
		};

		
		class HDF5Library: public FESystem::Base::InitializerBase
		{
		public:
			HDF5Library();
			
			virtual ~HDF5Library(); 
			
			virtual void initialize(); 
			
			virtual void close(); 
			
		protected: 
			
		};
		
		namespace Exception 
		{

			DeclareException0(InitializedBeforeCommandLineArgument , 
							  << "Cannot Set Command Line Argument After Initialization");
			
			DeclareException1(InitializerNotCreated, std::string ,
                              << "Initializer Not Created : " << Arg1);

			DeclareException1(DependencyLibraryAlreadyExists, std::string , 
							  << "Dependency Library Already Exists in Map : " << Arg1);
			
			DeclareException1(InitializerInsertionError, std::string , 
							  << "Initializer Could Not Be Inserted In Map : " << Arg1);

			DeclareException1(LibraryDidNotInitialize, std::string , 
							  << "Library Could Not Be Initialized : " << Arg1);

			DeclareException1(LibraryDidNotClose, std::string , 
							  << "Library Could Not Be Closed : " << Arg1);
		}
		
		
		std::auto_ptr<FESystem::Base::InitializerBase>
		createInitializer(const FESystem::Base::LibraryEnumerations& lib);
	
	}
}



#endif // __fesystem_fesystem_base_h__
