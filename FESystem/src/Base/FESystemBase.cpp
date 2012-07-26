
/*
 *  FESystemBase.cpp
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/4/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

// C++ includes
#include <memory>

// FESystem includes
#include "Base/FESystemBase.h"
#include "Utils/EnumName.h" 
#include "Utils/CommandLineArguments.h"
#include "Geom/RectangularCoordinateSystem.h"

// MPI includes
//#include "mpi.h"

// libmesh includes
//#include "base/libmesh.h"

// Petsc includes
//#include "petsc.h"
//
// Slepc includes
//#include "slepc.h"

//// HDF5 includes
//#include "hdf5.h"


FESystem::Base::InitializerBase::InitializerBase():
initialized(false),
command_line_arguments(NULL)
{
	
}


FESystem::Base::InitializerBase::~InitializerBase()
{
	// by default, iterate over the initialization sequence in reverse,  
	// close them, and delete the pointers
	FESystem::Base::InitializerBase::InitializerBaseVector::reverse_iterator it, end;
	it = this->initialization_sequence.rbegin();
	end = this->initialization_sequence.rend();
	for ( ; it != end; it++) 
	{
		FESystem::Base::InitializerBase& init_base = this->getDependencyLibrary(*it);
		if (init_base.isInitialized())
			init_base.close();
		delete &init_base;
	}
}


FESystemBoolean 
FESystem::Base::InitializerBase::isInitialized() const
{
	return this->initialized;
}


void 
FESystem::Base::InitializerBase::setCommandArguments
(const FESystem::Utility::CommandLineArguments& args)
{
	FESystemAssert0(!this->isInitialized(),
                  FESystem::Base::Exception::InitializedBeforeCommandLineArgument);
	
	this->command_line_arguments = &args;
}


const FESystem::Utility::CommandLineArguments& 
FESystem::Base::InitializerBase::getCommandArguments() const
{
	FESystemAssert0(this->command_line_arguments != NULL,
                  FESystem::Exception::NULLQuantity);
	
	return *(this->command_line_arguments);	
}


FESystemBoolean 
FESystem::Base::InitializerBase::dependencyLibraryCreated
(const FESystem::Base::LibraryEnumerations& lib) const
{
	// check if the library exists in the map
	return this->dependent_initializers_map.count(lib);
}



FESystemBoolean
FESystem::Base::InitializerBase::dependencyLibraryIsIntialized
(const FESystem::Base::LibraryEnumerations& lib) const
{
	// if the library does not exist in the map, return false, else
	// check if the library exists
	FESystemBoolean return_val = false;
	FESystem::Base::InitializerBase::InitializerBaseMap::const_iterator it, end;
	it = this->dependent_initializers_map.find(lib);
	end = this->dependent_initializers_map.end();
	
	// if something was found, then check for initialization
	if (it != end)
		return_val = it->second->isInitialized();
	
	return return_val;
}



FESystem::Base::InitializerBase& 
FESystem::Base::InitializerBase::getDependencyLibrary
(const FESystem::Base::LibraryEnumerations& lib) const
{
	// if the library does not exist in the map, return false, else
	// check if the library exists
	FESystem::Base::InitializerBase::InitializerBaseMap::const_iterator it, end;
	it = this->dependent_initializers_map.find(lib);
	end = this->dependent_initializers_map.end();
  
	// check for existence
	FESystemAssert1(it != end, 
                  FESystem::Base::Exception::InitializerNotCreated,
                  FESystem::Utility::libraryEnumToName(lib));
	
	// if something was found, then return this library
	return *(it->second);
}




FESystemBoolean
FESystem::Base::InitializerBase::addAndInitialize(const FESystem::Base::LibraryEnumerations& lib)
{
	// first, make sure that this library isn't already initialized. If it is not, then
	// create an instance of the library, and then initialize it
	FESystemAssert1(!(this->dependencyLibraryCreated(lib)),
                  FESystem::Base::Exception::DependencyLibraryAlreadyExists,
                  FESystem::Utility::libraryEnumToName(lib));
	
	// create the library
	std::auto_ptr<FESystem::Base::InitializerBase> 
	initializer(FESystem::Base::createInitializer(lib));
	
	std::pair<FESystem::Base::InitializerBase::InitializerBaseMap::iterator, bool> insert_return = 
	this->dependent_initializers_map.insert
	(FESystem::Base::InitializerBase::InitializerBaseMap::value_type(lib, initializer.release()));
	
	FESystemStrictAssert1(insert_return.second, 
                        FESystem::Base::Exception::InitializerInsertionError,
                        FESystem::Utility::libraryEnumToName(lib));
	
	// provide access of the command line arguments to the initializer
	insert_return.first->second->setCommandArguments(this->getCommandArguments());
	
	// now initialize the library
	insert_return.first->second->initialize();
	
	// record this sequence in the vector
	this->initialization_sequence.push_back(lib);
	
	// if we get this far, then all is well... keep working
	return true;
}




FESystem::Base::FESystemLibrary::FESystemLibrary():
FESystem::Base::InitializerBase()
{

}



FESystem::Base::FESystemLibrary::~FESystemLibrary()
{

}


void
FESystem::Base::FESystemLibrary::initialize()
{
	this->addAndInitialize(FESystem::Base::MPI);
  //	this->addAndInitialize(FESystem::Base::Petsc);
  //	this->addAndInitialize(FESystem::Base::Slepc);
  //	this->addAndInitialize(FESystem::Base::HDF5);
  //	this->addAndInitialize(FESystem::Base::libMesh);	
}


void
FESystem::Base::FESystemLibrary::close()
{
	
}



FESystem::Base::MPILibrary::MPILibrary():
FESystem::Base::InitializerBase()
{
	
}



FESystem::Base::MPILibrary::~MPILibrary()
{
	
}


void
FESystem::Base::MPILibrary::initialize()
{
	std::pair<FESystemUInt, char* const*> val = this->getCommandArguments().getArgumentsPair();
	FESystemInt num = val.first;
	char*** str = const_cast<char***>(&val.second);
	
//	MPI_Init(&(num), str);
}


void
FESystem::Base::MPILibrary::close()
{
//	MPI_Finalize();  
}



//FESystem::Base::libMeshLibrary::libMeshLibrary():
//FESystem::Base::InitializerBase()
//{
//	
//}
//
//
//
//FESystem::Base::libMeshLibrary::~libMeshLibrary()
//{
//	
//}


//void
//FESystem::Base::libMeshLibrary::initialize()
//{
//	std::pair<FESystemUInt, char* const*> val = this->getCommandArguments().getArgumentsPair();
//	FESystemInt num = val.first;
//	char** str = const_cast<char**>(val.second);
//	
//	this->libmesh_init.reset(new LibMeshInit(num, str));
//}


//void
//FESystem::Base::libMeshLibrary::close()
//{
//	this->libmesh_init.reset();
//}
//
//
//
//FESystem::Base::PetscLibrary::PetscLibrary():
//FESystem::Base::InitializerBase()
//{
//	
//}
//
//
//
//FESystem::Base::PetscLibrary::~PetscLibrary()
//{
//	
//}
//
//
//void
//FESystem::Base::PetscLibrary::initialize()
//{
//	std::pair<FESystemUInt, char* const*> val = this->getCommandArguments().getArgumentsPair();
//	FESystemInt num = val.first;
//	char*** str = const_cast<char***>(&val.second);
//
//	PetscErrorCode ierr = 0;
//	ierr = PetscInitialize(&num, str, PETSC_NULL, PETSC_NULL);
//	CHKERRABORT(PETSC_COMM_WORLD, ierr);	
//}
//
//
//void
//FESystem::Base::PetscLibrary::close()
//{
//	PetscErrorCode ierr = 0;
//	ierr = PetscFinalize();
//	CHKERRABORT(PETSC_COMM_WORLD, ierr);
//
//}
//
//
//
//FESystem::Base::SlepcLibrary::SlepcLibrary():
//FESystem::Base::InitializerBase()
//{
//	
//}
//
//
//
//FESystem::Base::SlepcLibrary::~SlepcLibrary()
//{
//	
//}
//
//
//void
//FESystem::Base::SlepcLibrary::initialize()
//{
//	std::pair<FESystemUInt, char* const*> val = this->getCommandArguments().getArgumentsPair();
//	FESystemInt num = val.first;
//	char*** str = const_cast<char***>(&val.second);
//
//	PetscErrorCode ierr = 0;
//	ierr = SlepcInitialize(&num, str, PETSC_NULL, PETSC_NULL);
//	CHKERRABORT(PETSC_COMM_WORLD, ierr);
//
//}
//
//
//void
//FESystem::Base::SlepcLibrary::close()
//{
//	PetscErrorCode ierr = 0;
//	ierr = SlepcFinalize();
//	CHKERRABORT(PETSC_COMM_WORLD, ierr);
//}



FESystem::Base::HDF5Library::HDF5Library():
FESystem::Base::InitializerBase()
{
	
}



FESystem::Base::HDF5Library::~HDF5Library()
{
	
}


void
FESystem::Base::HDF5Library::initialize()
{
//	FESystemUInt ierr = 0;
//	ierr = H5open();
//	FESystemStrictAssert1(ierr >= 0, 
//                        FESystem::Base::Exception::LibraryDidNotInitialize,
//                        "HDF5");
}


void
FESystem::Base::HDF5Library::close()
{
//	FESystemUInt ierr = 0;
//	ierr = H5close();
//	FESystemStrictAssert1(ierr >= 0, 
//                        FESystem::Base::Exception::LibraryDidNotClose,
//                        "HDF5");
}


std::auto_ptr<FESystem::Base::InitializerBase>
FESystem::Base::createInitializer(const FESystem::Base::LibraryEnumerations& lib)
{
	std::auto_ptr<FESystem::Base::InitializerBase> initializer;
	
	switch (lib) 
	{
		case FESystem::Base::FESystem:
			initializer.reset(new FESystem::Base::FESystemLibrary());
			break;
      
		case FESystem::Base::MPI:
			initializer.reset(new FESystem::Base::MPILibrary());
			break;
      
      //		case FESystem::Base::libMesh:
      //			initializer.reset(new FESystem::Base::libMeshLibrary());
      //			break;
      //
      //		case FESystem::Base::Petsc:
      //			initializer.reset(new FESystem::Base::PetscLibrary());
      //			break;
      //
      //		case FESystem::Base::Slepc:
      //			initializer.reset(new FESystem::Base::SlepcLibrary());
      //			break;
      
		case FESystem::Base::HDF5:
			initializer.reset(new FESystem::Base::HDF5Library());
			break;
			
		default:
			FESystemAssert1(false, FESystem::Exception::EnumerationNotHandled, lib);
			break;
	}
	
	return initializer;
}


