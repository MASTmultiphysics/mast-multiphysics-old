//
//  AutoPtrTable.h
//  FESystem
//
//  Created by Manav Bhatia on 4/6/12.
//  Copyright (c) 2012. All rights reserved.
//

//
//  Table.h
//  FESystem
//
//  Created by Manav Bhatia on 3/22/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem_autoptr_table_h__
#define __fesystem_autoptr_table_h__

// FESystem includes
#include "Utils/Table.h"


namespace FESystem
{
    namespace Utility
    {
        /*! 
         *   The AutoPtrTable class builds on the Table to handle a table/array of pointers. The pointers are non-shared, and will be deleted 
         *   once the AutoPtrTable object goes out of scope
         */
        template <typename ValType>
        class AutoPtrTable: public FESystem::Utility::Table<ValType*>
        {
        public: 
            /*!
             *   Default constructor
             */
            AutoPtrTable();
            
            virtual ~AutoPtrTable();

            /*!
             *   clears the data structures
             */
            virtual void clear();
            
            /*!
             *   Initializes the data container and sets all pointers to NULL
             */
            void reinit(const std::vector<FESystemUInt>& size);
            
            /*!
             *   sets the value of the element defined by \p el to \p val. It is an error to set the value if it has already assigned to a valid pointer. 
             *   In that case, use the resetVal method. 
             */
            void setVal(const std::vector<FESystemUInt>& el, ValType* val);

            /*!
             *   Sets the value of the element defined by \p el to \p val. Deletes the pointer if it already exists. 
             */
            void resetVal(const std::vector<FESystemUInt>& el, ValType* val);

            /*!
             *   Releases and returns the pointer at \p el. 
             */
            ValType* releaseVal(const std::vector<FESystemUInt>& el);

            /*!
             *   Returns the pointer at \p el. It thrown an exception if the pointer is NULL. 
             */
            ValType* getVal(const std::vector<FESystemUInt>& el);
            
            /*!
             *   Returns a constant pointer at \p el. It thrown an exception if the pointer is NULL. 
             */
            const ValType* getVal(const std::vector<FESystemUInt>& el) const;
            
            /*!
             *   Returns \p true if the value at \p el is a valid pointer. 
             */
            FESystemBoolean ifValExists(const std::vector<FESystemUInt>& el) const;
            
        protected:
            
        };
        
    }
}


template <typename ValType>
FESystem::Utility::AutoPtrTable<ValType>::AutoPtrTable():
FESystem::Utility::Table<ValType*>()
{

}

template <typename ValType>
FESystem::Utility::AutoPtrTable<ValType>::~AutoPtrTable()
{
    this->clear();
}



template <typename ValType>
void
FESystem::Utility::AutoPtrTable<ValType>::clear()
{
    // iterate over the pointers and delete them if they are not NULL
    for (FESystemUInt i=0; i<this->table_data.size(); i++)
        if (this->table_data[i] != NULL)
            delete this->table_data[i];
    
    FESystem::Utility::Table<ValType*>::clear();
}


template <typename ValType>
void
FESystem::Utility::AutoPtrTable<ValType>::reinit(const std::vector<FESystemUInt>& size)
{
    FESystem::Utility::Table<ValType*>::reinit(size);
    
    // iterate over the pointers and set them to NULL
    for (FESystemUInt i=0; i<this->table_data.size(); i++)
        this->table_data[i] = NULL;
}


template <typename ValType>
void
FESystem::Utility::AutoPtrTable<ValType>::setVal(const std::vector<FESystemUInt>& el, ValType* v)
{
    FESystemUInt index = this->getIndex(el);

    // make sure that the pointer is valid, and that the current pointer is not already set to something
    FESystemAssert0(v!=NULL, FESystem::Exception::NULLQuantity);
    FESystemAssert0(this->table_data[index] == NULL, FESystem::Exception::InvalidState);
    
    this->table_data[index] = v;
}



template <typename ValType>
void
FESystem::Utility::AutoPtrTable<ValType>::resetVal(const std::vector<FESystemUInt>& el, ValType* v)
{
    FESystemUInt index = this->getIndex(el);
    
    // make sure that the pointer is valid, and that the current pointer is not already set to something
    FESystemAssert0(v!=NULL, FESystem::Exception::NULLQuantity);
    
    if (this->table_data[index] != NULL)
        delete this->table_data[index];

    this->table_data[index] = v;
}



template <typename ValType>
ValType*
FESystem::Utility::AutoPtrTable<ValType>::releaseVal(const std::vector<FESystemUInt>& el)
{
    FESystemUInt index = this->getIndex(el);
    
    // make sure that the pointer is valid
    FESystemAssert0(this->table_data[index] != NULL, FESystem::Exception::InvalidState);
    
    ValType* v = this->table_data[index];
    this->table_data[index] = NULL;
    
    return v;
}


template <typename ValType>
ValType* 
FESystem::Utility::AutoPtrTable<ValType>::getVal(const std::vector<FESystemUInt>& el)
{
    FESystemUInt index = this->getIndex(el);

    // make sure that the pointer is valid
    FESystemAssert0(this->table_data[index] != NULL, FESystem::Exception::InvalidState);

    return this->table_data[index];
}




template <typename ValType>
const ValType* 
FESystem::Utility::AutoPtrTable<ValType>::getVal(const std::vector<FESystemUInt>& el) const
{
    FESystemUInt index = this->getIndex(el);
    
    // make sure that the pointer is valid
    FESystemAssert0(this->table_data[index] != NULL, FESystem::Exception::InvalidState);
    
    return this->table_data[index];
}



template <typename ValType>
FESystemBoolean
FESystem::Utility::AutoPtrTable<ValType>::ifValExists(const std::vector<FESystemUInt>& el) const
{
    FESystemUInt index = this->getIndex(el);
    
    return (this->table_data[index] != NULL);
}




#endif // __fesystem_autoptr_table_h__
