
//
//  Table.h
//  FESystem
//
//  Created by Manav Bhatia on 3/22/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem_table_h__
#define __fesystem_table_h__

// C++ includes
#include <vector>

// FESystem includes
#include "Base/FESystemTypes.h"
#include "Base/FESystemExceptions.h"


namespace FESystem
{
    namespace Utility
    {
        /*!
         *   This class serves as a data storage in a multidimensional array/table format. The table can have any number of elements in 
         *   each of the dimensions. This is only intended to be a storage container for arbitrary data types, and not for operations (/+-*) 
         *   on the elements. 
         */
        template <typename ValType> 
        class Table
        {
        public:
            /*!
             *  This is the default constructor. Once the table is instantiated, it should be reinitialized with the required dimensions. 
             */
            Table();
            
            virtual ~Table();
            
            /*!
             *    Clears the data structure
             */
            virtual void clear();
            
            /*!
             *   This will initialize the data structures for storing the data. The dimensions of the table is the size of the vector, 
             *   while the number of elements in each dimension is the value stored in the vector. Note that it is an error to have a zero 
             *   number of elements for a dimension.
             */
            void reinit(const std::vector<FESystemUInt>& size);
            
            /*!
             *    This sets the value of the index defined by the index \p el. The size of \p el should be equal to the dimension of this
             *    table.
             */
            void setVal(const std::vector<FESystemUInt>& el, ValType& v);
            
            /*!
             *    This returns a reference to the value of the index defined by the index \p el. The size of \p el should be equal to the dimension of this
             *    table.
             */
            ValType& getVal(const std::vector<FESystemUInt>& el);
            
            
            /*!
             *    This returns a constant reference to the value of the index defined by the index \p el. The size of \p el should be equal to the
             *    dimension of this table.
             */
            const ValType& getVal(const std::vector<FESystemUInt>& el) const;
            
            /*!
             *    This returns a reference to the vector of values used for internal storage. The user should avoid using this method to access 
             *    data unless absolutely necessary, since the other access methods like getVal and setVal ensure the correct indexing. 
             */
            std::vector<ValType>& getAllValues();
            
            
            /*!
             *    This returns a constant reference to the vector of values used for internal storage. The user should avoid using this method to access 
             *    data unless absolutely necessary, since the other access methods like getVal and setVal ensure the correct indexing. 
             */
            const std::vector<ValType>& getAllValues() const;
            
            /*!
             *    Returns the dimension of the table
             */
            FESystemUInt getDimension() const;
            
            /*!
             *    Returns the maximum number of elements in the dimension \p dim.
             */
            FESystemUInt getNElementsInDimension(FESystemUInt dim) const;
            
            /*!
             *    Returns a constant reference to the vector storing the dimensions
             */
            const std::vector<FESystemUInt>& getSizeVector() const;            
            
            
        protected:            
            
            /*!
             *   Returns the index of the entry in the local table_data vector for the element defined by \p el
             */
            FESystemUInt getIndex(const std::vector<FESystemUInt>& el) const;
            
            /*!
             *   The dimensions of the table is the size of the vector, while the number of elements in each dimension is the value stored 
             *   in the vector. Note that it is an error to have a zero number of elements for a dimension. 
             */
            std::vector<FESystemUInt> table_size;
            
            
            /*!
             *   The i^th element of this vector stores the product of the number of elements in all 0 ... (i-1)st dimensions. This is required 
             *   for calculation of element location in table_data. 
             */
            std::vector<FESystemUInt> table_plane_products;            
            
            
            /*!
             *   The data is stored in this vector. The data is stored in a column major format, and this is used to convert the table entry index
             *   to the location of the element in this vector. 
             */
            
            std::vector<ValType> table_data;
        };
        
                
        /*!
         *   Exception for zero give table size
         */ 
        DeclareException0(ZeroTableSize, << "Number of Elements in Table Cannot Be Zero. One or More Dimensions Have Zero Elements.");
        
    }
}


template <typename ValType>
FESystem::Utility::Table<ValType>::Table()
{
    
}

template <typename ValType>
FESystem::Utility::Table<ValType>::~Table()
{

}


template <typename ValType>
void
FESystem::Utility::Table<ValType>::clear()
{
    this->table_size.clear();
    this->table_plane_products.clear();
    this->table_data.clear();
}


template <typename ValType>
void
FESystem::Utility::Table<ValType>::reinit(const std::vector<FESystemUInt>& size)
{
    // clear before initializing the data
    this->table_size.clear();
    this->table_plane_products.clear();
    this->table_data.clear();
    
    // copy the size to the local vector
    this->table_size = size;
    
    // now calculate the total number of element
    FESystemUInt nvals = 1;
    for (FESystemUInt i=0; i<this->table_size.size(); i++)
        nvals *= this->table_size[i];
    
    // calculate the product and update it to the local vector
    {
        this->table_plane_products.resize(this->table_size.size());
        FESystemUInt prod_val=1;
        for (FESystemUInt i=0; i<this->table_size.size(); i++)
        {
            if (i > 0)
                prod_val *= this->table_size[i-1];
            this->table_plane_products[i] = prod_val;
        }
    }
    
    
    // the total number of elements cannot be zero
    FESystemAssert0(nvals > 0,  FESystem::Utility::ZeroTableSize);
    
    // resize the data vector to the required size
    this->table_data.resize(nvals);
}


template <typename ValType>
void
FESystem::Utility::Table<ValType>::setVal(const std::vector<FESystemUInt>& el, ValType& v)
{
    FESystemUInt index = this->getIndex(el);
    this->table_data[index] = v;
}


template <typename ValType>
ValType& 
FESystem::Utility::Table<ValType>::getVal(const std::vector<FESystemUInt>& el)
{
    FESystemUInt index = this->getIndex(el);
    return this->table_data[index];
}




template <typename ValType>
const ValType& 
FESystem::Utility::Table<ValType>::getVal(const std::vector<FESystemUInt>& el) const
{
    FESystemUInt index = this->getIndex(el);
    return this->table_data[index];
}



template <typename ValType>
std::vector<ValType>& 
FESystem::Utility::Table<ValType>::getAllValues()
{
    return this->table_data;
}


template <typename ValType>
const std::vector<ValType>& 
FESystem::Utility::Table<ValType>::getAllValues() const
{
    return this->table_data;
}



template <typename ValType>
FESystemUInt
FESystem::Utility::Table<ValType>::getDimension() const
{
    return this->table_size.size();
}



template <typename ValType>
FESystemUInt 
FESystem::Utility::Table<ValType>::getNElementsInDimension(FESystemUInt dim) const
{
    // make sure that the dimensions is within the dimensionality of this table
    FESystemAssert2(dim < this->table_size.size(), FESystem::Exception::IndexOutOfBound, 
                    dim, this->table_size.size());
    return this->table_size[dim];
}


template <typename ValType>
const std::vector<FESystemUInt>&
FESystem::Utility::Table<ValType>::getSizeVector() const 
{
    return this->table_size;
}


template <typename ValType>
FESystemUInt 
FESystem::Utility::Table<ValType>::getIndex(const std::vector<FESystemUInt>& el) const
{
    // make sure that the size and dimesions are correct
    FESystemAssert2(el.size() == this->table_size.size(), FESystem::Exception::DimensionsDoNotMatch, 
                    el.size(), this->table_size.size());
    for (FESystemUInt i=0; i<el.size(); i++)
        FESystemAssert2(el[i] < this->table_size[i], FESystem::Exception::IndexOutOfBound, 
                        el.size(), this->table_size[i]);
    
    // the elements are stored in a column major format, and is used to get the index value
    FESystemUInt index=0, tmp_val=0;
    for (FESystemUInt i=0; i<this->table_size.size(); i++)
    {
        tmp_val = this->table_size.size()-i-1;
        FESystemAssert2(el[tmp_val] < this->table_size[tmp_val], FESystem::Exception::IndexOutOfBound, 
                        el[tmp_val], this->table_size[tmp_val]);
        index += el[tmp_val] * this->table_plane_products[tmp_val];
    }
    
    return index;
}


#endif // __fesystem_auto_ptr_vector_h__


