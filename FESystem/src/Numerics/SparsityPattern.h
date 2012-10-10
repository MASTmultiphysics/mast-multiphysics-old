//
//  SparsityPattern.h
//  FESystem
//
//  Created by Manav Bhatia on 4/26/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem_sparsity_pattern_h__
#define __fesystem_sparsity_pattern_h__

// C++ includes
#include <vector>
#include <iostream>
#include <map>
#include <set>


// FESystem includes
#include "Base/FESystemTypes.h"

namespace FESystem
{
    // Forward decleration
    namespace Base {class DegreeOfFreedomMap;}
    
    namespace Numerics
    {
        /*!
         *    Stores the details about the nonzero columns for each node in the sparse matrix
         */
        class SparsityPattern
        {
        public:
            
            /*!
             *   Constructor
             */
            SparsityPattern();
            
            virtual ~SparsityPattern();
               
            
            /*!
             *    clears the data structure for reinitialization
             */
            virtual void clear();
            
            /*!
             *   clears the data structure nonzero_column_ids_per_row already specified
             */
            void reinit();

            
            /*!
             *   Creates the sparsity pattern by union of the this with \p s. 
             */
            void unionWithSparsityPattern(const FESystem::Numerics::SparsityPattern& s);
            
            /*!
             *   Sets the number of degrees of freedom in the sparsity pattern
             */
            void setNDofs(FESystemUInt n);            
            
            /*!
             *   Adds multiple nonzero columns for row. This method is used for initializing the sparsity pattern, and should only be used before initialization.
             */  
            void addNonZeroColumnsForRow(FESystemUInt row, const std::set<FESystemUInt>& col);

            /*!
             *    returns the number of degrees of freedom
             */
            FESystemUInt getNDOFs() const;
            
            /*!
             *   returns the number of nonzero columns in the entire matrix
             */
            FESystemUInt getNNonzeroValues() const;

            
            /*!
             *   returns the number of nonzero columns in row
             */
            FESystemUInt getNNonzeroColumnsInRow(const FESystemUInt row) const;
            
            
            /*!
             *   returns the nonzero columns for the specified row. The IDs are provided in a set, so they are in increasing order
             */
            const std::vector<std::pair<FESystemUInt, FESystemUInt> >& getAllNonzeroColumnsInRow(const FESystemUInt row) const;

//            /*!
//             *   returns the nonzero rows for the specified column. The IDs are provided in a set, so they are in increasing order
//             */
//            const std::vector<std::pair<FESystemUInt, FESystemUInt> >& getAllNonzeroRowsInColumn(const FESystemUInt column) const;
                        
            /*!
             *   Returns the ID of the specified row and column in the contiguous vector of all nonzero values. It is an error if the 
             *   value does not exist in the matrix, which can be checked by the method checkIfValueExists.
             */
            void getIDForRowAndColumn(const FESystemUInt row, const FESystemUInt column, FESystemBoolean& if_exists, FESystemUInt& id_val) const;
            
            /*!
             *   initializes the sparsity pattern for the LU factored matrices by adding the fill in valuess
             */
            void initSparsityPatternForLUFactorization(FESystem::Numerics::SparsityPattern& l_factorization, FESystem::Numerics::SparsityPattern& u_factorization) const;

            /*!
             *   initializes the sparsity pattern after removing the constrained degrees of freedom
             */
            void initSparsityPatternForNonConstrainedDOFs(const std::vector<FESystemUInt>& nonbc_dofs, const std::map<FESystemUInt, FESystemUInt>& old_to_new_id_map,
                                                          FESystem::Numerics::SparsityPattern& constrained_sparsity_pattern) const;

            
            /*!
             *   Increments the iterator \p it so that it points to the column \p val
             */
            void find(std::vector<std::pair<FESystemUInt, FESystemUInt> >::const_iterator& it, std::vector<std::pair<FESystemUInt, FESystemUInt> >::const_iterator& end, FESystemUInt val) const;

            /*!
             *   Increments the iterator \p it so that it points to the first column >= \p val
             */
            void lowerBound(std::vector<std::pair<FESystemUInt, FESystemUInt> >::const_iterator& it, std::vector<std::pair<FESystemUInt, FESystemUInt> >::const_iterator& end, FESystemUInt val) const;

            /*!
             *    writes the data structure details to the output
             */
            void write(std::ostream& output) const;
             
            
        protected:
            
            
            // this is a friend class
            friend class FESystem::Base::DegreeOfFreedomMap;

            /*!
             *   if the data structure has been initialized
             */
            FESystemBoolean if_initialized;

            /*!
             *   number of nonzero values in the matrix
             */
            FESystemUInt n_nonzero_values;

            /*!
             *   For each row, the set stores the columns for which the values are nonzero
             */
            std::vector<std::vector<std::pair<FESystemUInt, FESystemUInt> > >  nonzero_column_ids_per_row;

//            /*!
//             *   For each column, the set stores the rows for which the values are nonzero. This is post processed from the 
//             *   same data for rows.
//             */
//            std::vector<std::vector<std::pair<FESystemUInt, FESystemUInt> > > nonzero_row_ids_per_column;
            
            /*!
             *   a temporary vector of booleans to track the number of nonzeros in a row
             */
            std::vector<FESystemBoolean> nonzero_cols_tmp;
        };
        
    }
}

#endif // __fesystem_sparsity_pattern_h__


