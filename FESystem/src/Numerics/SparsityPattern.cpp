//
//  SparsityPattern.cpp
//  FESystem
//
//  Created by Manav Bhatia on 6/13/12.
//  Copyright (c) 2012. All rights reserved.
//

// C++ includes
#include <iomanip>
#include <map>
#include <algorithm>

// FESystem includes
#include "Numerics/SparsityPattern.h"
#include "Base/FESystemExceptions.h"

// Metis includes
#include "metis.h"


FESystem::Numerics::SparsityPattern::SparsityPattern():
if_initialized(false),
n_nonzero_values(0)
{
    
}

FESystem::Numerics::SparsityPattern::~SparsityPattern()
{
    
}
            

FESystemUInt
FESystem::Numerics::SparsityPattern::getNDOFs() const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    return this->nonzero_column_ids_per_row.size();
}


void
FESystem::Numerics::SparsityPattern::setNDofs(FESystemUInt n)
{
    FESystemAssert0(!this->if_initialized, FESystem::Exception::InvalidState);
    this->nonzero_column_ids_per_row.resize(n);
    this->nonzero_cols_tmp.resize(n);
}



void
FESystem::Numerics::SparsityPattern::unionWithSparsityPattern(const FESystem::Numerics::SparsityPattern& s)
{
    FESystemAssert0(!this->if_initialized, FESystem::Exception::InvalidState);
    FESystemAssert0(s.nonzero_column_ids_per_row.size() == this->nonzero_column_ids_per_row.size(), FESystem::Exception::InvalidValue);
    
    std::vector<std::pair<FESystemUInt, FESystemUInt> >::const_iterator it, end;
    
    if (this->nonzero_cols_tmp.size() != this->nonzero_column_ids_per_row.size())
        this->nonzero_cols_tmp.resize(this->nonzero_column_ids_per_row.size());
    
    FESystemUInt row_num=0, n=0;
    for ( ; row_num!=this->nonzero_column_ids_per_row.size(); row_num++)
    {
        std::fill(this->nonzero_cols_tmp.begin(), this->nonzero_cols_tmp.end(), false);

        // add nonzeros from this pattern
        n = this->nonzero_column_ids_per_row[row_num].size();
        it = this->nonzero_column_ids_per_row[row_num].begin();
        end = this->nonzero_column_ids_per_row[row_num].end();        
        for (; it!=end; it++)
            this->nonzero_cols_tmp[it->first] = true;

        // add nonzeros from the given pattern
        it = s.nonzero_column_ids_per_row[row_num].begin();
        end = s.nonzero_column_ids_per_row[row_num].end();        
        for (; it!=end; it++)
            if (!this->nonzero_cols_tmp[it->first])
            {
                this->nonzero_cols_tmp[it->first] = true;
                n++;
            }

        // now put the results back in the sparsity pattern
        this->nonzero_column_ids_per_row[row_num].resize(n);
        n=0;
        for (FESystemUInt j=0; j<this->nonzero_cols_tmp.size(); j++)
            if (this->nonzero_cols_tmp[j])
                this->nonzero_column_ids_per_row[row_num][n++].first = j;
    }
}



void 
FESystem::Numerics::SparsityPattern::addNonZeroColumnsForRow(FESystemUInt row, const std::set<FESystemUInt>& col)
{
    FESystemAssert0(!this->if_initialized, FESystem::Exception::InvalidState);
    FESystemAssert1(row < this->nonzero_column_ids_per_row.size(), FESystem::Exception::InvalidID, row);
    
    std::vector<std::pair<FESystemUInt, FESystemUInt> >::const_iterator it, end;
    
    if (this->nonzero_cols_tmp.size() != this->nonzero_column_ids_per_row.size())
        this->nonzero_cols_tmp.resize(this->nonzero_column_ids_per_row.size());
    std::fill(this->nonzero_cols_tmp.begin(), this->nonzero_cols_tmp.end(), false);

    FESystemUInt n=0, first_nonz=0, last_nonz=0;
    
    // add nonzeros from this pattern
    n = this->nonzero_column_ids_per_row[row].size();    
    if (n > 0)
    {
        it = this->nonzero_column_ids_per_row[row].begin();
        end = this->nonzero_column_ids_per_row[row].end();
        for (; it!=end; it++)
            this->nonzero_cols_tmp[it->first] = true;

        first_nonz = this->nonzero_column_ids_per_row[row].begin()->first;
        last_nonz = this->nonzero_column_ids_per_row[row].rbegin()->first;
    }
    
    // add nonzeros from the given pattern
    std::set<FESystemUInt>::const_iterator s_it = col.begin(), s_end = col.end();        
    for (; s_it!=s_end; s_it++)
        if (!this->nonzero_cols_tmp[*s_it])
        {
            this->nonzero_cols_tmp[*s_it] = true;
            if (*s_it < first_nonz)
                first_nonz = *s_it;
            if (*s_it > last_nonz)
                last_nonz = *s_it;
            n++;
        }
    
    // now put the results back in the sparsity pattern
    this->nonzero_column_ids_per_row[row].resize(n);
    n=0;
    for (FESystemUInt j=first_nonz; j<=last_nonz; j++)
        if (this->nonzero_cols_tmp[j])
            this->nonzero_column_ids_per_row[row][n++].first = j;
}


void
FESystem::Numerics::SparsityPattern::clear()
{
    this->nonzero_column_ids_per_row.clear();
    //this->nonzero_row_ids_per_column.clear();
    if_initialized = false;
    this->n_nonzero_values = 0;
}



FESystemUInt
FESystem::Numerics::SparsityPattern::getNNonzeroValues() const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    return this->n_nonzero_values;
}


FESystemUInt
FESystem::Numerics::SparsityPattern::getNNonzeroColumnsInRow(const FESystemUInt row) const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    FESystemAssert2(row < this->nonzero_column_ids_per_row.size(), FESystem::Exception::IndexOutOfBound, row, this->nonzero_column_ids_per_row.size());
    
    return this->nonzero_column_ids_per_row[row].size();
}




const std::vector<std::pair<FESystemUInt, FESystemUInt> >&
FESystem::Numerics::SparsityPattern::getAllNonzeroColumnsInRow(const FESystemUInt row) const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    FESystemAssert2(row < this->nonzero_column_ids_per_row.size(), FESystem::Exception::IndexOutOfBound, row, this->nonzero_column_ids_per_row.size());
    
    return this->nonzero_column_ids_per_row[row];
}



void
FESystem::Numerics::SparsityPattern::getIDForRowAndColumn(const FESystemUInt row, const FESystemUInt column, FESystemBoolean& if_exists, FESystemUInt& id_val) const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    FESystemAssert2(row < this->nonzero_column_ids_per_row.size(), FESystem::Exception::IndexOutOfBound, row, this->nonzero_column_ids_per_row.size());
    
    std::vector<std::pair<FESystemUInt, FESystemUInt> >::const_iterator it=this->nonzero_column_ids_per_row[row].begin(), end=this->nonzero_column_ids_per_row[row].end();
    this->find(it, end, column);
    if (it != end)
    {
        id_val = it->second;
        if_exists = true;
    }
    else 
        if_exists = false;
}



void
FESystem::Numerics::SparsityPattern::find(std::vector<std::pair<FESystemUInt, FESystemUInt> >::const_iterator& it, 
                                          std::vector<std::pair<FESystemUInt, FESystemUInt> >::const_iterator& end, FESystemUInt val) const
{
    while (it!=end)
        if (it->first == val)
            break;
        else
            it++;
}


void
FESystem::Numerics::SparsityPattern::lowerBound(std::vector<std::pair<FESystemUInt, FESystemUInt> >::const_iterator& it, 
                                                std::vector<std::pair<FESystemUInt, FESystemUInt> >::const_iterator& end, FESystemUInt val) const
{
    while (it!=end)
        if (it->first >= val)
            break;
        else
            it++;
}



void
FESystem::Numerics::SparsityPattern::reinit()
{
    FESystemAssert0(this->if_initialized == false, FESystem::Exception::InvalidState);
    FESystemAssert0(this->nonzero_column_ids_per_row.size() > 0, FESystem::Exception::InvalidState); // this structure must be initialized
    
    const FESystemUInt n_dofs = this->nonzero_column_ids_per_row.size();
    
    // now process the ID data: this ends up putting variables in a row major format
    std::vector<std::vector<std::pair<FESystemUInt, FESystemUInt> > >::iterator it=this->nonzero_column_ids_per_row.begin(), end=this->nonzero_column_ids_per_row.end();

    // boolean to store if a row and column is active: column major format
    std::vector<std::pair<FESystemUInt, FESystemUInt> > nonzeros_per_col(n_dofs);
    std::vector<std::pair<FESystemInt, FESystemInt> > col_row_limits(n_dofs);
    std::fill(col_row_limits.begin(), col_row_limits.end(), std::pair<FESystemUInt, FESystemUInt>(n_dofs+1,-1)); // keep easily recognizable limits to find if the value was set

    // set the internal ids for the values
    {        
        FESystemInt n=0;
        this->n_nonzero_values = 0;
        std::vector<std::pair<FESystemUInt, FESystemUInt> >::iterator s_it, s_end;
        for ( ; it != end; it++)
        {
            for ( s_it=it->begin(); s_it != it->end(); s_it++)
            {
                s_it->second = this->n_nonzero_values++;

                if (col_row_limits[s_it->first].first > n) col_row_limits[s_it->first].first=n; // defines the lower limit on row for a column
                if (col_row_limits[s_it->first].second < n) col_row_limits[s_it->first].second=n; // defines the upper limit on row for a column
            }
            n++;
        }
    }
    
    this->nonzero_cols_tmp.clear(); // not needed anymore, so can be cleared
    this->if_initialized = true;
}



void 
FESystem::Numerics::SparsityPattern::initSparsityPatternForLUFactorization(FESystem::Numerics::SparsityPattern& l_factorization, FESystem::Numerics::SparsityPattern& u_factorization) const
{
    std::vector<std::pair<FESystemUInt, FESystemUInt> >::const_iterator s_it, s_end;
    std::vector<FESystemUInt>::const_iterator it2, it2_end;
    FESystemUInt n_dofs = this->nonzero_column_ids_per_row.size();
    
    // create a copy, and then use that as the starting point
    std::vector<FESystemBoolean> combined_nonzero_column_bools(n_dofs);
    std::vector<std::vector<FESystemUInt> > combined_nonzero_ids(n_dofs);
    
    // the n^th row gets propagated down for each row that has a nonzero n^th column
    FESystemUInt n_nonzero_cols=0, first_col, last_col;
    for (FESystemUInt n=0; n<n_dofs; n++)
    {
        std::fill(combined_nonzero_column_bools.begin(), combined_nonzero_column_bools.end(), false);
        n_nonzero_cols = 0;
        first_col = 0;
        last_col = 0;
        
        // add the nonzeros from the original matrix
        s_it = this->nonzero_column_ids_per_row[n].begin();
        s_end = this->nonzero_column_ids_per_row[n].end();
        n_nonzero_cols = this->nonzero_column_ids_per_row[n].size();
        for ( ;s_it != s_end; s_it++)
            combined_nonzero_column_bools[s_it->first] = true;
        first_col = this->nonzero_column_ids_per_row[n].begin()->first;
        last_col = this->nonzero_column_ids_per_row[n].rbegin()->first;
        
        // if n^th row has a nonzero in the n2^th column, then it gets the nonzeros from n2
        for (FESystemUInt n2=this->nonzero_column_ids_per_row[n].begin()->first; n2<n; n2++)
            if (combined_nonzero_column_bools[n2])
            {
                it2 = combined_nonzero_ids[n2].begin();
                it2_end = combined_nonzero_ids[n2].end();
                it2 = std::find(it2, it2_end, n2);
                FESystemAssert1(it2!=it2_end, FESystem::Exception::InvalidID, n2);
                // increment the iterator it2 since the value after n2 will be added as nonzeros in the n^th row
                it2++;
                for (; it2!=it2_end; it2++)
                    if (!combined_nonzero_column_bools[*it2])
                    {
                        combined_nonzero_column_bools[*it2] = true;
                        if (*it2 < first_col)
                            first_col = *it2;
                        if (*it2 > last_col)
                            last_col = *it2;
                        n_nonzero_cols++;
                    }
            }
        // add these nonzeros to the data structure
        combined_nonzero_ids[n].resize(n_nonzero_cols);
        n_nonzero_cols = 0;
        for (FESystemUInt j=first_col; j<=last_col; j++)
            if (combined_nonzero_column_bools[j])
                combined_nonzero_ids[n][n_nonzero_cols++] = j;
    }
    
    // now distribute across l and u
    l_factorization.nonzero_column_ids_per_row.resize(n_dofs);
    u_factorization.nonzero_column_ids_per_row.resize(n_dofs);

    FESystemUInt n_lower_t, n_upper_t;
    for (FESystemUInt n=0; n<n_dofs; n++)
    {
        n_lower_t=1; // both atleast have a diagonal entry
        n_upper_t=1;
        // count the number of nonzeros
        for (FESystemUInt j=0; j<combined_nonzero_ids[n].size(); j++)
            if (combined_nonzero_ids[n][j]<n)
                n_lower_t++;
            else if (combined_nonzero_ids[n][j]>n)
                n_upper_t++;
        
        // add the diagonal terms
        l_factorization.nonzero_column_ids_per_row[n].resize(n_lower_t);
        u_factorization.nonzero_column_ids_per_row[n].resize(n_upper_t);
        
        for (FESystemUInt j=0; j<n_lower_t; j++)
            l_factorization.nonzero_column_ids_per_row[n][j].first = combined_nonzero_ids[n][j];

        for (FESystemUInt j=(n_lower_t-1); j<(n_lower_t+n_upper_t-1); j++)
            u_factorization.nonzero_column_ids_per_row[n][j-(n_lower_t-1)].first = combined_nonzero_ids[n][j];
        
        // clear this to save memory
        combined_nonzero_ids[n].clear();
    }

    // tell the two sparsity patterns to initialize themselves
    l_factorization.reinit();
    u_factorization.reinit();
}




void 
FESystem::Numerics::SparsityPattern::initSparsityPatternForNonConstrainedDOFs(const std::vector<FESystemUInt>& nonbc_dofs, const std::map<FESystemUInt, FESystemUInt>& old_to_new_id_map,
                                                                              FESystem::Numerics::SparsityPattern& constrained_sparsity_pattern) const
{
    // it is assumed that the vector nonbc_dofs is already sorted in ascending order
    constrained_sparsity_pattern.clear();
    constrained_sparsity_pattern.nonzero_column_ids_per_row.resize(nonbc_dofs.size());
    
    // first create a map of old ids to the new ones
    std::vector<FESystemUInt> nonzeros(nonbc_dofs.size());
    std::vector<std::pair<FESystemUInt, FESystemUInt> >::const_iterator nonz_it, nonz_end;
    std::vector<FESystemUInt>::const_iterator it1, end, col_it, col_it2, col_end;
    std::map<FESystemUInt, FESystemUInt>::const_iterator map_row_it, map_col_it, map_end = old_to_new_id_map.end();    
    
    FESystemUInt n_nonzeros=0, constrained_row_num=0;
    end = nonbc_dofs.end();
    
    for (it1=nonbc_dofs.begin(); it1 != end; it1++)
    {
        std::fill(nonzeros.begin(), nonzeros.begin()+n_nonzeros, 0);
        n_nonzeros=0;
        
        nonz_it = this->nonzero_column_ids_per_row[*it1].begin();
        nonz_end = this->nonzero_column_ids_per_row[*it1].end();

        col_it = nonbc_dofs.begin(); 
        col_end = nonbc_dofs.end();

        for ( ; nonz_it!=nonz_end; nonz_it++)
        {
            col_it2 = std::find(col_it, col_end, nonz_it->first);
            // if the column in source sparsity is in the requested set of dofs, then write it
            if (col_it2 != col_end) // dof was found in the requested set
            {
                map_col_it = old_to_new_id_map.find(*col_it2);
                FESystemAssert1(map_col_it != map_end, FESystem::Exception::InvalidID, *col_it); // make sure that the id exists in the map
                nonzeros[n_nonzeros++] = map_col_it->second; // the ID of this dof in the reduced set
                col_it = col_it2; // use this as the beginning iterator for the next search
            }
        }
        
        // now add the nonzeros to the constrained sparsity
        constrained_row_num = old_to_new_id_map.find(*it1)->second;
        constrained_sparsity_pattern.nonzero_column_ids_per_row[constrained_row_num].resize(n_nonzeros);
        for (FESystemUInt col_num=0; col_num<n_nonzeros; col_num++)
            constrained_sparsity_pattern.nonzero_column_ids_per_row[constrained_row_num][col_num].first = nonzeros[col_num];
    }

    constrained_sparsity_pattern.reinit();
}



            
void
FESystem::Numerics::SparsityPattern::write(std::ostream &output) const
{
    output << "Sparsity Pattern Output: " << std::endl;
    output << "Number of Non-zero Values = " << this->getNNonzeroValues() << std::endl;
    output << "Nonzero Columns In Row:" << std::endl;
    
    std::vector<std::vector<std::pair<FESystemUInt, FESystemUInt> > >::const_iterator it=this->nonzero_column_ids_per_row.begin(), end=this->nonzero_column_ids_per_row.end();
    std::vector<std::pair<FESystemUInt, FESystemUInt> >::const_iterator s_it, s_end;
     
    FESystemUInt i=0;
    // add all the nonzero columns per row
    for ( ; it != end; it++)
    {
        s_it = it->begin(); s_end = it->end();
        
        output << std::setw(10) << i << " : "  << std::setw(5) << it->size() << " :  ";
        for ( ; s_it != s_end; s_it++)
            output << std::setw(10) << s_it->first;
        output << std::endl;
        
        // increment the counter
        i++;
    }
}


void
FESystem::Numerics::SparsityPattern::calculateFillReducingOrdering(std::vector<FESystemUInt> &reordered_dofs)
{
    FESystemInt n_dofs = this->nonzero_column_ids_per_row.size(), adjncy_id=0;
    for (FESystemInt i=0; i<this->nonzero_column_ids_per_row.size(); i++)
        adjncy_id += (this->nonzero_column_ids_per_row[i].size()-1); // minus 1 to omit itself from the list
    
    // create the adjacency array
    int *xadj=new int[n_dofs+1], *adjncy=new int[adjncy_id], *perm=new int[n_dofs], *iperm=new int[n_dofs], *options=new int[METIS_NOPTIONS];

    for (FESystemInt i=0; i<n_dofs+1; i++) xadj[i] = 0;
    for (FESystemInt i=0; i<adjncy_id; i++) adjncy[i] = 0;
    for (FESystemInt i=0; i<n_dofs; i++) {perm[i] = 0; iperm[i] = 0;}

    xadj[0] = 0; adjncy_id=0;
    // initialize the adjacency array data
    for (FESystemInt i=1; i<n_dofs+1; i++)
    {
        xadj[i] = xadj[i-1]+this->nonzero_column_ids_per_row[i-1].size()-1;
        
        for (FESystemInt j=0; j<this->nonzero_column_ids_per_row[i-1].size(); j++)
            if (this->nonzero_column_ids_per_row[i-1][j].first != i-1)
                adjncy[adjncy_id++] = this->nonzero_column_ids_per_row[i-1][j].first;
    }
    
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_NUMBERING] = 0;
    
    FESystemInt metis_return = METIS_NodeND(&n_dofs, &xadj[0], &adjncy[0], NULL, NULL, &perm[0], &iperm[0]);
    FESystemAssert1(metis_return==METIS_OK, FESystem::Exception::InvalidID, metis_return);
    
    // copy the ids back to the return vector
    for (FESystemInt i=0; i<n_dofs; i++)
        reordered_dofs[i] = iperm[i];

    delete[] xadj;
    delete[] adjncy;
    delete[] perm;
    delete[] iperm;
    delete[] options;

}



