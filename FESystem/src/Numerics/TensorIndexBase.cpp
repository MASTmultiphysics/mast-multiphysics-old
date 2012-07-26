
/*
 *  TensorIndex.cpp
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/15/09.
 *  Copyright 2009 . All rights reserved.
 *
 */


// FESystem includes
#include "Numerics/TensorIndexBase.h"
#include "Base/FESystemExceptions.h"



FESystem::Numerics::TensorIndexBase::TensorIndexBase(const FESystemUInt d,
													 const FESystemUInt r):
dim(d),
rank(r)
{
	this->index.resize(r);
}



FESystem::Numerics::TensorIndexBase::~TensorIndexBase()
{
	this->index.clear();
}



void
FESystem::Numerics::TensorIndexBase::reset()
{
	std::fill(this->index.begin(), this->index.end(), 0); 
}


FESystemUInt
FESystem::Numerics::TensorIndexBase::getRank() const
{
	return this->rank;
}


FESystemUInt 
FESystem::Numerics::TensorIndexBase::getDimension() const
{
	return this->dim;
}


void 
FESystem::Numerics::TensorIndexBase::setIndex(const FESystemUInt r, const FESystemUInt d)
{
	FESystemAssert2(r < this->getRank(), 
					FESystem::Exception::IndexOutOfBound,
					r, this->getRank());
	FESystemAssert2(d < this->getDimension(), 
					FESystem::Exception::IndexOutOfBound,
					d, this->getDimension());
	this->index[r] = d;
}



void
FESystem::Numerics::TensorIndexBase::setIndex(const std::vector<FESystemUInt>& i)
{
	FESystemUInt given_rank = i.size();
	
	FESystemAssert2(given_rank < this->getRank(), 
					FESystem::Exception::IndexOutOfBound,
					given_rank, this->getRank());

	for (FESystemUInt rank_i=0; rank_i < given_rank; rank_i++)
	{
		FESystemAssert2(i[rank_i] < this->getDimension(), 
						FESystem::Exception::IndexOutOfBound,
						i[rank_i], this->getDimension());
		this->index[rank_i] = i[rank_i];	
	}	
}


const std::vector<FESystemUInt>&  
FESystem::Numerics::TensorIndexBase::getIndex() const
{
	return this->index;
}



FESystemUInt 
FESystem::Numerics::TensorIndexBase::getIndex(const FESystemUInt r) const
{
	FESystemAssert2(r < this->getRank(), 
					FESystem::Exception::IndexOutOfBound,
					r, this->getRank());
	return this->index[r];	
}

