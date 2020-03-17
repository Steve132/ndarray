#ifndef NDARRAY_ARRAY_HPP
#define NDARRAY_ARRAY_HPP

#include<utility>
#include<algorithm>
#include<stdexcept>
#include "src/ContiguousMemView.hpp"
#include "ColMajorOrder.hpp"


namespace nd
{
template<class VALUETYPE,unsigned int D,class AbstractStorageType=ColMajorOrder >
class Array: public AbstractStorageType::template Impl<VALUETYPE,D>
{
public:
	using storage_type=typename AbstractStorageType::template Impl<VALUETYPE,D>;
	using value_type=typename storage_type::value_type;
	using shape_type=typename storage_type::shape_type;
	using index_type=typename storage_type::index_type;
	using iterator_type=typename storage_type::iterator_type;
	using const_iterator_type=typename storage_type::const_iterator_type;
	static const unsigned int num_dimensions=D;
	using abstract_storage_type=AbstractStorageType; // type-id is vector<T, Alloc<T>>
public:
	storage_type& storage() { return *this; }
	const storage_type& storage() const { return *this; }

	using storage_type::operator[];
	using storage_type::storage_type;
	
	template<class ...IndexTail> 
	const value_type& operator()(const IndexTail&... tail) const { 
		return operator[](storage_type::ravel(tail...)); 
	}
	template<class ...IndexTail> 
	value_type& operator()(const IndexTail&... tail) { 
		return operator[](storage_type::ravel(tail...)); 
	}

	template<class BinaryOp,class Barray>
	auto elementWise(BinaryOp f,const Barray& b) const->
		Array<decltype(f(operator[](0),b[0])),D,abstract_storage_type>
	{
		if(this->shape()!=b.shape()) throw std::runtime_error("The two arrays must have the same shape to have an elementWise operation");
		Array<decltype(f(operator[](0),b[0])),D,abstract_storage_type> out(this->shape());
		if(storage_type::operator==(b.storage()))
		{
			for(size_t bi=0;bi<b.size();bi++)
			{
				out[bi]=f(operator[](bi),b[bi]);
			}
		}
		else
		{
			size_t bN=b.size();
			for(size_t bi=0;bi<bN;bi++)
			{
				size_t ai=storage_type::ravel(b.unravel(bi));
				out[ai]=f(operator[](ai),b[bi]);
			}
		}
		return out;
	}

	template<class BinaryOp,class ArrayB>
	Array& elementWiseInPlace(BinaryOp f,const ArrayB& b)
	{
		if(this->shape()!=b.shape()) throw std::runtime_error("The two arrays must have the same shape to have an elementWise operation");
		if(storage_type::operator==(b.storage()))
		{
			size_t bN=b.size();
			for(size_t i=0;i<bN;i++)
			{
				f(operator[](i),b[i]);
			}
		}
		else
		{
			size_t bN=b.size();
			for(size_t bi=0;bi<bN;bi++)
			{
				size_t ai=storage_type::ravel(b.unravel(bi));
				f(operator[](ai),b[bi]);
			}
		}
		return *this;
	}
	template<class UnaryOp>
	auto elementWise(UnaryOp f) const->
	Array<decltype(f(operator[](0))),D,abstract_storage_type>
	{
		Array<decltype(f(operator[](0))),D,abstract_storage_type> out(this->shape());
		size_t aN=storage_type::size();
		for(size_t ai=0;ai<aN;ai++)
		{
			out[ai]=f(operator[](ai));
		}
		return out;
	}
	template<class UnaryOp>
	Array& elementWiseInPlace(UnaryOp f)
	{
		size_t aN=storage_type::size();
		for(size_t ai=0;ai<aN;ai++)
		{
			f(operator[](ai));
		}
		return *this;
	}
	

};
}

#include "src/Operators.hpp"

#endif
