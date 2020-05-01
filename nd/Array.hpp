#ifndef NDARRAY_ARRAY_HPP
#define NDARRAY_ARRAY_HPP

#include<utility>
#include<algorithm>
#include<stdexcept>
#include "src/ContiguousMemView.hpp"
#include "RowMajorOrder.hpp"


namespace nd
{
template<class VALUETYPE,unsigned int D,class AbstractOrderType=RowMajorOrder>
class Array: public AbstractOrderType::template Impl<VALUETYPE,D>
{
public:
	using order_type=typename AbstractOrderType::template Impl<VALUETYPE,D>;
	using value_type=typename order_type::value_type;
	using shape_type=typename order_type::shape_type;
	using index_type=typename order_type::index_type;
	using iterator_type=typename order_type::iterator_type;
	using const_iterator_type=typename order_type::const_iterator_type;
	using abstract_order_type=AbstractOrderType;
	
	static const unsigned int num_dimensions=D;
	
public:
	order_type& storage() { return *this; }
	const order_type& storage() const { return *this; }

	using order_type::operator[];
	using order_type::order_type;
	
	template<class T2>
	Array(const Array<T2,D,AbstractOrderType>& arr2):Array(arr2.shape())
	{
		std::copy(arr2.begin(),arr2.end(),std::begin(*this));
	}
	
	
	template<class ...IndexTail> 
	const value_type& operator()(const IndexTail&... tail) const { 
		return operator[](order_type::ravel(tail...)); 
	}
	template<class ...IndexTail> 
	value_type& operator()(const IndexTail&... tail) { 
		return operator[](order_type::ravel(tail...)); 
	}

	template<class BinaryOp,class Barray>
	auto elementWise(BinaryOp f,const Barray& b) const->
		Array<decltype(f(operator[](0),b[0])),D,abstract_order_type>
	{
		if(this->shape()!=b.shape()) throw std::runtime_error("The two arrays must have the same shape to have an elementWise operation");
		Array<decltype(f(operator[](0),b[0])),D,abstract_order_type> out(this->shape());
		if(order_type::operator==(b.storage()))
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
				size_t ai=order_type::ravel(b.unravel(bi));
				out[ai]=f(operator[](ai),b[bi]);
			}
		}
		return out;
	}

	template<class BinaryOp,class ArrayB>
	Array& elementWiseInPlace(BinaryOp f,const ArrayB& b)
	{
		if(this->shape()!=b.shape()) throw std::runtime_error("The two arrays must have the same shape to have an elementWise operation");
		if(order_type::operator==(b.storage()))
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
				size_t ai=order_type::ravel(b.unravel(bi));
				f(operator[](ai),b[bi]);
			}
		}
		return *this;
	}
	template<class UnaryOp>
	auto elementWise(UnaryOp f) const->
	Array<decltype(f(operator[](0))),D,abstract_order_type>
	{
		Array<decltype(f(operator[](0))),D,abstract_order_type> out(this->shape());
		size_t aN=order_type::size();
		for(size_t ai=0;ai<aN;ai++)
		{
			out[ai]=f(operator[](ai));
		}
		return out;
	}
	template<class UnaryOp>
	Array& elementWiseInPlace(UnaryOp f)
	{
		size_t aN=order_type::size();
		for(size_t ai=0;ai<aN;ai++)
		{
			f(operator[](ai));
		}
		return *this;
	}
	
	void permuteInPlace(const shape_type& pvals) const
	{
		size_t aN=order_type::size();
		for(size_t ai=0;ai<aN;ai++)
		{
			index_type dex=order_type::unravel(ai);
			index_type new_dex;
			for(unsigned di=0;di<D;di++) { new_dex[di]=dex[pvals[di]]; }
			std::swap(operator[](order_type::ravel(new_dex)),dex);
		}
	}
};
}

#include "src/Operators.hpp"

#endif
