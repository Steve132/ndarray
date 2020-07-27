#ifndef NDARRAY_ARRAY_HPP
#define NDARRAY_ARRAY_HPP

#include<utility>
#include<algorithm>
#include<stdexcept>
#include "src/ContiguousOrderBase.hpp"
#include "RowMajorOrder.hpp"


namespace nd
{
template<class VALUETYPE,unsigned int D,class OrderType=RowMajorOrder>
class Array: 
	public OrderType::template Layout<D>, 
	public impl::ContiguousOrderBaseImpl<VALUETYPE>
{
public:
	using order_type=OrderType;
	using layout_type=typename order_type::template Layout<D>;
	using storage_type=impl::ContiguousOrderBaseImpl<VALUETYPE>;
	
	using shape_type=typename layout_type::shape_type;
	using coord_type=typename layout_type::coord_type;
	using index_type=typename layout_type::index_type;
	
	using value_type=typename storage_type::value_type;
	using iterator_type=typename storage_type::iterator_type;
	using const_iterator_type=typename storage_type::const_iterator_type;

	
	static const unsigned int num_dimensions=D;
	
public:
	layout_type& layout() { return *this; }
	const layout_type& layout() const { return *this; }

	using storage_type::operator[];

	Array(const shape_type& tshape={}):
			layout_type(tshape),
			storage_type(layout_type::num_elements())
	{}
	template<class ...StorageArgs>
	Array(const shape_type& tshape,StorageArgs ...sa):
		layout_type(tshape),
		storage_type(layout_type::num_elements(),sa...)
	{}
	
	template<class T2>
	Array(const Array<T2,D,OrderType>& arr2):Array(arr2.shape())
	{
		std::copy(arr2.begin(),arr2.end(),std::begin(*this));
	}
	
	
	template<class ...IndexTail> 
	typename std::enable_if<(sizeof...(IndexTail))==D,const value_type&>::type
	operator()(const IndexTail&... tail) const { 
		return operator[](layout_type::ravel(tail...)); 
	}
	template<class ...IndexTail>
	typename std::enable_if<(sizeof...(IndexTail))==D,value_type&>::type
	operator()(const IndexTail&... tail) { 
		return operator[](layout_type::ravel(tail...)); 
	}
	
	using layout_type::unravel;
	coord_type unravel(const iterator_type& at) const { //helper method for iterators
		return unravel(at-impl::ContiguousOrderBaseImpl<VALUETYPE>::begin());
	}
	
	template<class BinaryOp,class ArrayB>
	auto elementWise(BinaryOp f,const ArrayB& b) const->
		typename std::enable_if<std::is_same<layout_type,typename ArrayB::layout_type>::value,
			Array<decltype(f(operator[](0),b[0])),D,order_type>
		>::type
	{
		if(this->shape()!=b.shape()) throw std::runtime_error("The two arrays must have the same shape to have an elementWise operation");
		Array<decltype(f(operator[](0),b[0])),D,order_type> out(this->shape());
		size_t bN=b.size();
		for(size_t bi=0;bi<b.size();bi++)
		{
			out[bi]=f(operator[](bi),b[bi]);
		}
		return out;
	}
	template<class BinaryOp,class ArrayB>
	auto elementWise(BinaryOp f,const ArrayB& b) const->
		typename std::enable_if<!std::is_same<layout_type,typename ArrayB::layout_type>::value,
			Array<decltype(f(operator[](0),b[0])),D,order_type>
		>::type
	{
		if(this->shape()!=b.shape()) throw std::runtime_error("The two arrays must have the same shape to have an elementWise operation");
		Array<decltype(f(operator[](0),b[0])),D,order_type> out(this->shape());
		size_t bN=b.size();
		for(size_t bi=0;bi<bN;bi++)
		{
			size_t ai=layout_type::ravel(b.unravel(bi));
			out[ai]=f(operator[](ai),b[bi]);
		}
		return out;
	}


	template<class BinaryOp,class ArrayB>
	auto elementWiseInPlace(BinaryOp f,const ArrayB& b)->
		typename std::enable_if<std::is_same<layout_type,typename ArrayB::layout_type>::value,
			Array&
		>::type
	{
		if(this->shape()!=b.shape()) throw std::runtime_error("The two arrays must have the same shape to have an elementWise operation");
		size_t bN=b.size();
		for(size_t i=0;i<bN;i++)
		{
			f(operator[](i),b[i]);
		}
		return *this;
	}
	template<class BinaryOp,class ArrayB>
	auto elementWiseInPlace(BinaryOp f,const ArrayB& b)->
		typename std::enable_if<!std::is_same<layout_type,typename ArrayB::layout_type>::value,
			Array&
		>::type
	{
		if(this->shape()!=b.shape()) throw std::runtime_error("The two arrays must have the same shape to have an elementWise operation");
		size_t bN=b.size();
		for(size_t bi=0;bi<bN;bi++)
		{
			size_t ai=layout_type::ravel(b.unravel(bi));
			f(operator[](ai),b[bi]);
		}
		return *this;
	}

	template<class UnaryOp>
	auto elementWise(UnaryOp f) const->
	Array<decltype(f(operator[](0))),D,order_type>
	{
		Array<decltype(f(operator[](0))),D,order_type> out(this->shape());
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
	
	void permuteInPlace(const shape_type& pvals) const
	{
		size_t aN=storage_type::size();
		for(size_t ai=0;ai<aN;ai++)
		{
			coord_type dex=layout_type::unravel(ai);
			coord_type new_dex;
			for(unsigned di=0;di<D;di++) { new_dex[di]=dex[pvals[di]]; }
			std::swap(operator[](layout_type::ravel(new_dex)),dex);
		}
	}
};
}

#include "src/Operators.hpp"

#endif
