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
class Array
{
public:
	using value_type=VALUETYPE;
	using storage_type=typename AbstractStorageType::template Impl<D>;
	using shape_type=typename storage_type::shape_type;
	using index_type=typename storage_type::index_type;
	using iterator_type=VALUETYPE*;
	using const_iterator_type=const VALUETYPE*;
	static const unsigned int num_dimensions=D;
	using abstract_storage_type=AbstractStorageType; // type-id is vector<T, Alloc<T>>
	
protected:
	storage_type _storage;
	impl::ContiguousMemView<value_type> _data;

public:
	Array(VALUETYPE* ptr,const storage_type& tstorage):
		_storage(tstorage),
		_data(ptr,_storage.size())
	{}
	Array(const storage_type& tstorage,const VALUETYPE& fval=VALUETYPE()):
		_storage(tstorage),
		_data(_storage.size(),fval)
	{}
	Array(VALUETYPE* ptr,const shape_type& tshape):
	_storage(tshape),
	_data(ptr,_storage.size())
	{}
	Array(const shape_type& tshape={},const VALUETYPE& fval=VALUETYPE()):
	_storage(tshape),
	_data(_storage.size(),fval)
	{}
	
	storage_type& storage() { return _storage; }
	const storage_type& storage() const { return _storage; }
	
	size_t size() const { return _data._size; }
	
	iterator_type begin() { return data(); }
	const_iterator_type begin() const { return _data._ptr; }
	
	iterator_type end() { return data()+_data._size; }
	const_iterator_type end() const { return data()+_data._size; }
	
	const value_type& operator[](size_t dex) const { return _data._ptr[dex]; }
	value_type& operator[](size_t dex) { return _data._ptr[dex]; }

	value_type* data() { return _data._ptr; }
	const value_type* data() const { return _data._ptr; }
	
	
	template<typename... Args>
	size_t ravel(Args&&... args) const {
		return _storage.ravel(std::forward<Args>(args)...);
	}
	
	index_type unravel(size_t dex) const {
		return _storage.unravel(dex);
	}
	index_type unravel(const iterator_type& at) const { //helper method for iterators
		return _storage.unravel(at-begin());
	}
	template<typename... Args>
	auto shape(Args&&... args) const -> decltype(_storage.shape(std::forward<Args>(args)...)) {
		return _storage.shape(std::forward<Args>(args)...);
	}
	
	template<class ...IndexTail> 
	const value_type& operator()(const IndexTail&... tail) const { 
		return operator[](ravel(tail...)); 
	}
	template<class ...IndexTail> 
	value_type& operator()(const IndexTail&... tail) { 
		return operator[](ravel(tail...)); 
	}

	template<class BinaryOp,class Barray>
	auto elementWise(BinaryOp f,const Barray& b) const->
		Array<decltype(f(operator[](0),b[0])),D,abstract_storage_type>
	{
		if(this->shape()!=b.shape()) throw std::runtime_error("The two arrays must have the same shape to have an elementWise operation");
		Array<decltype(f(operator[](0),b[0])),D,abstract_storage_type> out(this->shape());
		if(_storage==b.storage())
		{
			std::transform(begin(),end(),b.begin(),out.begin(),f);
		}
		else
		{
			size_t bN=b.size();
			for(size_t bi=0;bi<bN;bi++)
			{
				size_t ai=ravel(b.unravel(bi));
				out[ai]=f(operator[](ai),b[bi]);
			}
		}
		return out;
	}

	template<class BinaryOp,class ArrayB>
	Array& elementWiseInPlace(BinaryOp f,const ArrayB& b)
	{
		if(this->shape()!=b.shape()) throw std::runtime_error("The two arrays must have the same shape to have an elementWise operation");
		if(_storage==b.storage())
		{
			size_t bN=b.size();
			auto aptr=data();auto bptr=b.data();
			for(size_t i=0;i<bN;i++)
			{
				f(aptr[i],bptr[i]);
			}
		}
		else
		{
			size_t bN=b.size();
			for(size_t bi=0;bi<bN;bi++)
			{
				size_t ai=ravel(b.unravel(bi));
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
		std::transform(begin(),end(),out.begin(),f);
		return out;
	}
	template<class UnaryOp>
	Array& elementWiseInPlace(UnaryOp f)
	{
		std::for_each(begin(),end(),f);
		return *this;
	}
	
	//TODO make sure there's handling for D==0 or D==1
	Array<value_type,D-1,abstract_storage_type> pop_order_ref(unsigned int d) const
	{
		using Am1=Array<value_type,D-1,abstract_storage_type>;
		typename Am1::shape_type shp;
		const shape_type& mshp=shape();
		std::copy(std::begin(mshp),std::begin(mshp)+D-1,std::begin(shp));
		typename Am1::shape_type zed={};
		zed[D-1]=d;
		const value_type* slicebegin=data()+ravel(zed);
		return Array<value_type,D-1,abstract_storage_type>(const_cast<value_type*>(slicebegin),shp);
	}
};
}

#include "src/Operators.hpp"

#endif
