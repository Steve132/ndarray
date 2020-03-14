#ifndef NDARRAY_HPP
#define NDARRAY_HPP

#include<array>
#include<numeric>
#include<functional>
#include<vector>
#include<type_traits>
#include<algorithm>
#include<iterator>
#include<iostream>

namespace nd
{
namespace impl
{
template<class VALUETYPE>
struct ptr_or_vector
{
	std::vector<VALUETYPE> _vector;
public:
	VALUETYPE* _ptr;
	size_t _size;
	ptr_or_vector(VALUETYPE* ptr,size_t sz=0):_ptr(ptr),_size(sz)
	{}
	template<typename... Args>
	explicit ptr_or_vector(Args&&... args):
		_vector(std::forward<Args>(args)...),
		_ptr(_vector.data()),
		_size(_vector.size())
	{}
	ptr_or_vector(const ptr_or_vector& o): _vector(o._vector),
		_ptr(_vector.size() ? _vector.data() : o._ptr),
		_size(o._size)
	{}
	ptr_or_vector(ptr_or_vector&& o): _vector(std::move(o._vector)),
		_ptr(_vector.size() ? _vector.data() : o._ptr),
		_size(o._size)
	{}
	ptr_or_vector& operator=(ptr_or_vector&& o)
	{
		_vector=std::move(o._vector);
		_ptr=_vector.size() ? _vector.data() : o._ptr;_size=o._size;
		return *this;
	}
	ptr_or_vector& operator=(const ptr_or_vector& o)
	{
		_vector=o._vector;
		_ptr=_vector.size() ? _vector.data() : o._ptr;_size=o._size;
		return *this;
	}
};
template<unsigned int D>
class StorageOrderBase
{
public:
	typedef std::array<size_t,D> size_type;

protected:
	size_type _shape;
	StorageOrderBase(const size_type& tshape):_shape(tshape)
	{}
public:
	const size_type& shape() const { return _shape; }
	size_t shape(unsigned int d) const { return _shape[d]; }
};
}

struct ColMajorOrder
{

template<unsigned int D>
class computed_impl: public impl::StorageOrderBase<D>
{
public:
	using typename impl::StorageOrderBase<D>::size_type;
protected:
	size_type _strides;
	computed_impl(const size_type& tshape={}):
		impl::StorageOrderBase<D>(tshape)
	{
		const size_type& shp=impl::StorageOrderBase<D>::shape();
		std::partial_sum(
			std::begin(shp),
			std::begin(shp)+D,
			std::begin(_strides),
			[](size_t a,size_t b){ return a*b; } 
		);
	}
public:

	template<class IndexType>
	size_t linearize(const IndexType& index) const {
		
		
		return index[0]+std::inner_product(std::begin(index)+1,std::begin(index)+D,std::begin(_strides),size_t(0));
	}
	template<class IndexHead1,class IndexHead2,class... IndexTail>
	size_t linearize(const IndexHead1& index1,const IndexHead2& index2,const IndexTail&... tail) const {
		return linearize(std::array<IndexHead1,2+sizeof...(tail)>{index1,index2,tail...});
	}
	size_t size() const { return _strides[D-1]; }
	size_type unlinearize(size_t index) const {
		size_type out;
		auto outiter=std::rbegin(out);
		auto rend=std::rbegin(_strides)+D;
		for(auto stiter=std::rbegin(_strides)+1;stiter!=rend;stiter++,outiter++){	
			*(outiter)=index/(*stiter);
			index=index%(*stiter);
		}
		return out;
	}
};

};

template<class VALUETYPE,unsigned int D,class StorageType=ColMajorOrder>
class Array: public StorageType::template computed_impl<D>
{
public:
	typedef VALUETYPE value_type;
	typedef StorageType storage_type;
	using typename StorageType::template computed_impl<D>::size_type;
	static const unsigned int num_dimensions=D;
	
	using StorageType::template computed_impl<D>::linearize;
	using StorageType::template computed_impl<D>::unlinearize;
protected:
	impl::ptr_or_vector<value_type> _data;
public:
	Array(VALUETYPE* ptr,const size_type& shape):
		StorageType::template computed_impl<D>(shape),
		_data(ptr,StorageType::template computed_impl<D>::size())
	{}
	Array(const size_type& shape={},const VALUETYPE& fval=VALUETYPE()):
		StorageType::template computed_impl<D>(shape),
		_data(StorageType::template computed_impl<D>::size(),fval)
	{}
	
	value_type* data() { return _data._ptr; }
	const value_type* data() const { return _data._ptr; }
	
	value_type* begin() { return data(); }
	const value_type* begin() const { return _data._ptr; }
	
	value_type* end() { return data()+_data._size; }
	const value_type* end() const { return data()+_data._size; }
	
	
	const value_type& operator[](size_t dex) const { return _data._ptr[dex]; }
	value_type& operator[](size_t dex) { return _data._ptr[dex]; }

	template<class ...IndexTail> 
	const value_type& operator()(const IndexTail&... tail) const { 
		return operator[](linearize(tail...)); 
	}
	template<class ...IndexTail> 
	value_type& operator()(const IndexTail&... tail) { 
		return operator[](linearize(tail...)); 
	}

public:
	template<class BinaryOp,class B,class StB>
	auto elementWise(BinaryOp f,const Array<B,D,StB>& b) const->
	typename std::enable_if<std::is_same<StorageType,StB>::value,
		Array<decltype(f(operator[](0),b[0])),D,StorageType>
	>::type
	{
		if(this->shape()!=b.shape()) throw std::runtime_error("The two arrays must have the same shape to have an elementWise operation");
		Array<decltype(f(operator[](0),b[0])),D,StorageType> out(this->shape());
		std::transform(data(),data()+this->size(),b.data(),out.data(),f);
		return out;
	}
	template<class BinaryOp,class B,class StB>
	auto elementWise(BinaryOp f,const Array<B,D,StB>& b) const->
	typename std::enable_if<!std::is_same<StorageType,StB>::value,
	Array<decltype(f(operator[](0),b[0])),D,StorageType>
	>::type
	{
		if(this->shape()!=b.shape()) throw std::runtime_error("The two arrays must have the same shape to have an elementWise operation");
		Array<decltype(f(operator[](0),b[0])),D,StorageType> out(this->shape());
		size_t bN=b.size();
		
		for(size_t bi=0;bi<bN;bi++)
		{
			size_t ai=linearize(b.unlinearize(bi));
			out[ai]=f(operator[](ai),b[bi]);
		}
		return out;
	}
	
	template<class UnaryOp>
	auto elementWise(UnaryOp f) const->
		Array<decltype(f(operator[](0))),D,StorageType>
	{
		Array<decltype(f(operator[](0))),D,StorageType> out(this->shape());
		std::transform(data(),data()+this->size(),out.data(),f);
		return out;
	}
	template<class UnaryOp>
	void elementWiseInPlace(UnaryOp f)
	{
		std::transform(data(),data()+this->size(),data(),f);
	}
};

template<class VALUETYPE,class StorageType>
std::ostream& operator<<(std::ostream& out,const Array<VALUETYPE,1,StorageType>& arr)
{
	for(size_t i=0;i<arr.shape(0);i++) out << arr(i) << " ";
	return out;
}
template<class VALUETYPE,class StorageType>
std::ostream& operator<<(std::ostream& out,const Array<VALUETYPE,2,StorageType>& arr)
{
	for(size_t i=0;i<arr.shape(0);i++) 
	{
		for(size_t j=0;j<arr.shape(1);j++) out << arr(i,j) << " ";
		out << "\n";
	}
	return out;
}
/*
template<class VALUETYPE,class StorageType>
std::ostream& operator<<(const Array<VALUETYPE,2,StorageType>& arr,std::ostream& out)
{
	for(size_t i=0;i<arr.shape(1);i++) out << arr << "\n";
	return out;
}
template<class VALUETYPE,unsigned int D,class StorageType>
std::ostream& operator<<(const Array<VALUETYPE,D,StorageType>& arr,std::ostream& out)
{
	for(size_t i=0;i<arr.shape(2);i++) 
	{
		out << arr << "\n";
	}
	return out;
}*/


template<class VALUETYPE,unsigned int D,class StorageType=ColMajorOrder>
class LinearArray: public Array<VALUETYPE,D,StorageType>
{
protected:

public:
	using typename Array<VALUETYPE,D,StorageType>::value_type;
	using typename Array<VALUETYPE,D,StorageType>::size_type; 
	using typename Array<VALUETYPE,D,StorageType>::storage_type; 
	using Array<VALUETYPE,D,StorageType>::Array;
	
	using Array<VALUETYPE,D,StorageType>::operator[];
	using Array<VALUETYPE,D,StorageType>::operator();
	
	
	template<class ArrayB>
	auto operator+(const ArrayB& ba) const
	->Array<decltype(operator[](0)+ba[0]),D,StorageType>
	{
		return Array<VALUETYPE,D,StorageType>::elementWise(
			[](const VALUETYPE& v1,const typename ArrayB::value_type& v2){ return v1+v2; },
		ba);
	}
	template<class ArrayB>
	auto operator-(const ArrayB& ba) const
	->Array<decltype(operator[](0)+ba[0]),D,StorageType>
	{
		return Array<VALUETYPE,D,StorageType>::elementWise(
			[](const VALUETYPE& v1,const typename ArrayB::value_type& v2){ return v1+v2; },
														   ba);
	}
};



}

#endif
