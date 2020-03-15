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

template<unsigned int D>
class ColMajorOrder: public impl::StorageOrderBase<D>
{
public:
	using typename impl::StorageOrderBase<D>::size_type;
protected:
	size_type _strides;
public:
	ColMajorOrder(const size_type& tshape={}):
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

	template<class IndexType>
	size_t linearize(const IndexType& index) const {
		return index[0]+std::inner_product(std::begin(index)+1,std::begin(index)+D,std::begin(_strides),size_t(0));
	}
	template<class IndexHead1,class IndexHead2,class... IndexTail>
	size_t linearize(const IndexHead1& index1,const IndexHead2& index2,const IndexTail&... tail) const {
		return linearize(std::array<IndexHead1,2+sizeof...(tail)>{index1,index2,tail...});
	}
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
	size_t size() const { return _strides[D-1]; }
};

template<class VALUETYPE,unsigned int D,template<unsigned int> class StorageType=ColMajorOrder >
class Array
{
public:
	typedef VALUETYPE value_type;
	typedef StorageType<D> storage_type;
	typedef typename storage_type::size_type size_type;
	static const unsigned int num_dimensions=D;
	
protected:
	storage_type _storage;
	impl::ptr_or_vector<value_type> _data;

public:
	Array(VALUETYPE* ptr,const storage_type& tstorage):
		_storage(tstorage),
		_data(ptr,_storage.size())
	{}
	Array(const storage_type& tstorage,const VALUETYPE& fval=VALUETYPE()):
		_storage(tstorage),
		_data(_storage.size(),fval)
	{}
	Array(VALUETYPE* ptr,const size_type& tshape):
	_storage(tshape),
	_data(ptr,_storage.size())
	{}
	Array(const size_type& tshape={},const VALUETYPE& fval=VALUETYPE()):
	_storage(tshape),
	_data(_storage.size(),fval)
	{}
	
	storage_type& storage() { return _storage; }
	const storage_type& storage() const { return _storage; }
	
	size_t size() const { return _data._size; }
	
	value_type* data() { return _data._ptr; }
	const value_type* data() const { return _data._ptr; }
	
	value_type* begin() { return data(); }
	const value_type* begin() const { return _data._ptr; }
	
	value_type* end() { return data()+_data._size; }
	const value_type* end() const { return data()+_data._size; }
	
	
	const value_type& operator[](size_t dex) const { return _data._ptr[dex]; }
	value_type& operator[](size_t dex) { return _data._ptr[dex]; }
	
	template<typename... Args>
	size_t linearize(Args&&... args) const {
		return _storage.linearize(std::forward<Args>(args)...);
	}
	template<typename... Args>
	size_type unlinearize(Args&&... args) const {
		return _storage.linearize(std::forward<Args>(args)...);
	}
	template<typename... Args>
	auto shape(Args&&... args) const -> decltype(_storage.shape(std::forward<Args>(args)...)) {
		return _storage.shape(std::forward<Args>(args)...);
	}

	template<class ...IndexTail> 
	const value_type& operator()(const IndexTail&... tail) const { 
		return operator[](linearize(tail...)); 
	}
	template<class ...IndexTail> 
	value_type& operator()(const IndexTail&... tail) { 
		return operator[](linearize(tail...)); 
	}

public:
	template<class BinaryOp,class Barray>
	auto elementWise(BinaryOp f,const Barray& b) const->
	typename std::enable_if<std::is_same<storage_type,typename Barray::storage_type>::value,
		Array<decltype(f(operator[](0),b[0])),D,StorageType>
	>::type
	{
		if(this->shape()!=b.shape()) throw std::runtime_error("The two arrays must have the same shape to have an elementWise operation");
		Array<decltype(f(operator[](0),b[0])),D,StorageType> out(this->shape());
		std::transform(begin(),end(),b.begin(),out.begin(),f);
		return out;
	}
	template<class BinaryOp,class Barray>
	auto elementWise(BinaryOp f,const Barray& b) const->
	typename std::enable_if<!std::is_same<storage_type,typename Barray::storage_type>::value,
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
	
	template<class BinaryOp,class Barray>
	auto elementWiseInPlace(BinaryOp f,const Barray& b)->
	typename std::enable_if<
		std::is_same<storage_type,typename Barray::storage_type>::value,Array&
	>::type
	{
		if(this->shape()!=b.shape()) throw std::runtime_error("The two arrays must have the same shape to have an elementWise operation");
		auto abegin=begin();auto aend=end();auto bbegin=b.begin();auto bend=b.end();
		while(abegin!=aend)
		{
			f(*abegin++,*bbegin++);
		}
		return *this;
	}
	template<class BinaryOp,class Barray>
	auto elementWiseInPlace(BinaryOp f,const Barray& b)->
	typename std::enable_if<
		!std::is_same<storage_type,typename Barray::storage_type>::value,Array&
	>::type
	{
		if(this->shape()!=b.shape()) throw std::runtime_error("The two arrays must have the same shape to have an elementWise operation");
		size_t bN=b.size();
		
		for(size_t bi=0;bi<bN;bi++)
		{
			size_t ai=linearize(b.unlinearize(bi));
			f(operator[](ai),b[bi]);
		}
		return *this;
	}
	template<class UnaryOp>
	auto elementWise(UnaryOp f) const->
	Array<decltype(f(operator[](0))),D,StorageType>
	{
		Array<decltype(f(operator[](0))),D,StorageType> out(this->shape());
		std::transform(begin(),end(),out.begin(),f);
		return out;
	}
	template<class UnaryOp>
	Array& elementWiseInPlace(UnaryOp f)
	{
		std::for_each(begin(),end(),f);
		return *this;
	}
};

namespace impl
{
//horrible.
template<template<class,unsigned int,template<unsigned int> class> class ArrayClass,class BinaryResult,unsigned int D,template<unsigned int> class Stg>
struct ElementWiseMeta{
	typedef ArrayClass<BinaryResult,D,Stg> type;
};
}

template<class ArrayClass,typename std::enable_if<ArrayClass::num_dimensions==1,int>::type = 0>
std::ostream& operator<<(std::ostream& out,const ArrayClass& arr)
{
	for(size_t i=0;i<arr.shape(0);i++) out << arr(i) << " ";
	return out;
}
template<class ArrayClass,typename std::enable_if<ArrayClass::num_dimensions==2,int>::type = 0>
std::ostream& operator<<(std::ostream& out,const ArrayClass& arr)
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

/*TODO: this should be an expression template*/
template<class VALUETYPE,unsigned int D,template<unsigned int> class StorageType=ColMajorOrder >
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
	LinearArray(ArrayB&& o):Array<VALUETYPE,D,StorageType>(std::move(o))
	{}
	template<class ArrayB>
	LinearArray(const ArrayB& o):Array<VALUETYPE,D,StorageType>(o)
	{}
	template<class ArrayB>
	LinearArray& operator=(ArrayB&& o)
	{
		Array<VALUETYPE,D,StorageType>::operator=(std::move(o));
		return *this;
	}
	template<class ArrayB>
	LinearArray& operator=(const ArrayB& o)
	{
		Array<VALUETYPE,D,StorageType>::operator=(o);
		return *this;
	}
};


//external values because SFNIAE for operator selection.
#define BINARY_OPERATOR_TEMPLATE( XXX ) \
template<class V,unsigned int D,template<unsigned int> class StG,class V2,template<unsigned int> class StG2> \
auto operator XXX(const LinearArray<V,D,StG>& a,const Array<V2,D,StG2>& b) \
	->typename impl::ElementWiseMeta<LinearArray,decltype(a[0] XXX b[0]),D,StG>::type \
{ \
	return a.elementWise([](const V& v1,const V2& v2){ return v1 XXX v2; },b); \
} \
template<class V,unsigned int D,template<unsigned int> class StG,class TypeB> \
auto operator XXX(const LinearArray<V,D,StG>& a,const TypeB& b) \
->typename impl::ElementWiseMeta<LinearArray,decltype(a[0] XXX b),D,StG>::type \
{ \
	return a.elementWise([&b](const V& v1){ return v1 XXX b; }); \
} \

/*
template<class V,unsigned int D,template<unsigned int> class StG,class TypeB> \
auto operator XXX(const TypeB& b,const LinearArray<V,D,StG>& a) \
->typename impl::ElementWiseMeta<LinearArray,decltype(b XXX a[0]),D,StG>::type \
{ \
	return a.elementWise([&b](const V& v1){ return b XXX v1; }); \
} \
*/


BINARY_OPERATOR_TEMPLATE(+)
BINARY_OPERATOR_TEMPLATE(-)
BINARY_OPERATOR_TEMPLATE(*)
BINARY_OPERATOR_TEMPLATE(/)
BINARY_OPERATOR_TEMPLATE(%)
/*


template<class V,unsigned int D,template<unsigned int> class StG,class ArrayB>
auto operator-(const LinearArray<V,D,StG>& a,const ArrayB& b)
->typename impl::ElementWiseMeta<LinearArray,decltype(a[0]-b[0]),D,StG>::type
{
	return a.elementWise([](const V& v1,const typename ArrayB::value_type& v2){ return v1-v2; },b);
}
template<class V,unsigned int D,template<unsigned int> class StG,class ArrayB>
auto operator/(const LinearArray<V,D,StG>& a,const ArrayB& b)
->typename impl::ElementWiseMeta<LinearArray,decltype(a[0]/b[0]),D,StG>::type
{
	return a.elementWise([](const V& v1,const typename ArrayB::value_type& v2){ return v1/v2; },b);
}
template<class V,unsigned int D,template<unsigned int> class StG,class ArrayB>
auto operator%(const LinearArray<V,D,StG>& a,const ArrayB& b)
->typename impl::ElementWiseMeta<LinearArray,decltype(a[0]%b[0]),D,StG>::type
{
	return a.elementWise([](const V& v1,const typename ArrayB::value_type& v2){ return v1%v2; },b);
}

template<class V,unsigned int D,template<unsigned int> class StG,class ArrayB>
auto operator/(const LinearArray<V,D,StG>& a,const ArrayB& b)
->typename impl::ElementWiseMeta<LinearArray,decltype(a[0]*b[0]),D,StG>::type
{
	return a.elementWise([](const typename V& v1,const typename ArrayB::value_type& v2){ return v1*v2; },b);
}
template<class V,unsigned int D,template<unsigned int> class StG>
auto operator+(const LinearArray<V,D,StG>& a)
->typename impl::ElementWiseMeta<LinearArray,decltype(+a[0]),D,StG>::type
{
	return a.elementWise([](const V& v1){ return +v1; });
}

template<class V,unsigned int D,template<unsigned int> class StG>
auto operator-(const LinearArray<V,D,StG>& a)
->typename impl::ElementWiseMeta<LinearArray,decltype(+a[0]),D,StG>::type
{
	return a.elementWise([](const V& v1){ return -v1; });
}*/


}

#endif