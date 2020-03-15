#ifndef NDARRAY_HPP
#define NDARRAY_HPP

#include<array>
#include<numeric>
#include<functional>
#include<memory>
#include<type_traits>
#include<algorithm>
#include<iterator>
#include<iostream>

namespace nd
{
namespace impl
{
template<class VALUETYPE>
class contiguous_memview
{
	std::unique_ptr<VALUETYPE[]> _vecdata;
public:
	VALUETYPE* _ptr;
	size_t _size;
	contiguous_memview(VALUETYPE* ptr,size_t sz=0):_ptr(ptr),_size(sz)
	{}
	contiguous_memview(size_t sz=0,const VALUETYPE& v={}):_ptr(nullptr),_size(sz)
	{
		if(sz) 
		{
			_vecdata.reset(new VALUETYPE[sz]);
			std::fill(_vecdata.get(),_vecdata.get()+sz,v);
			_ptr=_vecdata.get();
		}
	}
	contiguous_memview& operator=(contiguous_memview&& o)
	{
		_vecdata=std::move(o._vecdata);
		_ptr=_vecdata ? _vecdata.get() : o._ptr;_size=o._size;
		return *this;
	}
	contiguous_memview(contiguous_memview&& o): _vecdata(std::move(o._vecdata)),
		_ptr(_vecdata ? _vecdata.get() : o._ptr),
		_size(o._size)
	{}
	contiguous_memview& operator=(const contiguous_memview& o)
	{
		if(!o._vecdata)
		{
			_vecdata.release();
			_ptr=o._ptr;_size=o._size;
		}
		else
		{
			if(_size != o._size || !_vecdata)
			{
				_vecdata.reset(new VALUETYPE[o._size]);
			}
			std::copy(o._vecdata.get(),o._vecdata.get()+o._size,_vecdata.get());
			_ptr=_vecdata.get();_size=o._size;
		}
		return *this;
	}
	contiguous_memview(const contiguous_memview& o)
	{
		operator=(o);
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
class Impl: public impl::StorageOrderBase<D>
{
public:
	using typename impl::StorageOrderBase<D>::size_type;
protected:
	size_type _strides;
public:
	Impl(const size_type& tshape={}):
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
	bool operator==(const Impl& o) const { return true; }
	template<class OrderImpl>
	bool operator==(const OrderImpl& o) const { return false; }
};
};

template<class VALUETYPE,unsigned int D,class AbstractStorageType=ColMajorOrder >
class Array
{
public:
	using value_type=VALUETYPE;
	using storage_type=typename AbstractStorageType::template Impl<D>;
	using size_type=typename storage_type::size_type;
	using iterator_type=VALUETYPE*;
	static const unsigned int num_dimensions=D;
	using abstract_storage_type=AbstractStorageType; // type-id is vector<T, Alloc<T>>
	
protected:
	storage_type _storage;
	impl::contiguous_memview<value_type> _data;

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
	
	
	iterator_type begin() { return data(); }
	const iterator_type begin() const { return _data._ptr; }
	
	iterator_type end() { return data()+_data._size; }
	const iterator_type end() const { return data()+_data._size; }
	
	
	const value_type& operator[](size_t dex) const { return _data._ptr[dex]; }
	value_type& operator[](size_t dex) { return _data._ptr[dex]; }

	
	value_type* data() { return _data._ptr; }
	const value_type* data() const { return _data._ptr; }
	
	
	template<typename... Args>
	size_t linearize(Args&&... args) const {
		return _storage.linearize(std::forward<Args>(args)...);
	}
	
	size_type unlinearize(size_t dex) const {
		return _storage.unlinearize(dex);
	}
	size_type unlinearize(const iterator_type& at) const { //helper method for iterators
		return _storage.unlinearize(at-begin());
	}
	template<typename... Args>
	auto shape(Args&&... args) const -> decltype(_storage.shape(std::forward<Args>(args)...)) {
		return _storage.shape(std::forward<Args>(args)...);
	}
	
	value_type* data(const size_type& offset) { return _data._ptr+linearize(offset); }
	const value_type* data(const size_type& offset) const { return _data._ptr+linearize(offset); }

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
				size_t ai=linearize(b.unlinearize(bi));
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
			for(size_t bi=0;bi<bN;bi++)
			{
				size_t ai=linearize(b.unlinearize(bi));
				f(operator[](ai),b[bi]);
			}
		}
		else
		{
			size_t bN=b.size();
			auto aptr=data();auto bptr=b.data();
			for(size_t i=0;i<bN;i++)
			{
				f(aptr[i],bptr[i]);
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
};

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
template<class ArrayClass,typename std::enable_if<(ArrayClass::num_dimensions > 2),int>::type = 0>
std::ostream& operator<<(std::ostream& out,const ArrayClass& arr)
{
	for(size_t i=0;i<arr.shape(ArrayClass::num_dimensions-1);i++)
	{
		
		Array<typename ArrayClass::value_type,ArrayClass::num_dimensions-1,typename ArrayClass::abstract_storage_type> outmap;
	}
	return out;
}

//external values because SFNIAE for operator selection.
#define BINARY_OPERATOR_TEMPLATE( XXX ) \
template<class V,unsigned int D,class StG,class V2,class StG2> \
auto operator XXX(const Array<V,D,StG>& a,const Array<V2,D,StG2>& b) \
	->Array<decltype(a[0] XXX b[0]),D,StG> \
{ \
	return a.elementWise([](const V& v1,const V2& v2){ return v1 XXX v2; },b); \
} \
template<class V,unsigned int D,class StG,class TypeB> \
auto operator XXX(const Array<V,D,StG>& a,const TypeB& b) \
	-> Array<decltype(a[0] XXX b),D,StG> \
{ \
	return a.elementWise([&b](const V& v1){ return v1 XXX b; }); \
} \
template<class V,unsigned int D,class StG,class TypeB> \
auto operator XXX(const TypeB& b,const Array<V,D,StG>& a) \
-> Array<decltype(b XXX a[0]),D,StG> \
{ \
	return a.elementWise([&b](const V& v1){ return b XXX v1; }); \
} \

BINARY_OPERATOR_TEMPLATE(*)
BINARY_OPERATOR_TEMPLATE(/)
BINARY_OPERATOR_TEMPLATE(%)
BINARY_OPERATOR_TEMPLATE(+)
BINARY_OPERATOR_TEMPLATE(-)
BINARY_OPERATOR_TEMPLATE(<<)
BINARY_OPERATOR_TEMPLATE(>>)
BINARY_OPERATOR_TEMPLATE(<)
BINARY_OPERATOR_TEMPLATE(>)
BINARY_OPERATOR_TEMPLATE(<=)
BINARY_OPERATOR_TEMPLATE(>=)
BINARY_OPERATOR_TEMPLATE(==)
BINARY_OPERATOR_TEMPLATE(!=)
BINARY_OPERATOR_TEMPLATE(&)
BINARY_OPERATOR_TEMPLATE(|)
BINARY_OPERATOR_TEMPLATE(^)
BINARY_OPERATOR_TEMPLATE(&&)
BINARY_OPERATOR_TEMPLATE(||)

#define UNARY_OPERATOR_TEMPLATE( XXX ) \
template<class V,unsigned int D, class StG> \
auto operator XXX(const Array<V,D,StG>& a) \
->Array<decltype(XXX(a[0])),D,StG> \
{ \
	return a.elementWise([](const V& v1){ return XXX(v1); }); \
} \

UNARY_OPERATOR_TEMPLATE(~)
UNARY_OPERATOR_TEMPLATE(!)
UNARY_OPERATOR_TEMPLATE(-)
UNARY_OPERATOR_TEMPLATE(+)

#define INPLACE_OPERATOR_TEMPLATE( XXX ) \
template<class V,unsigned int D,class StG,class V2, class StG2> \
auto operator XXX(Array<V,D,StG>& a,const Array<V2,D,StG2>& b) \
->Array<V,D,StG>& \
{ \
	return a.elementWiseInPlace([](V& v1,const V2& v2){ v1 XXX v2; },b); \
} \
template<class V,unsigned int D,class StG,class TypeB> \
auto operator XXX(Array<V,D,StG>& a,const TypeB& b) \
-> Array<V,D,StG>& \
{ \
	return a.elementWiseInPlace([&b](V& v1){ v1 XXX b; }); \
} \

INPLACE_OPERATOR_TEMPLATE(*=)
INPLACE_OPERATOR_TEMPLATE(/=)
INPLACE_OPERATOR_TEMPLATE(+=)
INPLACE_OPERATOR_TEMPLATE(-=)
INPLACE_OPERATOR_TEMPLATE(>>=)
INPLACE_OPERATOR_TEMPLATE(<<=)
INPLACE_OPERATOR_TEMPLATE(&=)
INPLACE_OPERATOR_TEMPLATE(^=)
INPLACE_OPERATOR_TEMPLATE(|=)
//++,--










}

#endif
