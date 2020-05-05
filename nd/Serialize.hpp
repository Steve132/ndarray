#ifndef ND_SERIALIZE_HPP
#define ND_SERIALIZE_HPP

#include<iostream>
#include<type_traits>
#include<algorithm>


namespace nd
{

	template<class ArrayType,class PutFunctionType,
		typename std::enable_if<
			std::is_fundamental<typename ArrayType::value_type>::value, 
		int>::type = 0
	>
	void write(const ArrayType& arr,PutFunctionType put);
};

#if __cplusplus >= 202002L
#include<bit>
namespace nd { namespace impl {
	static constexpr bool is_big_endian=std::endian::native==std::endian::big;
}
}
#else
namespace nd
{
namespace impl
{
#if defined(__BYTE_ORDER) && __BYTE_ORDER == __BIG_ENDIAN || \
defined(__BIG_ENDIAN__) || \
defined(__ARMEB__) || \
defined(__THUMBEB__) || \
defined(__AARCH64EB__) || \
defined(_MIBSEB) || defined(__MIBSEB) || defined(__MIBSEB__)
static constexpr bool is_big_endian=true;
#elif defined(__BYTE_ORDER) && __BYTE_ORDER == __LITTLE_ENDIAN || \
defined(__LITTLE_ENDIAN__) || \
defined(__ARMEL__) || \
defined(__THUMBEL__) || \
defined(__AARCH64EL__) || \
defined(_MIPSEL) || defined(__MIPSEL) || defined(__MIPSEL__)
static constexpr bool is_big_endian=false;
#else
#error "I don't know what architecture this is!"
#endif
static constexpr uint32_t serialize_magic_number=0x0007ab1e;
}
}
#endif

namespace nd
{

namespace impl
{
	template<class T,class PutFunctionType> 
	inline void endian_put(const T& a,PutFunctionType& pf)
	{
		const uint8_t* fb=reinterpret_cast<const uint8_t*>(&a);
		if(is_big_endian)
		{
			for(unsigned i=0;i<sizeof(T);i++)
			{
				pf(fb+sizeof(T)-1-i,fb+sizeof(T)-i);
			}
		}
		else
		{
			pf(fb,fb+sizeof(T));
		}
	}
	
	template<class T,class PutFunctionType>
	inline void endian_put(const T* be,const T* ed,PutFunctionType& pf)
	{
		if(is_big_endian)
		{
			for(;be!=ed;be++) endian_put(*be,pf);
		}
		else
		{
			pf(be,ed);
		}
	}
	template<class T,class GetFunctionType> 
	inline void endian_get(T& a,GetFunctionType& pf)
	{
		uint8_t* fb=reinterpret_cast<uint8_t*>(&a);
		if(is_big_endian)
		{
			for(unsigned i=0;i<sizeof(T);i++)
			{
				pf(fb+sizeof(T)-1-i,fb+sizeof(T)-i);
			}
		}
		else
		{
			pf(fb,fb+sizeof(T));
		}
	}
	
	template<class T,class GetFunctionType>
	inline void endian_get(T* be,T* ed,GetFunctionType& pf)
	{
		if(is_big_endian)
		{
			for(;be!=ed;be++) endian_get(*be,pf);
		}
		else
		{
			pf(be,ed);
		}
	}
	/*template<class T>
	static inline void endian_get(const T* be,const T* end)
	{
		if(is_big_endian)
		{
			uint8_t* fb=reinterpret_cast<uint8_t*>(&a);
			uint8_t* lb=fb+sizeof(T);
			std::reverse(fb,lb);
		}
	}*/
	template<class T>
	uint32_t get_serialize_typecode() { return 0x80000000+sizeof(T); }
	
	template<> inline uint32_t get_serialize_typecode<int8_t>() { return 1; }
	template<> inline uint32_t get_serialize_typecode<uint8_t>() { return 2; }
	template<> inline uint32_t get_serialize_typecode<int16_t>() { return 3; }
	template<> inline uint32_t get_serialize_typecode<uint16_t>() { return 4; }
	template<> inline uint32_t get_serialize_typecode<int32_t>() { return 5; }
	template<> inline uint32_t get_serialize_typecode<uint32_t>() { return 6; }
	template<> inline uint32_t get_serialize_typecode<int64_t>() { return 7; }
	template<> inline uint32_t get_serialize_typecode<uint64_t>() { return 8; }
	template<> inline uint32_t get_serialize_typecode<float>() { return 9; }
	template<> inline uint32_t get_serialize_typecode<double>() { return 10; }
	template<> inline uint32_t get_serialize_typecode<long double>() { return 11; }
	
	template<class T>
	uint32_t get_storage_typecode() { return 0x80000000; }
	
	template<> inline uint32_t get_storage_typecode<nd::RowMajorOrder>(){ return 1;}
	//template<> inline uint32_t get_storage_typecode<nd::RowMajorOrder>(){ return 2;}
}

template<class ArrayType,class PutFunctionType,
typename std::enable_if<
std::is_fundamental<typename ArrayType::value_type>::value, 
int>::type = 0
>
void write(const ArrayType& arr,PutFunctionType& put)
{
	endian_put(impl::serialize_magic_number,put);
	endian_put(impl::get_serialize_typecode<typename ArrayType::value_type>(),put);
	endian_put(static_cast<uint16_t>(ArrayType::num_dimensions),put);
	for(size_t d:arr.shape())
	{
		endian_put(static_cast<uint64_t>(d),put);
	}
	endian_put(impl::get_storage_typecode<typename ArrayType::abstract_storage_type>(),put);
	endian_put(std::begin(arr),std::end(arr),put);
}
template<class ArrayType,class PutFunctionType,
typename std::enable_if<
std::is_fundamental<typename ArrayType::value_type>::value, 
int>::type = 0
>
ArrayType read(PutFunctionType& getf)
{
	uint32_t magicnum;
	endian_get(magicnum,getf);
	if(impl::serialize_magic_number != magicnum) throw std::runtime_error("Bad magic number!");
	
	uint32_t typecode;
	endian_get(typecode,getf);
	uint32_t etypecode=impl::get_serialize_typecode<typename ArrayType::value_type>();
	if(etypecode != typecode) throw std::runtime_error("Did not recognize typecode");
	
	uint16_t dims;
	endian_get(dims,getf);
	unsigned int edims=ArrayType::num_dimensions;
	if(edims != dims) throw std::runtime_error("Number of dimensions does not match");
	
	nd::Shape<ArrayType::num_dimensions> shapein;
	for(size_t di=0;di<ArrayType::num_dimensions;di++)
	{
		uint64_t tdim;
		endian_get(tdim,getf);
		shapein[di]=static_cast<size_t>(tdim);
	}
	
	uint32_t abstype;
	endian_get(abstype,getf);
	uint32_t eabstype=impl::get_storage_typecode<typename ArrayType::abstract_storage_type>();
	if(abstype != eabstype) throw std::runtime_error("Bad abstract_storage_type");
	
	ArrayType arr(shapein);
	
	endian_get(std::begin(arr),std::end(arr),getf);
	return arr;
}
}
#endif
