#ifndef NDARRAY_ZORDER_HPP
#define NDARRAY_ZORDER_HPP

#include<numeric>
#include<iterator>
#include "src/ContiguousOrderBase.hpp"
#include<bitset>
#include<stdexcept>
#include<x86intrin.h>

namespace nd
{
struct ZOrder
{
template<class VALUETYPE,unsigned int D>
class Impl: public impl::ContiguousOrderBaseImpl<VALUETYPE,D>
{
public:
	using value_type=typename impl::ContiguousOrderBaseImpl<VALUETYPE,D>::value_type;
	using shape_type=typename impl::ContiguousOrderBaseImpl<VALUETYPE,D>::shape_type;
	using index_type=typename impl::ContiguousOrderBaseImpl<VALUETYPE,D>::index_type;
	using iterator_type=typename impl::ContiguousOrderBaseImpl<VALUETYPE,D>::iterator_type;
	using const_iterator_type=typename impl::ContiguousOrderBaseImpl<VALUETYPE,D>::const_iterator_type;
protected:
	shape_type _masks;
	static inline std::uint64_t pdep(std::uint64_t x,std::uint64_t msk)
	{
		return _pdep_u64(x,msk);
	}
	static inline std::uint32_t pdep(std::uint32_t x,std::uint32_t msk)
	{
		return _pdep_u32(x,msk);
	}
	static inline std::uint64_t pext(std::uint64_t x,std::uint64_t msk)
	{
		return _pext_u64(x,msk);
	}
	static inline std::uint32_t pext(std::uint32_t x,std::uint32_t msk)
	{
		return _pext_u32(x,msk);
	}
	static inline unsigned int popcount(size_t a)
	{
		return std::bitset<sizeof(size_t)*8>(a).count(); 
	}
	static inline size_t ctz(size_t x)
	{
		return popcount((x & -x)-1);
	}
	
	static constexpr size_t mask_compute(unsigned int startpoint,unsigned int maxpoint, unsigned int Delta)
	{
		size_t out=0;
		for(unsigned int ci=startpoint;ci<maxpoint;ci+=Delta)
		{
			out|=size_t(1) << ci;
		}
		return out;
	}
	
	size_t computeMasksHypercube(const shape_type& shp)
	{
		size_t sz=std::accumulate(std::begin(shp),std::begin(shp),size_t(1),[](size_t a,size_t b) { return a*b; });
		//confirm that all dimensions are the same power of two:  
		size_t ored_sizes=std::accumulate(
			std::begin(shp),std::begin(shp)+D,
			[](const std::size_t& s1,const std::size_t& s2) { return s1 | s2; },
			0
		);
		if(popcount(ored_sizes)!=1)
		{
			throw std::invalid_argument("A Z-Order array must (currently) be a power-of-two hypercube (all dimensions are the same power of 2)");
		}
		for(size_t d=0;d<D;d++)
		{
			_masks[d]=mask_compute(d,sizeof(size_t)*8,D) & (shp[0]-1);
		}
		return sz;
	}
	
	size_t computeMasksOuter(const shape_type& shp)
	{
		size_t sz=std::accumulate(std::begin(shp),std::begin(shp)+D,size_t(1),[](size_t a,size_t b) { return a*b; });
		
		if(!std::all_of(std::begin(shp),std::begin(shp)+D,[](size_t a){ return popcount(a) == 1; }))
		{
			throw std::invalid_argument("Each element of a Z-Order array must (currently) be a power-of-two");
		}
		
		
		shape_type sorted_shape=shp;
		std::sort(std::begin(sorted_shape),std::begin(sorted_shape)+D);
		for(unsigned int d=0;d<D;d++)
		{
			unsigned int ci=d;
			size_t outmask=0;
			for(unsigned nzc=0;sorted_shape[nzc]==0 && nzc < D;nzc++) ci--;
			for(size_t current_dim_bit=0;(size_t(1)<<current_dim_bit)<shp[d];current_dim_bit++)
			{
				unsigned int ciDid=0;
				while((size_t(1)<<current_dim_bit) >= sorted_shape[ciDid]) ciDid++; 
				ci+=(D-ciDid);
				outmask|=size_t(1) << ci;
			}
			_masks[d]=outmask >> D;
		}
		return sz;
	}
	size_t setMasks(const shape_type& shp)
	{
		return computeMasksOuter(shp);
	}
public:
	
	const shape_type& masks() const { return _masks; }
	
	Impl(const shape_type& tshape={},const value_type& v={}):
		impl::ContiguousOrderBaseImpl<VALUETYPE,D>(setMasks(tshape),tshape,v)
	{}
	
	template<class IndexType>
	size_t ravel(const IndexType& index) const {
		size_t out=0;
		for(unsigned int d=0;d<D;d++)
		{
			out|=pdep(index[d],_masks[d]);
		}
	}

	template<class IndexHead1,class IndexHead2,class... IndexTail>
	size_t ravel(const IndexHead1& index1,const IndexHead2& index2,const IndexTail&... tail) const {
		return ravel(std::array<IndexHead1,2+sizeof...(tail)>{index1,index2,tail...});
	}
	index_type unravel(size_t index) const {
		index_type out;
		for(unsigned int d=0;d<D;d++)
		{
			out[d]=pext(index,_masks[d]);
		}
		return out;
	}
	
	bool operator==(const Impl& o) const { return true; }
	template<class OrderImpl>
	bool operator==(const OrderImpl& o) const { return false; }
};
};

}

#endif
