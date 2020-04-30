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
		
	
		std::array<unsigned int,D> sorted_shape_indices;
		std::iota(std::begin(sorted_shape_indices),std::end(sorted_shape_indices),0);
		std::sort(std::begin(sorted_shape_indices),std::end(sorted_shape_indices),
			[&shp](unsigned int ai,unsigned int bi) { return shp[ai] < shp[bi]; }
		);
		
		shape_type sorted_shape=shp;
		for(unsigned int ai=0;ai<D;ai++) sorted_shape[ai]=shp[sorted_shape_indices[ai]];
		//std::sort(std::begin(sorted_shape),std::begin(sorted_shape)+D);
		unsigned int cistart=0;
		size_t cisize=1;
		std::fill(std::begin(_masks),std::begin(_masks)+D,0);
		
		for(unsigned int sdi=0;sdi<D;sdi++)
		{
			if(cisize != sorted_shape[sdi])
			{
				size_t maxci=0;
				for(unsigned int sdi2=sdi;sdi2<D;sdi2++)
				{
					unsigned int d=sorted_shape_indices[sdi2];
					size_t maskout=0;
					unsigned int ci=cistart+sdi2-sdi;
					for(size_t current_dim_bit=0;(cisize<<current_dim_bit)<sorted_shape[sdi];current_dim_bit++)
					{
						if(ci > maxci) maxci=ci;
						maskout|=size_t(1) << ci;
						ci+=D-sdi;
					}
					_masks[d]|=maskout;
				}
				cisize=sorted_shape[sdi];
				cistart=maxci+1;
			}
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
	
	size_t ravel_neighbor_offset(size_t rdex,size_t dim,index_type off) const
	{
		size_t msk=_masks[dim];
		size_t nmsk=~_masks[dim];
		size_t result=(rdex & msk)|nmsk+pdep(off,msk);
		return rdex & nmsk + result & msk;
	}
	
	bool operator==(const Impl& o) const { return true; }
	template<class OrderImpl>
	bool operator==(const OrderImpl& o) const { return false; }
};
};

}

#endif
