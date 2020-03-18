#ifndef NDARRAY_COL_MAJOR_ORDER_HPP
#define NDARRAY_COL_MAJOR_ORDER_HPP

#include<numeric>
#include<iterator>
#include "src/ContiguousOrderBase.hpp"

namespace nd
{
template<class V,unsigned int D, class StG> class Array;

struct ColMajorOrder
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
	shape_type _strides;
private:
	size_t setShape(const shape_type& tshape)
	{
		std::partial_sum(
			std::begin(tshape),
			std::begin(tshape)+D,
			std::begin(_strides),
			[](size_t a,size_t b){ return a*b; } 
		);
		return _strides[D-1];
	}
public:
	Impl(VALUETYPE* ptr,const shape_type& tshape):
		impl::ContiguousOrderBaseImpl<VALUETYPE,D>(setShape(tshape),ptr,tshape)
	{}
	
	Impl(const shape_type& tshape={},const VALUETYPE& fval=VALUETYPE()):
		impl::ContiguousOrderBaseImpl<VALUETYPE,D>(setShape(tshape),tshape,fval)
	{}
	
	template<class IndexType>
	size_t ravel(const IndexType& index) const {
		return index[0]+std::inner_product(std::begin(index)+1,std::begin(index)+D,std::begin(_strides),size_t(0));
	}
	template<class IndexHead1,class IndexHead2,class... IndexTail>
	size_t ravel(const IndexHead1& index1,const IndexHead2& index2,const IndexTail&... tail) const {
		return ravel(std::array<IndexHead1,2+sizeof...(tail)>{index1,index2,tail...});
	}
	index_type unravel(size_t index) const {
		index_type out;
		auto outiter=std::rbegin(out);
		auto rend=std::rbegin(_strides)+D;
		for(auto stiter=std::rbegin(_strides)+1;stiter!=rend;stiter++,outiter++){	
			*(outiter)=index/(*stiter);
			index=index%(*stiter);
		}
		*outiter=index;
		return out;
	}
	index_type unravel(const iterator_type& at) const { //helper method for iterators
		return unravel(at-impl::ContiguousOrderBaseImpl<VALUETYPE,D>::begin());
	}
	
	bool operator==(const Impl& o) const { return true; }
	template<class OrderImpl>
	bool operator==(const OrderImpl& o) const { return false; }
	
	using impl::ContiguousOrderBaseImpl<VALUETYPE,D>::operator[];
	
	//TODO make sure there's handling for D==0 or D==1
	//TODO: this doesnt' work with any data ordering that's not Rowise or Colwise (Zorder and non-contig don't work)
	//Honestly this is why we should go back to inheritance from storage order.  Different accessors for different storage orders
	//this also allows SNFIAE
	Array<value_type,D-1,ColMajorOrder> pop_order_ref(unsigned int d) const
	{
		using Am1=Array<value_type,D-1,ColMajorOrder>;
		typename Am1::shape_type shp;
		const shape_type& mshp= impl::StorageOrderBase<D>::shape();
		std::copy(std::begin(mshp),std::begin(mshp)+D-1,std::begin(shp));
		typename Am1::shape_type zed={};
		zed[D-1]=d;
		const value_type* slicebegin=impl::ContiguousOrderBaseImpl<VALUETYPE,D>::data()+ravel(zed);
		return Array<value_type,D-1,ColMajorOrder>(const_cast<value_type*>(slicebegin),shp);
	}
};
};

}

#endif
