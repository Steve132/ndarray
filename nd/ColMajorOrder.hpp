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
		_strides[0]=1;
		for(unsigned int di=0;di<(D-1);di++)
		{
			_strides[di+1]=_strides[di]*tshape[di];
		}
		return _strides[D-1]*tshape[D-1];
	}
public:
	Impl(VALUETYPE* ptr,const shape_type& tshape):
		impl::ContiguousOrderBaseImpl<VALUETYPE,D>(setShape(tshape),ptr,tshape)
	{}
	
	Impl(const shape_type& tshape={},const VALUETYPE& fval=VALUETYPE()):
		impl::ContiguousOrderBaseImpl<VALUETYPE,D>(setShape(tshape),tshape,fval)
	{}
	
	size_t ravel(const index_type& index) const {
		return std::inner_product(std::begin(index),std::begin(index)+D,std::begin(_strides),size_t(0));
	}
	
	size_t ravel_neighbor_offset(size_t rdex,size_t dim,index_type off) const
	{
		return rdex+off*_strides[dim];
	}
	
	template<class IndexHead1,class IndexHead2,class... IndexTail>
	size_t ravel(const IndexHead1& index1,const IndexHead2& index2,const IndexTail&... tail) const {
		return ravel(index_type{(ptrdiff_t)index1,(ptrdiff_t)index2,static_cast<ptrdiff_t>(tail)...});
	}
	index_type unravel(size_t index) const {
		index_type out;
		auto outiter=D-1;
		auto rend=0;
		for(unsigned di=D-1;di!=0;di--){	
			size_t st=_strides[di];
			out[di]=index/st;
			index=index%st;
		}
		out[0]=index;
		return out;
	}
	index_type unravel(const iterator_type& at) const { //helper method for iterators
		return unravel(at-impl::ContiguousOrderBaseImpl<VALUETYPE,D>::begin());
	}
	
	bool operator==(const Impl& o) const { return true; }
	template<class OrderImpl>
	bool operator==(const OrderImpl& o) const { return false; }
	
	using impl::ContiguousOrderBaseImpl<VALUETYPE,D>::operator[];
	
	/*TODO: This doesn't work for row major, even
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
	}*/
};
};

}

#endif
