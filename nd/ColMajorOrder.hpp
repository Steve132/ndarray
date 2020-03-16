#ifndef NDARRAY_COL_MAJOR_ORDER_HPP
#define NDARRAY_COL_MAJOR_ORDER_HPP

#include<numeric>
#include<iterator>
#include "src/StorageOrderBase.hpp"

namespace nd
{
struct ColMajorOrder
{
template<unsigned int D>
class Impl: public impl::StorageOrderBase<D>
{
public:
	using shape_type=typename impl::StorageOrderBase<D>::shape_type;
	using index_type=typename impl::StorageOrderBase<D>::index_type;
protected:
	shape_type _strides;
public:
	Impl(const shape_type& tshape={}):
		impl::StorageOrderBase<D>(tshape)
	{
		const shape_type& shp=impl::StorageOrderBase<D>::shape();
		std::partial_sum(
			std::begin(shp),
						 std::begin(shp)+D,
						 std::begin(_strides),
						 [](size_t a,size_t b){ return a*b; } 
		);
	}
	
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
	size_t size() const { return _strides[D-1]; }
	bool operator==(const Impl& o) const { return true; }
	template<class OrderImpl>
	bool operator==(const OrderImpl& o) const { return false; }
};
};

}

#endif
