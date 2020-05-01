#ifndef NDARRAY_ROW_MAJOR_ORDER_HPP
#define NDARRAY_ROW_MAJOR_ORDER_HPP

#include<numeric>
#include<iterator>
#include "src/LayoutBase.hpp"
#include<type_traits>

namespace nd
{

struct RowMajorOrder
{
	template<unsigned int D>
	class Layout: public impl::LayoutBase<D>
	{
	public:
		using shape_type=typename impl::LayoutBase<D>::shape_type;
		using coord_type=typename impl::LayoutBase<D>::coord_type;
	protected:
		shape_type _strides;
	private:
		void setShape(const shape_type& tshape)
		{
			_strides[D-1]=1;
			for(unsigned int di=(D-1);di>0;di--)
			{
				_strides[di-1]=_strides[di]*tshape[di];
			}
		}
	protected:
		size_t num_elements() const { return _strides[0]*impl::LayoutBase<D>::shape(0); }
	public:
		const shape_type& stride() const { return _strides; }
		const size_t& stride(unsigned int di) const { return _strides[di]; }
		
		Layout(const shape_type& tshape):
			impl::LayoutBase<D>(tshape)
		{
			setShape(tshape);
		}
		
		size_t ravel(const coord_type& index) const {
			return std::inner_product(std::begin(index),std::begin(index)+D,std::begin(_strides),size_t(0));
		}
		
		size_t ravel_neighbor_offset(size_t rdex,size_t dim,coord_type off) const
		{
			return rdex+off*_strides[dim];
		}
		
		template<class IndexHead1,class IndexHead2,class... IndexTail>
		typename std::enable_if<(2+sizeof...(IndexTail))==D,size_t>::type
		ravel(const IndexHead1& index1,const IndexHead2& index2,const IndexTail&... tail) const {
			return ravel(coord_type{(ptrdiff_t)index1,(ptrdiff_t)index2,static_cast<ptrdiff_t>(tail)...});
		}
		coord_type unravel(size_t index) const {
			coord_type out;
			auto outiter=D-1;
			auto rend=0;
			
			for(unsigned di=D-1;di>0;di--){	
				size_t st=_strides[(D-1)-di];
				out[di]=index/st;
				index=index%st;
			}
			out[0]=index;
			return out;
		}
	};
};

}


#endif
