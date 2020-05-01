#ifndef NDARRAY_TRIANGLEORDER_HPP
#define NDARRAY_TRIANGLEORDER_HPP

#include<numeric>
#include<iterator>
#include "src/ContiguousOrderBase.hpp"

#include<stdexcept>


namespace nd
{
	struct TriangleOrder
	{
		template<class VALUETYPE,unsigned int D>
		class Impl: public impl::ContiguousOrderBaseImpl<VALUETYPE,D>
		{
		public:
			using value_type=typename impl::ContiguousOrderBaseImpl<VALUETYPE,D>::value_type;
			using shape_type=typename impl::ContiguousOrderBaseImpl<VALUETYPE,D>::shape_type;
			using coord_type=typename impl::ContiguousOrderBaseImpl<VALUETYPE,D>::coord_type;
			using iterator_type=typename impl::ContiguousOrderBaseImpl<VALUETYPE,D>::iterator_type;
			using const_iterator_type=typename impl::ContiguousOrderBaseImpl<VALUETYPE,D>::const_iterator_type;
		protected:
			shape_type tshape(size_t n)
			{
				shape_type shpout;
				for(unsigned int d=0;d<D;d++) shpout[d]=n;
				return shpout;
			}
			
			constexpr static std::array<size_t,D> get_denoms()
			{
				std::array<size_t,D> out;
				out[0]=1;
				for(unsigned int i=1;i<D;i++)
				{
					out[i]=out[i-1]*(i+1);
				}
				return out;
			}
			constexpr static std::array<size_t,D> _denoms=get_denoms();
			size_t tsize(size_t n)
			{
				size_t nfact=1;
				size_t dfact=1;
				for(unsigned int d=0;d<D;d++)
				{
					nfact*=(n+d);
					dfact*=(1+d);
				}
				return nfact/_denoms[D];
			}
		public:
			Impl(size_t n=0,const value_type& v={}):
				impl::ContiguousOrderBaseImpl<VALUETYPE,D>(tsize(n),tshape(n),v)
			{}
			Impl(VALUETYPE* ptr,size_t n):
				impl::ContiguousOrderBaseImpl<VALUETYPE,D>(tsize(n),ptr,tshape(n))
			{}
			
			template<class IndexType>
			size_t ravel(const IndexType& index) const {
				size_t out=0;
				for(unsigned int d=0;d<D;d++)
				{
					size_t sval=index[d];
					size_t x=sval;
					for(unsigned int i=0;i<d;i++)
					{
						sval*=++x;
					}
					out+=sval/_denoms[d];
				}
				return out;
			}
			
			template<class IndexHead1,class IndexHead2,class... IndexTail>
			size_t ravel(const IndexHead1& index1,const IndexHead2& index2,const IndexTail&... tail) const {
				return ravel(std::array<IndexHead1,2+sizeof...(tail)>{index1,index2,tail...});
			}
			coord_type unravel(size_t index) const {
				coord_type out;
				/*for(unsigned int d=0;d<D;d++)
				{
					out[d]=pext(index,_masks[d]);
				}*/
				return out;
			}
			
			size_t ravel_neighbor_offset(size_t rdex,size_t dim,coord_type off) const
			{
				size_t out;
				return out;
			}
			
			bool operator==(const Impl& o) const { return true; }
			template<class OrderImpl>
			bool operator==(const OrderImpl& o) const { return false; }
		};
	};
	
}

#endif

