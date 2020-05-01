#ifndef NDARRAY_SRC_CONTIGUOUS_ORDER_BASE_HPP
#define NDARRAY_SRC_CONTIGUOUS_ORDER_BASE_HPP

#include "StorageOrderBase.hpp"

namespace nd
{
	namespace impl
	{
		template<class VALUETYPE,unsigned int D>
		class ContiguousOrderBaseImpl: public impl::StorageOrderBase<D>
		{
		public:
			using value_type=VALUETYPE;
			using shape_type=typename impl::StorageOrderBase<D>::shape_type;
			using index_type=typename impl::StorageOrderBase<D>::index_type;
			using iterator_type=VALUETYPE*;
			using const_iterator_type=const VALUETYPE*;
		protected:
			impl::ContiguousMemView<value_type> _data;
			
		public:
			ContiguousOrderBaseImpl(const size_t& sz, VALUETYPE* ptr,const shape_type& tshape):
			impl::StorageOrderBase<D>(tshape),
			_data(ptr,sz)
			{}
			
			ContiguousOrderBaseImpl(const size_t& sz,const shape_type& tshape={},const VALUETYPE& fval=VALUETYPE()):
			impl::StorageOrderBase<D>(tshape),
			_data(sz,fval)
			{}
			size_t size() const { return _data._size; }
			const value_type& operator[](size_t dex) const { return _data._ptr[dex]; }
			value_type& operator[](size_t dex) { return _data._ptr[dex]; }
			
			iterator_type begin() { return data(); }
			const_iterator_type begin() const { return _data._ptr; }
			
			iterator_type end() { return data()+_data._size; }
			const_iterator_type end() const { return data()+_data._size; }
			
			value_type* data() { return _data._ptr; }
			const value_type* data() const { return _data._ptr; }
		};		
	}
}

#endif
