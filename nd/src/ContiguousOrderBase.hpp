#ifndef NDARRAY_SRC_CONTIGUOUS_ORDER_BASE_HPP
#define NDARRAY_SRC_CONTIGUOUS_ORDER_BASE_HPP
#include "ContiguousMemView.hpp"
namespace nd
{
	namespace impl
	{
		template<class VALUETYPE>
		class ContiguousOrderBaseImpl
		{
		public:
			using value_type=VALUETYPE;
			using iterator_type=VALUETYPE*;
			using const_iterator_type=const VALUETYPE*;
		protected:
			impl::ContiguousMemView<value_type> _data;
			
		public:
			ContiguousOrderBaseImpl(const size_t& sz, VALUETYPE* ptr):
				_data(ptr,sz)
			{}
			
			ContiguousOrderBaseImpl(const size_t& sz,const VALUETYPE& fval=VALUETYPE()):
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
