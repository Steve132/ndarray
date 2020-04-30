#ifndef NDARRAY_SRC_CONTIGUOUS_MEMVIEW_HPP
#define NDARRAY_SRC_CONTIGUOUS_MEMVIEW_HPP

#include<memory>
#include<algorithm>

namespace nd
{
namespace impl
{
template<class VALUETYPE>
class ContiguousMemView
{
	std::unique_ptr<VALUETYPE[]> _vecdata;
public:
	VALUETYPE* _ptr;
	size_t _size;
	ContiguousMemView(VALUETYPE* ptr,size_t sz=0):_ptr(ptr),_size(sz)
	{}
	ContiguousMemView(size_t sz=0,const VALUETYPE& v={}):_ptr(nullptr),_size(sz)
	{
		if(sz) 
		{
			_vecdata.reset(new VALUETYPE[sz]);
			std::fill(_vecdata.get(),_vecdata.get()+sz,v);
			_ptr=_vecdata.get();
		}
	}
	ContiguousMemView& operator=(ContiguousMemView&& o)
	{
		_vecdata=std::move(o._vecdata);
		_ptr=_vecdata ? _vecdata.get() : o._ptr;_size=o._size;
		return *this;
	}
	ContiguousMemView(ContiguousMemView&& o): _vecdata(std::move(o._vecdata)),
		_ptr(_vecdata ? _vecdata.get() : o._ptr),
		_size(o._size)
	{}
	
	ContiguousMemView& operator=(const ContiguousMemView& o)
	{
		if(!o._vecdata)
		{
			_vecdata.release();
			_ptr=o._ptr;_size=o._size;
		}
		else
		{
			if(_size != o._size || !_vecdata)
			{
				_vecdata.reset(new VALUETYPE[o._size]);
			}
			std::copy(o._vecdata.get(),o._vecdata.get()+o._size,_vecdata.get());
			_ptr=_vecdata.get();_size=o._size;
		}
		return *this;
	}
	ContiguousMemView(const ContiguousMemView& o)
	{
		operator=(o);
	}
};

}
}

#endif
