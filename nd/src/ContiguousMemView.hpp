#ifndef NDARRAY_SRC_CONTIGUOUS_MEMVIEW_HPP
#define NDARRAY_SRC_CONTIGUOUS_MEMVIEW_HPP

#include<memory>
#include<algorithm>
#include<iostream>

namespace nd
{
namespace impl
{
template<class VALUETYPE>
class ContiguousMemView
{
private:
	bool _shared;
	void deleter()
	{
		if(_ptr && !_shared)
		{
			delete [] _ptr;
			_ptr=nullptr;
		}
		_ptr=nullptr;
		_size=0;
	}
	
	void copy(const ContiguousMemView& o)
	{
		if(!_ptr)
		{
			if(o._ptr) 
			{
				_ptr=new VALUETYPE[o._size];
				_size=o._size;
				std::copy(o._ptr,o._ptr+o._size,_ptr);
			}
		}
		else
		{
			if(!o._ptr) 
			{
				deleter();
				_shared=false;
			}
			else 
			{
				if(_size != o._size)
				{
					deleter();
					_ptr=new VALUETYPE[o._size];
					_size=o._size;
					_shared=false;
				}
				std::copy(o._ptr,o._ptr+o._size,_ptr);
			}
		}
	}
public:
	VALUETYPE* _ptr;
	size_t _size;
	ContiguousMemView(VALUETYPE* ptr,size_t sz=0):_ptr(sz==0 ? nullptr : ptr),_size(ptr ? sz : 0),_shared(true)
	{


	}
	ContiguousMemView(size_t sz=0,const VALUETYPE& v={}):_ptr(nullptr),_size(sz),_shared(false)
	{
		if(sz) 
		{
			_ptr=new VALUETYPE[sz];
			std::fill(_ptr,_ptr+sz,v);
		}
	}
	ContiguousMemView& operator=(ContiguousMemView&& o)
	{
		deleter();
		_ptr=o._ptr;_size=o._size;_shared=o._shared;
		o._ptr=nullptr;o._size=0;o._shared=false;
		return *this;
	}
	ContiguousMemView(ContiguousMemView&& o)
	{
		_ptr=o._ptr;_size=o._size;_shared=o._shared;
		o._ptr=nullptr;o._size=0;o._shared=false;
	}
	
	ContiguousMemView& operator=(const ContiguousMemView& o)
	{
		copy(o);
		return *this;
	}
	ContiguousMemView(const ContiguousMemView& o):_ptr(nullptr),_size(0),_shared(false)
	{
		copy(o);
	}
	~ContiguousMemView()
	{
		deleter();
	}
};

}
}

#endif
