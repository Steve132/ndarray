

#ifndef NDARRAY_SRC_STORAGE_ORDER_BASE_HPP
#define NDARRAY_SRC_STORAGE_ORDER_BASE_HPP

#include<array>
#include<cstddef>

namespace nd
{
namespace impl
{

template<unsigned int D>
class StorageOrderBase
{
public:
	using shape_type=std::array<size_t,D>;
	using index_type=std::array<ptrdiff_t,D>;
protected:
	shape_type _shape;
	StorageOrderBase(const shape_type& tshape):_shape(tshape)
	{}
public:
	const shape_type& shape() const { return _shape; }
	size_t shape(unsigned int d) const { return _shape[d]; }
};

}
}

#endif
