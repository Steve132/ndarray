#ifndef NDARRAY_SRC_LAYOUT_BASE_HPP
#define NDARRAY_SRC_LAYOUT_BASE_HPP

#include<array>
#include<cstddef>

namespace nd
{
namespace impl
{

template<unsigned int D>
class LayoutBase
{
public:
	using index_type=ptrdiff_t;
	using shape_type=std::array<size_t,D>;
	using coord_type=std::array<index_type,D>;
protected:
	shape_type _shape;
	LayoutBase(const shape_type& tshape):_shape(tshape)
	{}
	size_t num_elements() const;
public:
	const shape_type& shape() const { return _shape; }
	size_t shape(unsigned int d) const { return _shape[d]; }
};

}
}

#endif
