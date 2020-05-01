#ifndef ND_SHAPE_HPP
#define ND_SHAPE_HPP

#include<cstdlib>
#include<array>

namespace nd
{
template<unsigned D>
using Shape=std::array<size_t,D>;

template<unsigned D>
using Coord=std::array<ptrdiff_t,D>;
}

#endif
