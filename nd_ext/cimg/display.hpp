#ifndef ND_CIMG_DISPLAY_HPP
#define ND_CIMG_DISPLAY_HPP

#include<nd/ndarray.hpp>
#include<string>

namespace nd
{
namespace cimg
{

template<class T>
std::array<unsigned int,3> display(const nd::Array<T,3,nd::RowMajorOrder>& arr,const std::string& title="Image");

template<class T>
std::array<unsigned int,3> display(const nd::Array<T,2,nd::RowMajorOrder>& arr,const std::string& title="Image");

template<class T>
nd::Coord<2> select_point(const nd::Array<T,3,nd::RowMajorOrder>& arr,const std::string& title="Image");

template<class T>
nd::Coord<2> select_point(const nd::Array<T,2,nd::RowMajorOrder>& arr,const std::string& title="Image");

}

}

#endif
