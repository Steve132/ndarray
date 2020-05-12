#ifndef ND_CIMG_DISPLAY_HPP
#define ND_CIMG_DISPLAY_HPP

#include<nd/ndarray.hpp>
#include<string>
#include<memory>

namespace cimg_library
{
	struct CImgDisplay;
}

namespace nd
{
namespace cimg
{

typedef std::shared_ptr<cimg_library::CImgDisplay> DisplayHandle;
	
template<class T>
void display(const nd::Array<T,3,nd::RowMajorOrder>& arr,const std::string& title="Image");

template<class T>
void display(const nd::Array<T,2,nd::RowMajorOrder>& arr,const std::string& title="Image");

template<class T>
nd::Coord<2> select_point(const nd::Array<T,3,nd::RowMajorOrder>& arr,const std::string& title="Image");

template<class T>
nd::Coord<2> select_point(const nd::Array<T,2,nd::RowMajorOrder>& arr,const std::string& title="Image");

template<class T>
std::array<nd::Coord<2>,2> select_segment(const nd::Array<T,3,nd::RowMajorOrder>& arr,const std::string& title="Image");

template<class T>
std::array<nd::Coord<2>,2> select_segment(const nd::Array<T,2,nd::RowMajorOrder>& arr,const std::string& title="Image");

}

}

#endif
