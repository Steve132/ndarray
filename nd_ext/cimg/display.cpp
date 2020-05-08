#include "display.hpp"
#include <CImg.h>
using namespace cimg_library;

namespace nd {
namespace cimg {
//template<class T>    const CImg<T>& display(const char *const title=0, const bool display_info=true, unsigned int *const XYZ=0,
//const bool exit_on_anykey=false) const {
template<class T>
std::array<unsigned int,3> display(const nd::Array<T,3,nd::RowMajorOrder>& arr,const std::string& title)
{
	CImg<T> ref(arr.data(),arr.shape(2),arr.shape(1),arr.shape(0),1,true);
	std::array<unsigned int,3> arrout;
	ref.get_permute_axes("yzcx").display(title.c_str(),true,&arrout[0]);
	return arrout;
}

template<class T>
std::array<unsigned int,3> display(const nd::Array<T,2,nd::RowMajorOrder>& arr,const std::string& title)
{
	CImg<T> ref(arr.data(),1,arr.shape(1),arr.shape(0),1,true);
	std::array<unsigned int,3> arrout;
	ref.get_permute_axes("yzcx").display(title.c_str(),true,&arrout[0]);
	return arrout;
}

template std::array<unsigned int,3> display(const nd::Array<uint8_t,2,nd::RowMajorOrder>& arr,const std::string& title); 
template std::array<unsigned int,3> display(const nd::Array<uint8_t,3,nd::RowMajorOrder>& arr,const std::string& title);

template std::array<unsigned int,3> display(const nd::Array<float,2,nd::RowMajorOrder>& arr,const std::string& title); 
template std::array<unsigned int,3> display(const nd::Array<float,3,nd::RowMajorOrder>& arr,const std::string& title);

template std::array<unsigned int,3> display(const nd::Array<uint32_t,2,nd::RowMajorOrder>& arr,const std::string& title); 
template std::array<unsigned int,3> display(const nd::Array<uint32_t,3,nd::RowMajorOrder>& arr,const std::string& title);

template<class T>
nd::Coord<2> select_point(const nd::Array<T,3,nd::RowMajorOrder>& arr,const std::string& title)
{
	CImg<T> ref(arr.data(),arr.shape(2),arr.shape(1),arr.shape(0),1,true);
	CImg<uint32_t> result=ref.get_permute_axes("yzcx").get_select(title.c_str(),0);
	nd::Coord<2> out{result[1],result[0]};
	return out;
}
template<class T>
nd::Coord<2> select_point(const nd::Array<T,2,nd::RowMajorOrder>& arr,const std::string& title)
{
	CImg<T> ref(arr.data(),1,arr.shape(1),arr.shape(0),1,true);
	CImg<uint32_t> result=ref.get_permute_axes("yzcx").get_select(title.c_str(),0);
	nd::Coord<2> out{result[1],result[0]};
	return out;
}

template nd::Coord<2> select_point(const nd::Array<uint8_t,3,nd::RowMajorOrder>& arr,const std::string& title);
template nd::Coord<2> select_point(const nd::Array<uint8_t,2,nd::RowMajorOrder>& arr,const std::string& title);


template<class T>
std::array<nd::Coord<2>,2> select_segment(const nd::Array<T,3,nd::RowMajorOrder>& arr,const std::string& title)
{
	CImg<T> ref(arr.data(),arr.shape(2),arr.shape(1),arr.shape(0),1,true);
	CImg<uint32_t> result=ref.get_permute_axes("yzcx").get_select(title.c_str(),1);
	std::array<nd::Coord<2>,2> out{nd::Coord<2>{result[1],result[0]},nd::Coord<2>{result[4],result[3]}};  
	return out;
}
template<class T>
std::array<nd::Coord<2>,2> select_segment(const nd::Array<T,2,nd::RowMajorOrder>& arr,const std::string& title)
{
	CImg<T> ref(arr.data(),1,arr.shape(1),arr.shape(0),1,true);
	CImg<uint32_t> result=ref.get_permute_axes("yzcx").get_select(title.c_str(),1);
	std::array<nd::Coord<2>,2> out{nd::Coord<2>{result[1],result[0]},nd::Coord<2>{result[4],result[3]}};
	return out;
}

template std::array<nd::Coord<2>,2>  select_segment(const nd::Array<uint8_t,3,nd::RowMajorOrder>& arr,const std::string& title);
template std::array<nd::Coord<2>,2>  select_segment(const nd::Array<uint8_t,2,nd::RowMajorOrder>& arr,const std::string& title);

}
}
