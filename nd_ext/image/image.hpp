#ifndef ND_IMAGE_HPP
#define ND_IMAGE_HPP

#include<nd/ndarray.hpp>
#include<string>

namespace nd
{
namespace image
{

template<class T>
Array<T,3,RowMajorOrder> load_image(const std::string& filename,unsigned int desired_channels=0);
template<> Array<uint8_t,3,RowMajorOrder> load_image(const std::string& filename,unsigned int desired_channels);
template<> Array<float,3,RowMajorOrder> load_image(const std::string& filename,unsigned int desired_channels);
template<> Array<uint16_t,3,RowMajorOrder> load_image(const std::string& filename,unsigned int desired_channels);
template<class T>
Array<T,3,RowMajorOrder> load_image(const std::string& filename,unsigned int desired_channels)
{
	auto extloc=filename.find_last_of(".");
	if(extloc != std::string::npos)
	{
		if(filename[extloc+1]=='h' || filename[extloc+1]=='H')
		{
			return load_image<float>(filename);
		}
	}
	return load_image<uint8_t>(filename);
}


template<class T>
void save_image(const Array<T,3,RowMajorOrder>& img,const std::string& filename);
template<> void save_image(const Array<uint8_t,3,RowMajorOrder>& img,const std::string& filename);
template<> void save_image(const Array<float,3,RowMajorOrder>& img,const std::string& filename);
//template<> void save_image(const Array<uint16_t,3,RowMajorOrder>& img,const std::string& filename);

template<class T>
void save_image(const Array<T,3,RowMajorOrder>& img,const std::string& filename)
{
	auto extloc=filename.find_last_of(".");
	if(extloc != std::string::npos)
	{
		if(filename[extloc+1]=='h' || filename[extloc+1]=='H')
		{
			save_image<float>(img,filename);
			return;
		}
	}
	save_image<uint8_t>(img,filename);
	return;
}

}

}

#endif
