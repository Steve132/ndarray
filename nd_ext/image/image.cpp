#include "image.hpp"
using namespace nd;

#define STB_IMAGE_IMPLEMENTATION
#include "deps/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "deps/stb_image_write.h"

template<> Array<uint8_t,3,ColMajorOrder> nd::image::load_image(const std::string& filename,unsigned int desired_channels)
{
	int x,y,c;
	uint8_t* dat=stbi_load(filename.c_str(),&x,&y,&c,desired_channels);
	if(desired_channels!=0) c=desired_channels;
	Array<uint8_t,3,ColMajorOrder> img({(size_t)c,(size_t)x,(size_t)y});
	size_t N=x;
	N*=y;N*=c;
	memcpy(img.data(),dat,N*sizeof(uint8_t));
	stbi_image_free(dat);
	return img;
}
template<> Array<float,3,ColMajorOrder> nd::image::load_image(const std::string& filename,unsigned int desired_channels)
{
	int x,y,c;
	float* dat=stbi_loadf(filename.c_str(),&x,&y,&c,desired_channels);
	if(desired_channels!=0) c=desired_channels;
	Array<float,3,ColMajorOrder> img({(size_t)c,(size_t)x,(size_t)y});
	size_t N=x;
	N*=y;N*=c;
	memcpy(img.data(),dat,N*sizeof(float));
	stbi_image_free(dat);
	return img;
}
template<> Array<uint16_t,3,ColMajorOrder> nd::image::load_image(const std::string& filename,unsigned int desired_channels)
{
	int x,y,c;
	uint16_t* dat=stbi_load_16(filename.c_str(),&x,&y,&c,desired_channels);
	if(desired_channels!=0) c=desired_channels;
	Array<uint16_t,3,ColMajorOrder> img({(size_t)c,(size_t)x,(size_t)y});
	size_t N=x;
	N*=y;N*=c;
	memcpy(img.data(),dat,N*sizeof(uint16_t));
	stbi_image_free(dat);
	return img;
}
static char get_ext(const std::string& filename)
{
	auto extloc=filename.find_last_of(".");
	if(extloc != std::string::npos)
	{
		return filename[extloc+1];
	}
	return ' ';
}
template<> void nd::image::save_image(const Array<uint8_t,3,ColMajorOrder>& img,const std::string& filename)
{
	char ext=get_ext(filename);
	if(ext=='h' && ext=='H')
	{
		Array<float,3,ColMajorOrder> cimg;
		cimg=img;
		save_image(cimg,filename);
		return;
	}
	int result;
	switch(ext)
	{
		case 'p':
		case 'P':
			result=stbi_write_png(filename.c_str(),img.shape(0),img.shape(1),img.shape(2),img.data(),img.shape(0)*img.shape(1)*sizeof(uint8_t));
			break;
		case 't':
		case 'T':
			result=stbi_write_tga(filename.c_str(),img.shape(0),img.shape(1),img.shape(2),img.data());
			break;
		case 'j':
		case 'J':
			result=stbi_write_jpg(filename.c_str(),img.shape(0),img.shape(1),img.shape(2),img.data(),0);
			break;
		case 'b':
		case 'B':
			result=stbi_write_bmp(filename.c_str(),img.shape(0),img.shape(1),img.shape(2),img.data());
			break;
		default:
			throw std::runtime_error(std::string("Unrecognized Extension in file:")+filename);
	};
	if(result != 0)
	{
		throw std::runtime_error(std::string("Failed to save HDR file:")+filename);
	}
}
template<> void nd::image::save_image(const Array<float,3,ColMajorOrder>& img,const std::string& filename)
{
	char ext=get_ext(filename);
	if(ext!='h' && ext!='H')
	{
		Array<uint8_t,3,ColMajorOrder> cimg;
		cimg=img;
		save_image(cimg,filename);
		return;
	}
	int result=stbi_write_hdr(filename.c_str(),img.shape(0),img.shape(1),img.shape(2),img.data());
	if(result != 0)
	{
		throw std::runtime_error(std::string("Failed to save HDR file:")+filename);
	}
}
/*template<> void nd::image::save_image(const Array<uint16_t,3,ColMajorOrder>& img,const std::string& filename)
{
	
}*/
