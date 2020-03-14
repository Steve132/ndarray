#include "ndarray.hpp"

int main()
{
	nd::LinearArray<double,3> a({3,3,3});
	nd::LinearArray<double,3> b({3,3,3});
	auto G=a+b;
	nd::Array<double,2> data({3,3});
	data(1,0)=3.5;
	std::cout << data << std::endl;
	return 0;
}
