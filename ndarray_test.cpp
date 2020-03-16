#include "ndarray.hpp"
#include<random>





int main()
{
	nd::Array<int,3> a({50,50,50});
	nd::Array<double,3> b({3,3,3});
	//auto G=a+b;
	nd::Array<double,2> data({3,3});
	data(1,0)=3.5;
	
	
	data+=3.0;
	
	nd::Array<double,3> A({100,100,3});
	nd::Array<double,3> B(A.shape());
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(0.0,1.0);
	A.elementWiseInPlace([&distribution,&generator](double& v){ v=distribution(generator); });
	B.elementWiseInPlace([&distribution,&generator](double& v){ v=distribution(generator); });
	
	
	std::vector<nd::Array<double,3>> allBlocks;

	
	nd::Array<double,3> u=(A+B)/2.0;
	A-=u;
	B-=u;
	
	
	std::cout << data << std::endl;
	std::cout << "nx\n" << std::endl;
	std::cout << b << std::endl;
	return 0;
}
