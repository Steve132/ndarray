

#ifndef NDARRAY_SRC_OPERATORS_HPP
#define NDARRAY_SRC_OPERATORS_HPP

#include<type_traits>
#include<iostream>


namespace nd
{

template<class ArrayClass,typename std::enable_if<ArrayClass::num_dimensions==1,int>::type = 0>
std::ostream& operator<<(std::ostream& out,const ArrayClass& arr)
{
	for(size_t i=0;i<arr.shape(0);i++) out << arr(i) << " ";
	return out;
}
template<class ArrayClass,typename std::enable_if<ArrayClass::num_dimensions==2,int>::type = 0>
std::ostream& operator<<(std::ostream& out,const ArrayClass& arr)
{
	for(size_t i=0;i<arr.shape(0);i++) 
	{
		for(size_t j=0;j<arr.shape(1);j++) out << arr(i,j) << " ";
		out << "\n";
	}
	return out;
}
template<class ArrayClass,typename std::enable_if<(ArrayClass::num_dimensions > 2),int>::type = 0>
std::ostream& operator<<(std::ostream& out,const ArrayClass& arr)
{
	for(size_t i=0;i<arr.shape(ArrayClass::num_dimensions-1);i++)
	{
		auto outmap=arr.pop_order_ref(i);
		out << "Dim[" << ArrayClass::num_dimensions-1 << "]=" << i << "\n";
		out << outmap;
	}
	return out;
}

//external values because SFNIAE for operator selection.
#define BINARY_OPERATOR_TEMPLATE( XXX ) \
template<class V,unsigned int D,class StG,class V2,class StG2> \
auto operator XXX(const Array<V,D,StG>& a,const Array<V2,D,StG2>& b) \
	->Array<decltype(a[0] XXX b[0]),D,StG> \
{ \
	return a.elementWise([](const V& v1,const V2& v2){ return v1 XXX v2; },b); \
} \
template<class V,unsigned int D,class StG,class TypeB> \
auto operator XXX(const Array<V,D,StG>& a,const TypeB& b) \
	-> Array<decltype(a[0] XXX b),D,StG> \
{ \
	return a.elementWise([&b](const V& v1){ return v1 XXX b; }); \
} \
template<class V,unsigned int D,class StG,class TypeB> \
auto operator XXX(const TypeB& b,const Array<V,D,StG>& a) \
-> Array<decltype(b XXX a[0]),D,StG> \
{ \
	return a.elementWise([&b](const V& v1){ return b XXX v1; }); \
} \

BINARY_OPERATOR_TEMPLATE(*)
BINARY_OPERATOR_TEMPLATE(/)
BINARY_OPERATOR_TEMPLATE(%)
BINARY_OPERATOR_TEMPLATE(+)
BINARY_OPERATOR_TEMPLATE(-)
BINARY_OPERATOR_TEMPLATE(<<)
BINARY_OPERATOR_TEMPLATE(>>)
BINARY_OPERATOR_TEMPLATE(<)
BINARY_OPERATOR_TEMPLATE(>)
BINARY_OPERATOR_TEMPLATE(<=)
BINARY_OPERATOR_TEMPLATE(>=)
BINARY_OPERATOR_TEMPLATE(==)
BINARY_OPERATOR_TEMPLATE(!=)
BINARY_OPERATOR_TEMPLATE(&)
BINARY_OPERATOR_TEMPLATE(|)
BINARY_OPERATOR_TEMPLATE(^)
BINARY_OPERATOR_TEMPLATE(&&)
BINARY_OPERATOR_TEMPLATE(||)

#define UNARY_OPERATOR_TEMPLATE( XXX ) \
template<class V,unsigned int D, class StG> \
auto operator XXX(const Array<V,D,StG>& a) \
->Array<decltype(XXX(a[0])),D,StG> \
{ \
	return a.elementWise([](const V& v1){ return XXX(v1); }); \
} \

UNARY_OPERATOR_TEMPLATE(~)
UNARY_OPERATOR_TEMPLATE(!)
UNARY_OPERATOR_TEMPLATE(-)
UNARY_OPERATOR_TEMPLATE(+)

#define INPLACE_OPERATOR_TEMPLATE( XXX ) \
template<class V,unsigned int D,class StG,class V2, class StG2> \
auto operator XXX(Array<V,D,StG>& a,const Array<V2,D,StG2>& b) \
->Array<V,D,StG>& \
{ \
	return a.elementWiseInPlace([](V& v1,const V2& v2){ v1 XXX v2; },b); \
} \
template<class V,unsigned int D,class StG,class TypeB> \
auto operator XXX(Array<V,D,StG>& a,const TypeB& b) \
-> Array<V,D,StG>& \
{ \
	return a.elementWiseInPlace([&b](V& v1){ v1 XXX b; }); \
} \

INPLACE_OPERATOR_TEMPLATE(*=)
INPLACE_OPERATOR_TEMPLATE(/=)
INPLACE_OPERATOR_TEMPLATE(+=)
INPLACE_OPERATOR_TEMPLATE(-=)
INPLACE_OPERATOR_TEMPLATE(>>=)
INPLACE_OPERATOR_TEMPLATE(<<=)
INPLACE_OPERATOR_TEMPLATE(&=)
INPLACE_OPERATOR_TEMPLATE(^=)
INPLACE_OPERATOR_TEMPLATE(|=)
//TODO: ++,--

}

#endif
