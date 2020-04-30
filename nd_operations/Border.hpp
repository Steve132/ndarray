
#ifndef ND_OPERATIONS_BORDER_HPP
#define ND_OPERATIONS_BORDER_HPP

#include <nd/Array.hpp>

namespace nd
{
	namespace Border
	{
		struct Wrap
		{
			template<class VFunc>
			static auto sample(const VFunc& v,size_t n,ptrdiff_t k) -> decltype(v(k))
			{
				if(k >= n) k-=n;
				else if(k < 0) k+=n;
				return v(k);
			}
		};
		struct Clamp
		{
			template<class V>
			static auto sample(const VFunc& v,size_t n,ptrdiff_t k) -> decltype(v(k))
			{
				k=std::max(k,0);
				k=std::min(k,n-1);
				return v(k);
			}
		};
				
		struct Zero
		{
			template<class V>
			static auto sample(const VFunc& v,size_t n,ptrdiff_t k) -> decltype(v(k))
			{
				return (k >= n || k < 0) ? {} : v(k);
			}
		};
		
		struct Mirror
		{
			template<class V>
			static auto sample(const VFunc& v,size_t n,ptrdiff_t k) -> decltype(v(k))
			{
				if(k>=n) k=n-k;
				else if(k < 0) k=-k;
				return v(k);
			}
		};
	}
}
