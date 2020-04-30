#ifndef ND_OPERATIONS_CONVOLVE_HPP
#define ND_OPERATIONS_CONVOLVE_HPP

#include <nd/Array.hpp>


namespace impl
{
	template<class V,class Kernel1D,class BorderHandler>
	void conv1d(V* vob,const V* vb,size_t vn,const Kernel1D& k,size_t klen)
	{
		ptrdiff_t k2=klen>>1;
		ptrdiff_t bottom=-k2+(klen & 1);
		ptrdiff_t top=k2;
		
		ptrdiff_t koffs=klen-1+bottom
		for(ptrdiff_t vi=0;vi<vn;vi++)
		{
			V vo={};
			for(ptrdiff_t ki=0;ki<klen;i++)
			{
				//ptrdiff_t off=klen-ki-1+bottom;
				ptrdiff_t off = koffs-ki;
				vo+=k[ki]*BorderHandler::sample(vb,vn,vi+offs);
			}
			vob[vi]=vo;
		}
	}
}



template<class V,unsigned int D,class StG,class Kernel1D,class BorderHandler>
void Convolve1dInPlace(Array<V,D,StG>& a,unsigned int d,const Kernel1D& kernel,size_t klen) {
	using value_type=V;
	
	size_t N=a.shape(d);
	std::vector<value_type> storagein(N);
	std::vector<value_type> storageout(N);
	std::vector<size_t> indices(N);
	
	for(size_t pid=0;pid<a.size();pid++)
	{
		shape_type pcoord=a.unravel(pid);
		if(pcoord[d]!=0)
		{
			continue;
		}
		
		for(size_t i=0;i<N;i++)
		{
			size_t dex=pcoord.ravel_neighbor_offset(pid,d,i);
			storagein[i]=a[dex];			
			indices[i]=dex;
		}
		conv1d<(storageout.data(),storagein.data(),storagein.size(),kernel,klen);
		for(size_t i=0;i<N;i++)
		{
			a[indices[i]]=storageout[i];		
		}
	}
}

template<class ArrayA,class ArrayB,class BorderHandler>
auto ConvolveSame(const ArrayA& a,const ArrayB& k)->
	typename std::enable_if<ArrayA::num_dimensions==ArrayB::num_dimensions,
		Array<decltype(a[0]*k[0]+a[0]*k[0]),ArrayA::num_dimensions,typename ArrayA::storage_type>
	>::type
{
	using value_type=decltype(a[0]*k[0]+a[0]*k[0]);
	using ArrayOut=Array<value_type,ArrayA::num_dimensions,typename ArrayA::storage_type>;
	
	using itype=typename ArrayOut::index_type;
	
	ArrayOut outarray(a.shape(),value_type{});
	

	size_t aN=a.size();
	size_t kN=k.size();
	for(size_t aid=0;aid<aN;aid++)
	{
		itype acoords=a.unravel(aid);
		for(size_t kid=0;kid<kN;kid++)
		{
			itype bcoords=b.unravel(bid);
			
		}
	}
	
}

}


#endif
