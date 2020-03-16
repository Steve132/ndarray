# ndarray
ndArray C++:  A new header-only C++11 N-dimensional Array


## Why should you use this

Think np.array(),(or FORTRAN, or MATLAB or whatever you want), except in C++.

## But there are dozens of image classes?

None that work in 3D or 5D as well as 2D or 1D, and can use templated data members.  Or arbitrary layouts.

## Yes there are:

Blitz++ and boost.multi_array are good options.  But Blitz++ is *very large* and provides a lot of features that you dont' really want to use.  Plus you have to *build it* and have the prerequisites installed.

boost.multi_array requires **boost**.  'nuff said.

This is header-only and it always will be.  No hassles.

Example code:

    #include<ndarray.hpp>
    
    int main(int argc,char** argv)
    {
        nd::Array<double,3> something({2,3,3},30.0); //an 2x2x3 array filled with the value 30.0
        typedef pixel std::array<unsigned char,3>;
        nd::Array<pixel,2,nd::ColumnMajorOrder> image({640,480}); //a 640x480 image of RGB pixels 
    }
    
    
