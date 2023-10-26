#ifndef BCKG_COMMON_H
#define BCKG_COMMON_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <vector>
#include <typeinfo>
#include <limits>
#include <stdlib.h>
#include <stdio.h>
#include <stdexcept>
#include <string.h>
#include <math.h>

#if   defined(_WIN32)
    #include <direct.h>
    #define GetCurrentDir _getcwd
#elif defined(_WIN64)
    #include <direct.h>
    #define GetCurrentDir _getcwd
#else
    #include <unistd.h>
    #define GetCurrentDir getcwd
#endif

#define FILENAME_UMAX 1024


namespace com
{
/*******************************************************************************************/

inline std::string get_current_dir(void)
{
    char buff[FILENAME_UMAX];
    char* tmp = GetCurrentDir( buff, FILENAME_UMAX );

    std::string current_working_dir(buff);
    return current_working_dir;
}

/*******************************************************************************************/

inline char slash(void)
{
#if   defined(_WIN32)
    return '\\';
#elif defined(_WIN64)
    return '\\';
#else
    return '/';
#endif
}

/*******************************************************************************************/

inline char ascii_toupper_char(char c)
{
    return ('a' <= c && c <= 'z') ? c^0x20 : c;
}

/*******************************************************************************************/

size_t strtoupper_autovec(char* dst, size_t lendst, const char* src)
{
    size_t lensrc = strlen(src);
    size_t len = (lensrc <= lendst) ? lensrc : lendst;      //std::cout<<"len src = "<<lensrc<<", len dst = "<<lendst<<", len  = "<<len<<"\n\n";

    if ( len < lensrc )
    {
        int out = printf("\nError: insufficient length of the destiny string.\a\a\a\n\n");
        exit(0);
    }

    for (size_t i=0; i<len; ++i)   dst[i] = ascii_toupper_char(src[i]);

    if ( len < lendst )   {   for (size_t i=len; i<lendst; ++i)  dst[i] = '\0';   }
    return len;
}

/*******************************************************************************************/
/*******************************************************************************************/

template <typename T> inline T Pi(void) //!<pi number
{
    return static_cast<T>(3.1415926535897932384626433832795028841971693993751);
}

/*******************************************************************************************/

template <typename T> inline T Eps(void) //!<numerical "zero"
{
    return std::numeric_limits<T>::epsilon();
}

/*******************************************************************************************/

template <typename T> inline T Inf(void) //!<numerical "+infinity"
{
    if(std::numeric_limits<T>::has_infinity) return std::numeric_limits<T>::infinity();
    else                                     return std::numeric_limits<T>::max();
}

/*******************************************************************************************/

template <typename T> inline T NaN(void) //!<numerical "NaN"
{
    return std::numeric_limits<T>::quiet_NaN();
}

/*******************************************************************************************/
/*******************************************************************************************/

template <typename T> inline T Mp(void) //!<proton mass, SI
{
    return static_cast<T>(1.6726219e-27);
}

/*******************************************************************************************/

template <typename T> inline T Me(void) //!<electron mass, SI
{
    return static_cast<T>(9.10938356e-31);
}

/*******************************************************************************************/

template <typename T> inline T ElC(void) //!<elementary charge, SI
{
    return static_cast<T>(1.6021766208989898989898989898989898989898989898e-19);
}

/*******************************************************************************************/

template <typename T> inline T Kb(void) //!<Boltzmann constant, SI
{
    return static_cast<T>(1.3806485279797979797979797979797979797979797979e-23);
}

/*******************************************************************************************/

template <typename T> inline T LS(void) //!<light speed in vacuum, SI
{
    return static_cast<T>( 299792458.0 );
}

/*******************************************************************************************/

template <typename T> inline T MUo(void) //!<vacuum permeability, SI
{
    return static_cast<T>( 4.0e-7*Pi<T>() );
}

/*******************************************************************************************/

template <typename T> inline T Eo(void) //!<vacuum permittivity, SI
{
    return static_cast<T>( 1.0 / MUo<T>() / LS<T>() / LS<T>() );
}

/*******************************************************************************************/
/*******************************************************************************************/

template <typename T> inline T abs(const T & x)
{
    return (x > 0) ? x : -x;
}

/*******************************************************************************************/

/// round numerical value to some number of digits behind the comma
template <typename T> inline T round_decdig(T x, size_t dig=15)
{
    return  static_cast<T>( round(x*pow(10.0,dig)) / pow(10.0,dig) );
}

/*******************************************************************************************/

/// round numerical value to some number of valuable digits
template <typename T> T round_valdig(T x, int dig=3)
{
    int j = 0; int dj = 0;
    double y;

    if      ( abs<T>(x) == 0.0 )  return  x;
    else if ( abs<T>(x) >  1.0 )  dj = -1;
    else if ( abs<T>(x) <  1.0 )  dj = +1;
    else                          j += (dig-1);

    if ( dj > 0 )
    {
        do { j += dj;   y = static_cast<double>( floor( abs<T>(x)*pow(10.0,j) ) );
        }   while( abs<double>(y) < 0.5 );
        j += dig-1;
    }

    else if ( dj < 0 )
    {
        do { j += dj;   y = static_cast<double>( floor( abs<T>(x)*pow(10.0,j) ) );
        }   while( abs<double>(y) > 0.5 );
        j += dig;
    }

    return  static_cast<T>( round(x*pow(10.0,j)) / pow(10.0,j) );
}

/*******************************************************************************************/
/*******************************************************************************************/

/// create 3D array
template <typename T> T*** array3_generator(size_t dim0, size_t dim1, size_t dim2)
{
    T*** arr = new T**[dim0];
    for(size_t i0 = 0; i0 < dim0; i0++) {
        arr[i0] = new T*[dim1];
        for(size_t i1 = 0; i1 < dim1; i1++) {
            arr[i0][i1] = new T[dim2];
            for(size_t i2 = 0; i2 < dim2; i2++) {
                arr[i0][i1][i2] = 0;  //std::cout << std::setw(4) << arr[i0][i1][i2];
            }
        }
    }
    return arr;
}

/*******************************************************************************************/

/// destroy 3D array
template <typename T> void array3_destroyer(T*** arr, size_t dim0, size_t dim1)
{
    for (size_t i0 = 0; i0 < dim0; i0++) {
        for (size_t i1 = 0; i1 < dim1; i1++) {
            delete [] arr[i0][i1];
        }
        delete [] arr[i0];
    }
    delete [] arr;
}

/*******************************************************************************************/

/// write 3D array in a binary grid-file
template <typename T> void array3_grid_file_write(size_t dim0, size_t dim1, size_t dim2,
                          T* X1, T* X2, T*** VAL, std::string file_name)
{
    int tmp = std::remove(file_name.c_str());     /// delete old file if any

    std::ofstream file;
    file.open(file_name.c_str(), std::ofstream::out | std::ofstream::binary);

    tmp = (int)dim0; file.write((char *)&tmp, sizeof(int));         /// rank of state vector

    tmp = 2;         file.write((char *)&tmp, sizeof(int));         /// space dimension

    tmp = (int)dim1; file.write((char *)&tmp, sizeof(int));         /// length x
    tmp = (int)dim2; file.write((char *)&tmp, sizeof(int));         /// length z

    file.write((char *)X1,  dim1 * sizeof(T));                      /// x-coordinates
    file.write((char *)X2,  dim2 * sizeof(T));                      /// z-coordinates

    for (size_t i2=0; i2<dim2; i2++){
        for (size_t i1=0; i1<dim1; i1++){
            for (size_t i0=0; i0<dim0; i0++){
                file.write((char *)&(VAL[i0][i1][i2]), sizeof(T));  /// data array
            }
        }
    }
    file.close();
    std::cout << "\n3d array is written to:\t\t" << file_name << "\n";
    std::cout << "size of a number:\t\t" << sizeof(T) << " B\n\n";
}

/*******************************************************************************************/
/*******************************************************************************************/

/// make uniform 1D vector
template <typename T> void make_uniform_vector(size_t dimu, T umin, T umax, T du, T* u)
{
    if ( dimu < 1 ) {
        int out = printf("Error: inconsistent vector length: %d\a\a\a\n", (int)dimu);
        exit(0);
    }

    u[0] = umin;   if ( dimu == 1 )   return;

    for (size_t j=1; j<dimu; j++) { u[j] = u[j-1]+du; }

    u[dimu-1] = umax;         /// make sure that right boundary is exact
}

/*******************************************************************************************/
/*******************************************************************************************/

} /// end of namespace

#endif /// BCKG_COMMON_H
