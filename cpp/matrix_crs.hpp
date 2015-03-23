// matrix_crs.hpp
//
// James Folberth
// Spring 2015

#ifndef _MATRIX_CRS_HPP_
#define _MATRIX_CRS_HPP_

#include <algorithm>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <valarray>
#include <vector>

#include "utils.hpp"
#include "matrix_coo.hpp"

using namespace std;

// Forward declare to make function declaration possible
template<typename T>
class matrix_crs;

template<typename T>
class matrix_coo;

template<typename T>
ostream& operator<<(ostream& os, const matrix_crs<T>& mat); 

template <typename T>
class matrix_crs {
   public:
      using value_type = T;

      vector<unsigned> row_ptr;
      vector<unsigned> col_ind;
      vector<T> val;
      size_t m;          // number of rows
      size_t n;          // number of cols

      //////////////////////////////////
      // Construction and Destruction //
      //////////////////////////////////
      matrix_crs<T>() = default;

      matrix_crs<T>(vector<unsigned>& init_row_ind,
                 vector<unsigned>& init_col_ind,
                 vector<T>& init_val,
                 size_t init_m=0, size_t init_n=0, unsigned flag=0);

      // copy, move, destruct
      matrix_crs<T>(const matrix_crs<T>& ) = default;
      matrix_crs<T>& operator=(const matrix_crs<T>& ) = default;
      matrix_crs<T>(matrix_crs<T>&& ) = default;
      matrix_crs<T>& operator=(matrix_crs<T>&& ) = default;
      ~matrix_crs<T>() = default;

      void clean(void);
      void deepclean(void); // calls COO clean; hopefully shouldn't use

      ////////////////
      // Operations //
      ////////////////
      // scalar
      // Didn't implement +=, + as that would obliterate sparsity
      matrix_crs<T>& operator*=(const T& value);
      matrix_crs<T>& operator/=(const T& value);

      // matrix add/sub
      matrix_crs<T>& operator+=(const matrix_crs<T>& B);
      matrix_crs<T>& operator-=(const matrix_crs<T>& B);


      /////////////////////
      // Type conversion //
      /////////////////////
      matrix_coo<T> to_coo(void);
      matrix_crs<T> to_crs(void);

      ////////////
      // Output //
      ////////////
      // Defined below, not in matrix_crs.cpp
      friend ostream& operator<< <>(ostream& os, const matrix_crs& mat); 

      void print_full(void);
 
};

///////////////////////////
// Non-member Operations //
///////////////////////////
// scalar
template<typename T>
matrix_crs<T> operator*(const matrix_crs<T>& lhs, const T& rhs);

template<typename T>
matrix_crs<T> operator*(const T& lhs, const matrix_crs<T>& rhs);

template<typename T>
matrix_crs<T> operator/(const matrix_crs<T>& lhs, const T& rhs);


// matrix add/sub
template<typename T>
matrix_crs<T> operator+(const matrix_crs<T>& lhs, const matrix_crs<T>& rhs);

template<typename T>
matrix_crs<T> operator-(const matrix_crs<T>& lhs, const matrix_crs<T>& rhs);


// matrix-vector product
template<typename T>
valarray<T> operator*(const matrix_crs<T>& A, const valarray<T>& x);



//////////////////////////
// Special Constructors //
//////////////////////////
template<typename T>
matrix_crs<T> eye_crs(unsigned m, unsigned n);

template<typename T>
matrix_crs<T> kron(const matrix_crs<T>& A, const matrix_crs<T>& B);



////////////
// Output //
////////////
template <typename T>
ostream& operator<<(ostream& os, const matrix_crs<T>& mat) {

   os << "Compressed Row Storage (CRS) (rows = " << mat.m << ", cols = " << 
      mat.n << ", nnz = " << mat.val.size() << ")" << endl; // << endl;

   if ( _DEBUG_ >= 1 ) {
      os << "debug 1: CRS ostream printing" << endl;
      os << "debug 1: row_ptr: ";
      for (unsigned i=0; i<mat.row_ptr.size(); ++i)
         os << mat.row_ptr[i] << "  ";
      os << endl;
      os << "debug 1: col_ind: ";
      for (unsigned i=0; i<mat.col_ind.size(); ++i)
         os << mat.col_ind[i] << "  ";
      os << endl;
      os << "debug 1: val: ";
      for (unsigned i=0; i<mat.val.size(); ++i)
         os << mat.val[i] << "  ";
      os << endl;
   }

   if ( mat.val.size() == 0 ) os << "Empty matrix";

   for (unsigned i=0; i<mat.m; ++i) {
      for (unsigned j=mat.row_ptr[i]; j<mat.row_ptr[i+1]; ++j) {
         cout << "  (" << i << ", " << mat.col_ind[j] << ") -> " 
            << setprecision(_PRINT_SPARSE_PREC_) << mat.val[j] << endl;
      }
   }

   return os;
}

#endif
