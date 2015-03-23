// matrix_coo.hpp
//
// James Folberth
// Spring 2015

#ifndef _MATRIX_COO_HPP_
#define _MATRIX_COO_HPP_

#include <algorithm>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <vector>

#include "utils.hpp"
#include "matrix_crs.hpp"

using namespace std;

// Forward declare to make function declaration possible
template<typename T>
class matrix_coo;

template<typename T>
class matrix_crs;

template<typename T>
ostream& operator<<(ostream& os, const matrix_coo<T>& mat); 

template <typename T>
class matrix_coo {
   public:
      using value_type = T;

      vector<unsigned> row_ind;
      vector<unsigned> col_ind;
      vector<T> val;
      size_t m;          // number of rows
      size_t n;          // number of cols

      //////////////////////////////////
      // Construction and Destruction //
      //////////////////////////////////
      matrix_coo<T>() = default;

      matrix_coo<T>(vector<unsigned>& init_row_ind,
                 vector<unsigned>& init_col_ind,
                 vector<T>& init_val,
                 size_t init_m=0, size_t init_n=0);

      // copy, move, destruct
      matrix_coo<T>(const matrix_coo<T>& ) = default;
      matrix_coo<T>& operator=(const matrix_coo<T>& ) = default;
      matrix_coo<T>(matrix_coo<T>&& ) = default;
      matrix_coo<T>& operator=(matrix_coo<T>&& ) = default;
      ~matrix_coo<T>() = default;

      // clean up storage (sort inds, deal with duplicates,
      // remove elements below _ELEMENT_ZERO_TOL_)
      void clean(void);


      ////////////////
      // Operations //
      ////////////////
      // scalar
      matrix_coo<T>& operator*=(const T& value);
      matrix_coo<T>& operator/=(const T& value);

      // matrix add/sub
      matrix_coo<T>& operator+=(const matrix_coo<T>& B);
      matrix_coo<T>& operator-=(const matrix_coo<T>& B);


      /////////////////////
      // Type conversion //
      /////////////////////
      matrix_crs<T> to_crs(void);


      ////////////
      // Output //
      ////////////
      // Defined below, not in matrix_coo.cpp
      friend ostream& operator<< <>(ostream& os, const matrix_coo& mat); 

      void print_full(void);
 
};


///////////////////////////
// Non-member Operations //
///////////////////////////
// scalar
template<typename T>
matrix_coo<T> operator*(const matrix_coo<T>& lhs, const T& rhs);

template<typename T>
matrix_coo<T> operator*(const T& lhs, const matrix_coo<T>& rhs);

template<typename T>
matrix_coo<T> operator/(const matrix_coo<T>& lhs, const T& rhs);


// matrix add/sub
template<typename T>
matrix_coo<T> operator+(const matrix_coo<T>& lhs, const matrix_coo<T>& rhs);

template<typename T>
matrix_coo<T> operator-(const matrix_coo<T>& lhs, const matrix_coo<T>& rhs);



//////////////////////////
// Special Constructors //
//////////////////////////
template<typename T>
matrix_coo<T> eye_coo(unsigned m, unsigned n);


////////////
// Output //
////////////
template <typename T>
ostream& operator<<(ostream& os, const matrix_coo<T>& mat) {

   os << "Coordinate Sparse (COO) (rows = " << mat.m << ", cols = " << 
      mat.n << ", nnz = " << mat.val.size() << ")" << endl; // << endl;

   if ( _DEBUG_ >=1 ) {
      os << "debug 1: COO ostream printing" << endl;
      os << "debug 1: row_ind: ";
      for (unsigned i=0; i < mat.row_ind.size(); ++i)
         cout << mat.row_ind[i] << "  ";
      os << endl;
      os << "debug 1: col_ind: ";
      for (unsigned i=0; i < mat.col_ind.size(); ++i)
         cout << mat.col_ind[i] << "  ";
      os << endl;
      os << "debug 1: val: ";
      for (unsigned i=0; i < mat.val.size(); ++i)
         cout << setprecision(_PRINT_SPARSE_PREC_) << mat.val[i] << "  ";
      os << endl;
   }

   if ( mat.val.size() == 0 ) os << "Empty matrix";

   for ( unsigned i=0; i < mat.row_ind.size(); ++i) {
      os << "  (" << mat.row_ind[i] << ", " << mat.col_ind[i]
         << ") -> " << setprecision(_PRINT_SPARSE_PREC_)
         << mat.val[i] << endl;
   }

   return os;
}

#endif
