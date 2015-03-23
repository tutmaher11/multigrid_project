// matrix_crs.cpp
//
// James Folberth
// Spring 2015

#include "matrix_crs.hpp"

using namespace std;

//////////////////////////////////
// Construction and Destruction //
//////////////////////////////////
template<typename T>
matrix_crs<T>::matrix_crs(
      vector<unsigned>& init_row_ind,
      vector<unsigned>& init_col_ind,
      vector<T>& init_val,
      size_t init_m, size_t init_n, unsigned flag) {

   // if a (complete) CRS matrix is passed in, build it
   if ( flag == 1 ) {
      row_ptr = init_row_ind; // actually row pointers
      col_ind = init_col_ind;
      val = init_val;
      m = init_m;
      n = init_n;
      return;
   }
   
// Construct a CRS matrix from vectors of row index, column index, and values
// We first build a COO matrix to sort and deal with duplicate entries, and
// also determine the size of the matrix.
// Then we build the vector of row pointers, taking into account zero rows

   // Build COO matrix (which will sort indexes)
   matrix_coo<T> coomat(init_row_ind, init_col_ind, init_val, init_m, init_n); 

   if ( _DEBUG_ >= 2) {
      cout << "debug 2: CRS matrix construction" << endl;
      cout << coomat << endl;
   }

   col_ind = coomat.col_ind;
   val = coomat.val;
   m = coomat.m;
   n = coomat.n;

   // Build vector of row pointers
   row_ptr.resize(m+1);
   row_ptr[m] = val.size();
   unsigned row_ptr_ind = 0;

   // Test for matrix of zeros
   if (val.size() == 0) {
      fill(row_ptr.begin(), row_ptr.end(), 0);
      return;
   }

   // Handle initial rows of zeros
   for (unsigned j=0; j<coomat.row_ind[0]; ++j)
      row_ptr[row_ptr_ind + j] = 0;
   row_ptr_ind = coomat.row_ind[0];

   for (unsigned i=1; i<col_ind.size(); ++i) {

      // nothing to do; keep indexing along row
      if (coomat.row_ind[i-1] == coomat.row_ind[i]) {

         // unless it's the last row
         if (coomat.row_ind[i] == m) {
            row_ptr[row_ptr_ind+1] = i;
         }
         continue;
      }

      // jumping down to next column
      else if (coomat.row_ind[i-1] == coomat.row_ind[i] - 1) {
         row_ptr[row_ptr_ind+1] = i;
         row_ptr_ind += 1;
      }

      // we have a row or rows of zeros
      // Repeat last index of col_ind for each row we skip
      // Also put in the next row_ptr for the next non-zero row
      else if (coomat.row_ind[i-1] < coomat.row_ind[i] - 1 ) {
         for (unsigned j=1; j<coomat.row_ind[i]-row_ptr_ind; ++j) {
            row_ptr[row_ptr_ind + j] = i;
         }
         row_ptr_ind = coomat.row_ind[i];
         row_ptr[row_ptr_ind] = i;
      }

      else {
         cerr << "error: matrix_crs.cpp::matrix_crs: error constructing "
              << "CRS matrix" << endl;
         exit(-1);
      }
   }

   // Handle rows of zeros at the end of the matrix
   for (unsigned r=*(coomat.row_ind.end()-1)+1; r<m; ++r) {
      row_ptr[r] = col_ind.size();
   }
}


// Remove zero entries.  This assumes a proper CRS format was input
template<typename T>
void matrix_crs<T>::clean(void) {
   unsigned ind=0;

   if ( val.size() == 0 ) return; // nothing to do

   for (unsigned row=0; row < m; ++row) {
      while (val.size() > 0 && ind < row_ptr[row+1]) {
         if ( abs(val[ind]) < _ELEMENT_ZERO_TOL_ ) {
            col_ind.erase(col_ind.begin() + ind);
            val.erase(val.begin() + ind);

            // fix row pointers
            for (unsigned r=row+1; r <= m; ++r) {
               row_ptr[r] -= 1;
            }
            continue;
         }
         ++ind;
      }
   }
}

// Hopefully don't have to use this 'cuz everything is out of order
template<typename T>
inline void matrix_crs<T>::deepclean(void) {
   *this = (this->to_coo()).to_crs();
}



////////////////
// Operations //
////////////////
// scalar
template<typename T>
matrix_crs<T>& matrix_crs<T>::operator*=(const T& value) {
   for (auto it=val.begin(); it != val.end(); ++it) {
      *it *= value;
   }
   return *this;
}

template<typename T>
inline matrix_crs<T> operator*(const matrix_crs<T>& lhs, const T& rhs) {
   return matrix_crs<T>(lhs) *= rhs;
}

template<typename T>
inline matrix_crs<T> operator*(const T& lhs, const matrix_crs<T>& rhs) {
   return matrix_crs<T>(rhs) *= lhs;
}


template<typename T>
matrix_crs<T>& matrix_crs<T>::operator/=(const T& value) {
   for (auto it=val.begin(); it != val.end(); ++it) {
      *it /= value;
   }
   return *this;
}

template<typename T>
inline matrix_crs<T> operator/(const matrix_crs<T>& lhs, const T& rhs) {
   return matrix_crs<T>(lhs) /= rhs;
}



// matrix add/sub
// *this += B;
template<typename T>
matrix_crs<T>& matrix_crs<T>::operator+=(const matrix_crs<T>& B) {
   //unsigned this_col_ind;
   unsigned ptrA, stopA, ptrB, stopB, ptrS;
   unsigned colA, colB;

   vector<unsigned> row_ptrS(m+1);
   vector<unsigned> col_indS(this->val.size() + B.val.size());
   vector<T> valS(this->val.size() + B.val.size());
   // Note: this->val.size() + B.val.size() is an upper bound on the number of
   // nonzeros in the result A+B.  We'll resize at the end;

   // check sizes
   if ( m != B.m || n != B.n ) {
      cerr << "error: matrix_crs:+=: matrix sizes do not match - ("
           << m << "," << n << ") vs. ("
           << B.m << "," << B.n << ")" << endl;
      exit(-1);
   }

   colA = 0; colB = 0; // to get -Wall to shut up
   ptrS = 0;

   // similar to Julia's version
   // This doesn't work in place; it returns a new matrix
   // It's still a lot faster than the one I wrote!
   // TODO: this does some dumb copy stuff.  Better to define operator+ instead
   for (unsigned row=0; row < m; ++row) {
      ptrA = row_ptr[row];
      stopA = row_ptr[row+1];
      ptrB = B.row_ptr[row];
      stopB = B.row_ptr[row+1];

      while ( ptrA < stopA && ptrB < stopB ) {
         colA = col_ind[ptrA];
         colB = B.col_ind[ptrB];

         // entry in *this (A) comes first
         if ( colA < colB ) {
            col_indS[ptrS] = colA;
            valS[ptrS] = val[ptrA];
            ptrS += 1;
            ptrA += 1;
         }

         // entry in B comes first
         else if ( colA > colB) {
            col_indS[ptrS] = colB;
            valS[ptrS] = B.val[ptrB];
            ptrS += 1;
            ptrB += 1;
         }

         else {
            col_indS[ptrS] = colA;
            valS[ptrS] = val[ptrA] + B.val[ptrB];

            // if nonzero element, then add it; if not, it'll get overwritten
            // or put past the end of the value vector
            if ( abs(valS[ptrS]) > _ELEMENT_ZERO_TOL_ )
               ptrS += 1;

            ptrA += 1;
            ptrB += 1;
         }
      }

      while ( ptrA < stopA ) {
         col_indS[ptrS] = col_ind[ptrA];
         valS[ptrS] = val[ptrA];
         ptrS += 1;
         ptrA += 1;
      }

      while ( ptrB < stopB ) {
         col_indS[ptrS] = B.col_ind[ptrB];
         valS[ptrS] = B.val[ptrB];
         ptrS += 1;
         ptrB += 1;
      }
      
      row_ptrS[row+1] = ptrS;
   }

   col_indS.resize(ptrS);
   valS.resize(ptrS);

   row_ptr = row_ptrS;
   col_ind = col_indS;
   val = valS;

   return *this;
}

template<typename T>
matrix_crs<T> operator+(const matrix_crs<T>& lhs, const matrix_crs<T>& rhs) {
   return matrix_crs<T>(lhs) += rhs;
}

template<typename T>
matrix_crs<T>& matrix_crs<T>::operator-=(const matrix_crs<T>& B) {
   return (*this += -1.0*B);
}

template<typename T>
matrix_crs<T> operator-(const matrix_crs<T>& lhs, const matrix_crs<T>& rhs) {
   return matrix_crs<T>(lhs) -= rhs;
}


// matrix-vector product
template<typename T>
valarray<T> operator*(const matrix_crs<T>& A, const valarray<T>& x) {
  
   // check sizes
   if ( A.n != x.size() ) {
      cerr << "error: matrix_crs:*: matrix-vector product dimension mismatch."
           << endl;
      exit(-1);
   }

   // allocate and fill with zeros
   valarray<T> b(static_cast<T>(0.),A.m);

   for (size_t row=0; row < A.m; ++row) {
      for (size_t ptr=A.row_ptr[row]; ptr<A.row_ptr[row+1]; ++ptr) {
         b[row] += A.val[ptr] * x[A.col_ind[ptr]];
      }
   }

   return b;
}




/////////////////////
// Type conversion //
/////////////////////
template<typename T>
matrix_coo<T> matrix_crs<T>::to_coo(void) {
   vector<unsigned> new_row_ind(col_ind.size());

   for (unsigned i=0; i<row_ptr.size()-1; ++i) {
      fill(new_row_ind.begin()+row_ptr[i],
           new_row_ind.begin()+row_ptr[i+1],
           i);
   }

   return matrix_coo<T>(new_row_ind, col_ind, val, m, n);
}

template<typename T>
inline matrix_crs<T> matrix_crs<T>::to_crs(void) {
   return *this;
}

//////////////////////////
// Special Constructors //
//////////////////////////
template<typename T>
matrix_crs<T> eye_crs(unsigned m, unsigned n) {
   vector<unsigned> rind(MIN(m,n)), cind(MIN(m,n));
   vector<T> val(MIN(m,n), static_cast<T>(1)); // fill val with 1s on construct

   for (unsigned i=0; i < rind.size(); ++i) {
      rind[i] = i;
      cind[i] = i;
   }

   return matrix_crs<T>(rind, cind, val, m, n);
}

template<typename T>
matrix_crs<T> kron(const matrix_crs<T>& A, const matrix_crs<T>& B) {
   unsigned mA,nA,mB,nB,mK,nK,nnzK;
   unsigned startA, stopA, ptrA, lA;
   unsigned startB, stopB, ptrB, lB;
   unsigned ptrK, row, ptrK_ran;

   mA = A.m; nA = A.n; mB = B.m; nB = B.n;
   mK = mA*mB; nK = nA*nB;
   nnzK = A.val.size()*B.val.size();

   vector<unsigned> row_ptrK(mK+1);
   vector<unsigned> col_indK(nnzK);
   vector<T> valK(nnzK);// this is the exact size

   // this is similar to Julia's
   row = 0;
   row_ptrK[0] = 0;

   for (unsigned j=0; j<mA; ++j) {
      startA = A.row_ptr[j];
      stopA = A.row_ptr[j+1];
      lA = stopA-startA;

      for (unsigned i=0; i<mB; ++i) {
         startB = B.row_ptr[i];
         stopB = B.row_ptr[i+1];
         lB = stopB-startB;

         ptrK_ran = row_ptrK[row];

         for (ptrA = startA; ptrA < stopA; ++ptrA) {
            ptrB = startB;

            for (ptrK = ptrK_ran; ptrK < ptrK_ran+lB; ++ptrK) {
               col_indK[ptrK] = A.col_ind[ptrA]*nB + B.col_ind[ptrB];
               valK[ptrK] = A.val[ptrA] * B.val[ptrB];
               ptrB += 1;
            }
            
            ptrK_ran += lB;
         }

         row_ptrK[row+1] = row_ptrK[row] + lA*lB;
         row += 1;
      }
   }

   return matrix_crs<T>(row_ptrK, col_indK, valK, mK, nK, 1);
}



////////////
// Output //
////////////
template<typename T>
void matrix_crs<T>::print_full(void) {
// Print sparse matrix as a dense array

   unsigned ind=0;
   for (unsigned r=0; r<m; ++r) {
      for (unsigned c=0; c<n; ++c) {

         // we've hit a nonzero in this row
         if (col_ind.size() > 0 &&
             ind < row_ptr[r+1] &&
             ind < col_ind.size() && // to appease valgrind memcheck
             c == col_ind[ind] ) { 
            cout << ' ' << setw(_PRINT_FULL_WIDTH_) << setfill(' ')
                 << setprecision(_PRINT_FULL_PREC_)
                 << static_cast<double>(val[ind]);
            ind += 1;
         }

         // print a zero
         else {
            cout << ' ' << setw(_PRINT_FULL_WIDTH_) << setfill(' ')
                 << setprecision(_PRINT_FULL_PREC_)
                 << 0e0;
         }

         // end of row, so print endl
         if (c == n-1) cout << endl;
      }
   }

   cout << endl;
}

// TODO is there a better way?
////////////////////////////////////////////
// Force instantiation for specific types //
////////////////////////////////////////////
// double
/////////
template class matrix_crs<double>;
template matrix_crs<double> eye_crs<double>(unsigned,unsigned);
template matrix_crs<double> kron<double>(const matrix_crs<double>&,
                                         const matrix_crs<double>&);

template matrix_crs<double> operator*(const matrix_crs<double>&,const double&);
template matrix_crs<double> operator*(const double&,const matrix_crs<double>&);
template matrix_crs<double> operator/(const matrix_crs<double>&,const double&);

template matrix_crs<double> operator+(const matrix_crs<double>&,
                                      const matrix_crs<double>&);
template matrix_crs<double> operator-(const matrix_crs<double>&,
                                      const matrix_crs<double>&);

template valarray<double> operator*(const matrix_crs<double>&,
                                      const valarray<double>&);

