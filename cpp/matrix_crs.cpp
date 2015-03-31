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
      size_t init_m, size_t init_n, int crsflag, int sortflag) {

   // if a (complete) CRS matrix is passed in, build it
   if ( crsflag == 1 ) {
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
   matrix_coo<T> coomat(init_row_ind, init_col_ind, init_val, init_m, init_n,
         sortflag); 

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
// TODO: it might be faster to just create new row_ptr, col_ind, val vectors
template<typename T>
void matrix_crs<T>::clean(void) {

   if ( val.size() == 0 ) return; // nothing to do

   //unsigned ind = 0;
   //for (unsigned row=0; row < m; ++row) {
   //   while (val.size() > 0 && ind < row_ptr[row+1]) {
   //      if ( abs(val[ind]) < _ELEMENT_ZERO_TOL_ ) {
   //         col_ind.erase(col_ind.begin() + ind);
   //         val.erase(val.begin() + ind);

   //         // fix row pointers (TODO this should be accomplished via a num_removed counter)
   //         for (unsigned r=row+1; r <= m; ++r) {
   //            row_ptr[r] -= 1;
   //         }
   //         continue;
   //      }
   //      ++ind;
   //   }
   //}

   // this is better, but could still use some work
   //unsigned ind = 0;
   //for (unsigned row=0; row < m; ++row) {
   //   unsigned num_removed = 0;
   //   unsigned ind_max = row_ptr[row+1];
   //   while ( val.size() > 0 && ind < ind_max ) {
   //      if ( abs(val[ind]) < _ELEMENT_ZERO_TOL_ ) {
   //         col_ind.erase(col_ind.begin() + ind);
   //         val.erase(val.begin() + ind);
   //         ++num_removed;
   //         --ind_max;
   //         continue;
   //      }
   //      else ++ind;
   //   }

   //   // Fix row pointers
   //   for (unsigned r = row+1; r <= m; ++r) {
   //      row_ptr[r] -= num_removed;
   //   }
   //}

   // This copies the non-zero elements of *this to new vectors
   vector<unsigned> rp(row_ptr.size(), 0), ci;
   vector<T> v;
   T biggest = MAX(-*min_element(val.begin(), val.end()), *max_element(val.begin(), val.end()));
   T zero_cutoff = _ELEMENT_ZERO_TOL_ * biggest;

   rp[0] = 0;
   for (unsigned i = 0; i < m; ++i) {
      rp[i+1] = rp[i];

      for (unsigned jp = row_ptr[i]; jp < row_ptr[i+1]; ++jp) {
         if ( abs(val[jp]) > zero_cutoff ) {
            ++rp[i+1];
            ci.push_back(col_ind[jp]);
            v.push_back(val[jp]);
         }
      }
   }

   row_ptr = rp;
   col_ind = ci;
   val = v;
}

// Hopefully don't have to use this 'cuz everything is out of order
template<typename T>
inline void matrix_crs<T>::deepclean(void) {
   *this = (this->to_coo()).to_crs();
}

// Sort the column indexes for each row
// Used after Gustavson's algorithm for sparse mat-mat
template<typename T>
void matrix_crs<T>::sort_inds(void) {

   //struct sort_pair {
   //   unsigned col;
   //   T v;
   //};

   //vector<sort_pair> pair(n); // worst case is n entries in a row

   //unsigned ind = 0;
   //for (unsigned i=0; i < m; ++i) {
   // 
   //   // Accumulate pairs of (column index,val) for row i
   //   ind = 0;
   //   pair.resize(n); // this shouldn't cause any reallocation
   //   for (unsigned j=row_ptr[i]; j < row_ptr[i+1]; ++j) {
   //      pair[ind].col = col_ind[j];
   //      pair[ind].v   = val[j];
   //      ind += 1; 
   //   }
   //   if ( ind > 0 ) pair.resize(ind);
   //   else continue; // nothing to sort!

   //   // Sort the pair
   //   sort(pair.begin(), pair.end(),
   //         [&] (const sort_pair& lhs, const sort_pair& rhs) -> bool {
   //         return (lhs.col < rhs.col);
   //      });

   //   //  Place the sorted pairs back in the matrix
   //   ind = 0;
   //   for (unsigned j=row_ptr[i]; j < row_ptr[i+1]; ++j) {
   //      col_ind[j] = pair[ind].col;
   //      val[j]     = pair[ind].v;
   //      ind += 1; 
   //   }
   //}

   // https://github.com/JuliaLang/julia/blob/master/base/sparse/sparsematrix.jl
   // http://stackoverflow.com/questions/17074324/how-can-i-sort-two-vectors-in-the-same-way-with-criteria-that-uses-only-one-of
   //unsigned ind = 0;
   vector<unsigned> index(n, 0), col_ind_wrk(n, 0); // this isn't exactly
                                                    // memory efficient...
   vector<T> val_wrk(n, 0);

   for (unsigned i = 0; i < m; ++i) {

      unsigned num = row_ptr[i+1]-row_ptr[i];
      // easy cases
      if ( num <= 1 ) continue; // zero or one elem row
      else if ( num == 2 ) {
         unsigned f = row_ptr[i], s = f+1;
         if ( col_ind[f] > col_ind[s] ) { // swap the two elements
            unsigned utmp = col_ind[s];
            col_ind[s] = col_ind[f];
            col_ind[f] = utmp;

            T Ttmp = val[s];
            val[s] = val[f];
            val[f] = Ttmp;
         }
         continue;
      }

      copy(col_ind.begin() + row_ptr[i], col_ind.begin() + row_ptr[i+1],
            col_ind_wrk.begin());
      copy(val.begin() + row_ptr[i], val.begin() + row_ptr[i+1],
            val_wrk.begin());

      iota(index.begin(), index.begin() + num, 0);
      sort(index.begin(), index.begin() + num, 
            [&] (unsigned pi, unsigned pj) -> bool 
            { return *(col_ind_wrk.begin() + pi) 
            < *(col_ind_wrk.begin() + pj); });

      unsigned ind = 0;
      for (unsigned jp = row_ptr[i]; jp < row_ptr[i+1]; ++jp) {
         col_ind[jp] = col_ind_wrk[index[ind]];
         val[jp] = val_wrk[index[ind]];

         ++ind;
      }
   }

   //cout << "C (post-sort) = " << endl;
   //cout << *this << endl;
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


// matrix transpose
// inspired by
// http://www.cise.ufl.edu/research/sparse/CSparse/CSparse/Source/cs_transpose.c
template<typename T>
void transpose(matrix_crs<T>& dest, const matrix_crs<T>& src) {

   dest.m = src.n;
   dest.n = src.m;
   dest.row_ptr.resize(src.n+1);
   dest.col_ind.resize(src.col_ind.size());
   dest.val.resize(src.val.size());

   fill(dest.row_ptr.begin(), dest.row_ptr.end(), 0);
   for (unsigned p = 0; p < src.val.size(); ++p) {
      dest.row_ptr[src.col_ind[p]+1] += 1;
   }

   unsigned cumsum = dest.row_ptr[0];
   for (unsigned i = 1; i < dest.row_ptr.size(); ++i) {
      cumsum += dest.row_ptr[i];
      dest.row_ptr[i] = cumsum;
   }

   vector<unsigned> w = dest.row_ptr; // copy
  
   unsigned q;
   for (unsigned i = 0; i < src.m; ++i) {
      for (unsigned p = src.row_ptr[i]; p < src.row_ptr[i+1]; ++p) {
         q = w[src.col_ind[p]]++; // increment subtleties :)
         dest.col_ind[q] = i;
         dest.val[q] = src.val[p];
      }
   }
}

template<typename T>
matrix_crs<T> transpose(const matrix_crs<T>& A) {
   
   matrix_crs<T> At;
   transpose(At,A);
   return At;
   //return matrix_crs<T>(At.row_ptr, At.col_ind, At.val, At.m, At.n, 1);
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


// matrix-matrix product
// See "Two Fast Algorithms for Sparse Matrices: Multiplication
//      and Permuted Transposition" - Gustavson, 1978 
template<typename T>
matrix_crs<T> operator*(const matrix_crs<T>& A, const matrix_crs<T>& B) {

   unsigned Am = A.m, An = A.n, Bm = B.m, Bn = B.n;

   if ( An != Bm ) {
      cerr << "error: matrix_crs.cpp:operator*: matrix-matrix product "
           << "dimension mismatch." << endl;
      exit(-1);
   }

   vector<unsigned> xb(Bn,-1);
   vector<T> x(Bn,0.);

   // guess size for C variables
   unsigned Cnnz = MIN(Am*Bn, A.val.size() + B.val.size());

   vector<unsigned> Crow_ptr, Ccol_ind;
   vector<T> Cval;
   Crow_ptr.resize(Am+1);
   Ccol_ind.resize(Cnnz);
   Cval.resize(Cnnz);
  
   // Do the Gustavson algorithm
   unsigned ip = 0;

   for (unsigned i=0; i < Am; ++i) {
      Crow_ptr[i] = ip;

      // check if we need more storage
      if ( ip > Cnnz - Bm) {
         Ccol_ind.resize(Cnnz+MIN(Cnnz,Am));
         Cval.resize(Cnnz+MIN(Cnnz,Am));
         Cnnz = Cval.size();
      }
     
      // Do the work of mat-mat
      for (unsigned jp = A.row_ptr[i]; jp < A.row_ptr[i+1]; ++jp) {
         unsigned j = A.col_ind[jp];
         T nzA = A.val[jp];

         for (unsigned kp = B.row_ptr[j]; kp < B.row_ptr[j+1]; ++kp) {
            unsigned k = B.col_ind[kp];
            T nzC = nzA * B.val[kp];


            if ( xb[k] != i ) {
               Ccol_ind[ip] = k;
               ip += 1;
               xb[k] = i;
               x[k] = nzC;
            }

            else {
               x[k] += nzC;
            }
         }
      }

      // copy result over to Cval
      for (unsigned vp = Crow_ptr[i]; vp < ip; ++vp) {
         Cval[vp] = x[Ccol_ind[vp]];
      }
 
   }
   Crow_ptr[Am] = ip;

   Ccol_ind.resize(*(Crow_ptr.end()-1));
   Cval.resize(*(Crow_ptr.end()-1));

   // C is not guaranteed to be sorted
   matrix_crs<T> C(Crow_ptr, Ccol_ind, Cval, Am, Bn, 1); // Build CRS directly

   //cout << "C (pre-sort) = " << endl;
   //cout << C << endl;
   //C.print_full();

   C.sort_inds();

   return C;
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

   return matrix_coo<T>(new_row_ind, col_ind, val, m, n, 0);
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
matrix_crs<T> zeros(unsigned m, unsigned n) {
   vector<unsigned> row_ptr(m+1,0);
   vector<unsigned> col_ind(0);
   vector<T> val(0);
   
   return matrix_crs<T>(row_ptr, col_ind, val, m, n, 1);
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


template<typename T>
valarray<T> diag(const matrix_crs<T>& A) {
   //unsigned N = MIN(A.m, A.n);
   //valarray<T> d(0.,N);

   //for (unsigned i = 0; i < N; ++i) {
   //   for (unsigned jp = A.row_ptr[i]; jp < A.row_ptr[i+1]; ++jp) {
   //      unsigned j = A.col_ind[jp];
   //      if ( j < i ) continue;
   //      if ( j == i) d[i] = A.val[jp];
   //      if ( j > i ) break; // diagonal of A is zero; d initialized to zero
   //   }
   //}

   // This assumes that the matrix A is in sorted CRS format
   unsigned N = MIN(A.m, A.n);
   valarray<T> d(0.,N);

   for (unsigned i = 0; i < N; ++i) {
      auto it_jp = lower_bound(A.col_ind.begin() + A.row_ptr[i],
            A.col_ind.begin() + A.row_ptr[i+1], i);
      unsigned j = *it_jp;
      if ( it_jp != A.col_ind.begin() + A.row_ptr[i+1] && j == i ) d[i] = A.val[distance(A.col_ind.begin(),it_jp)];
      else break; // diagonal of A is zero; d initialized to zero
   }
 
   return d;
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

template<typename T>
void matrix_crs<T>::print_full(void) const {
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
template void transpose<double>(matrix_crs<double>& dest,
                                              const matrix_crs<double>& src);
template matrix_crs<double> transpose<double>(const matrix_crs<double>& A);

template matrix_crs<double> eye_crs<double>(unsigned, unsigned);
template matrix_crs<double> zeros<double>(unsigned, unsigned);
template matrix_crs<double> kron<double>(const matrix_crs<double>&,
                                         const matrix_crs<double>&);

template valarray<double> diag<double>(const matrix_crs<double>&);

template matrix_crs<double> operator*(const matrix_crs<double>&,const double&);
template matrix_crs<double> operator*(const double&,const matrix_crs<double>&);
template matrix_crs<double> operator/(const matrix_crs<double>&,const double&);

template matrix_crs<double> operator+(const matrix_crs<double>&,
                                      const matrix_crs<double>&);
template matrix_crs<double> operator-(const matrix_crs<double>&,
                                      const matrix_crs<double>&);

template valarray<double> operator*(const matrix_crs<double>&,
                                      const valarray<double>&);

template matrix_crs<double> operator*(const matrix_crs<double>&,
                                      const matrix_crs<double>&);

