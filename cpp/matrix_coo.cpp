// matrix_coo.cpp
//
// James Folberth
// Spring 2015

#include "matrix_coo.hpp"

using namespace std;

//////////////////////////////////
// Construction and destruction //
//////////////////////////////////
template<typename T>
matrix_coo<T>::matrix_coo(
      vector<unsigned>& init_row_ind,
      vector<unsigned>& init_col_ind,
      vector<T>& init_val,
      size_t init_m, size_t init_n) {

   row_ind = init_row_ind;
   col_ind = init_col_ind;
   val = init_val;

   if (val.size() == 0) {
      if (init_m == 0 && init_n == 0) {
         cerr << "error: matrix_coo: Can't construct an empty (zero) matrix "
              << "without specifying m and n!" << endl;
         exit(-1);
      }
      m = init_m; n = init_n;
      return;
   }

   unsigned max_row = 1+*max_element(row_ind.begin(), row_ind.end());
   unsigned max_col = 1+*max_element(col_ind.begin(), col_ind.end());
   m = (init_m < max_row) ? max_row : init_m;
   n = (init_n < max_col) ? max_col : init_n;

   assert(row_ind.size() == col_ind.size());
   assert(col_ind.size() == val.size());

   this->clean();
}


// Insertion sort
// O(n^2) worst case, O(n*k) for sorted w/ "bandwidth" k
// done in place, O(1) extra memory
template<typename T>
void insertion_sort(vector<unsigned>& rind, vector<unsigned>& cind,
      vector<T>& v) {

   auto compare = [&] (const unsigned lhs_r, const unsigned lhs_c,
      const unsigned rhs_r, const unsigned rhs_c) -> bool {
      if ( lhs_r < rhs_r ) {
         return 1;
      }
      else if ( lhs_r == rhs_r ) {
         if ( lhs_c < rhs_c ) return 1;
         else return 0; // ignore repeats
      }
      else {
         return 0;
      }
   };
   
   unsigned rtemp, ctemp;
   T vtemp;
   unsigned i,j;
   for (i = 1; i < rind.size(); ++i) {
      rtemp = rind[i];
      ctemp = cind[i];
      vtemp = v[i];
      j = i;

      // rotate
      while ( j > 0 && compare(rtemp, ctemp, rind[j-1], cind[j-1]) ) {
         rind[j] = rind[j-1];
         cind[j] = cind[j-1];
         v[j] = v[j-1];
         j -= 1;
      }

      // insert
      rind[j] = rtemp;
      cind[j] = ctemp;
      v[j] = vtemp;
   }             
}               


// TODO should accept a combine function (default to lambda add) like
// Julia's CSC
template<typename T>
void matrix_coo<T>::clean(void) {
// Reorder the row/col/vals so that they are in row-major order

   //// Insertion sort
   //// O(n^2) worst case, O(n*k) for sorted w/ "bandwidth" k
   //// done in place, O(1) extra memory
   // For 2D model problem, this is crazy slow
   //insertion_sort(row_ind, col_ind, val);
   //
   //cout << "Done sorting (" << row_ind.size() << ")" << endl;

   ////cout << *this << endl;

   //// Combine duplicates and remove zero elements
   //unsigned num_dups = 0;
   //if ( abs(val[0]) < _ELEMENT_ZERO_TOL_ ) {
   //   row_ind.erase(row_ind.begin());
   //   col_ind.erase(col_ind.begin());
   //   val.erase(val.begin());

   //   num_dups += 1;
   //}

   //for (unsigned i = 1; i < row_ind.size(); ++i) {

   //   if ( abs(val[i]) < _ELEMENT_ZERO_TOL_ ) {
   //      row_ind.erase(row_ind.begin() + i);
   //      col_ind.erase(col_ind.begin() + i);
   //      val.erase(val.begin() + i);

   //      num_dups += 1;
   //   }

   //   // if duplicate
   //   if ( col_ind[i] == col_ind[i-1] && row_ind[i] == row_ind[i-1] ) {

   //      // combine then erase duplicate
   //      val[i-1] += val[i];
   //      row_ind.erase(row_ind.begin() + i);
   //      col_ind.erase(col_ind.begin() + i);
   //      val.erase(val.begin() + i);
   //      
   //      num_dups += 1;
   //   }
   //}

   //if (_DEBUG_ >= 1 && num_dups > 0) {
   //   cerr << "warning: matrix_coo:clean: duplicate entries found or "
   //        << "zeros; combining entries by adding and removing zero entries"
   //        << endl;
   //}

   //// shrink container to fit size, which may cause reallocation
   //row_ind.shrink_to_fit();
   //col_ind.shrink_to_fit();
   //val.shrink_to_fit();


   // TODO this is slow
   // Idea: make a class with pointers to row_ind, col_ind, val where the iterator is a "facade" that points to a tuple of (i,j,v) or something.  I don't know how to do this
   // Idea: sort (i,j,n), where n = 0:row_ind.size()-1; then sort val using ordering of n.  This will save a bit of memory usage
   // OLD SORT
   //
   // Copy row inds, col inds, vals to vector of (i,j,val) tuples
   // Then sort vector of tuples with stdlib sort
   // Then overwrite values and remove duplicates
   // I don't know how fast this is
   
   struct sort_tuple {
      unsigned i;
      unsigned j;
      T val;
   };

   vector<sort_tuple> sort_me;
   sort_me.resize(row_ind.size());
   for (unsigned i=0; i<sort_me.size(); ++i) {
      sort_me[i].i = row_ind[i];
      sort_me[i].j = col_ind[i];
      sort_me[i].val = val[i];
   }

   sort(sort_me.begin(), sort_me.end(), 
         [&] (const sort_tuple& lhs, const sort_tuple& rhs) -> bool {
            if (lhs.i < rhs.i) {
               return 1;
            }
            else if (lhs.i == rhs.i) {
               if (lhs.j < rhs.j) return 1;
               else return 0; // ignore repeats
            }
            else return 0;
         });

   //cout << "Done sorting (" << row_ind.size() << ")" << endl;
   
   // assign to class member and deal with duplicate entries by adding
   unsigned old_row=-1, old_col=-1, num_dups=0;
   for (unsigned i=0; i<sort_me.size(); ++i) {

      if ( abs(sort_me[i].val) < _ELEMENT_ZERO_TOL_ ) {
         num_dups += 1;
         continue;
      }

      // if not a duplicate (case i=0 always passes, since it's first)
      if (old_row != sort_me[i].i || old_col != sort_me[i].j) {

         row_ind[i-num_dups] = sort_me[i].i;
         col_ind[i-num_dups] = sort_me[i].j;
         val[i-num_dups] = sort_me[i].val;

         old_row = sort_me[i].i;
         old_col = sort_me[i].j;
      }
      
      else {
         num_dups += 1;
         val[i-num_dups] += sort_me[i].val; // TODO combine function
         continue;
      }
   }

   if (_DEBUG_ >= 1 && num_dups > 0) {
      cerr << "warning: matrix_coo:clean: duplicate entries found or "
           << "zeros; combining entries by adding and removing zero entries"
           << endl;
   }

   row_ind.resize(row_ind.size()-num_dups);
   col_ind.resize(col_ind.size()-num_dups);
   val.resize(val.size()-num_dups);
   // END OLD SORT
}



// TODO maybe move to seperate file?
// special forms
template<typename T>
matrix_coo<T> eye_coo(unsigned m, unsigned n) {
   vector<unsigned> rind(MIN(m,n)), cind(MIN(m,n));
   vector<T> val(MIN(m,n), static_cast<T>(1)); // fill val with 1s on construct

   for (unsigned i=0; i < rind.size(); ++i) {
      rind[i] = i;
      cind[i] = i;
   }

   return matrix_coo<T>(rind, cind, val, m, n);
}


////////////////
// Operations //
////////////////
// scalar
template<typename T>
matrix_coo<T>& matrix_coo<T>::operator*=(const T& value) {
   for (auto it=val.begin(); it != val.end(); ++it) {
      *it *= value;
   }
   return *this;
}

template<typename T>
inline matrix_coo<T> operator*(const matrix_coo<T>& lhs, const T& rhs) {
   return matrix_coo<T>(lhs) *= rhs;
}

template<typename T>
inline matrix_coo<T> operator*(const T& lhs, const matrix_coo<T>& rhs) {
   return matrix_coo<T>(rhs) *= lhs;
}


template<typename T>
matrix_coo<T>& matrix_coo<T>::operator/=(const T& value) {
   for (auto it=val.begin(); it != val.end(); ++it) {
      *it /= value;
   }
   return *this;
}

template<typename T>
inline matrix_coo<T> operator/(const matrix_coo<T>& lhs, const T& rhs) {
   return matrix_coo<T>(lhs) /= rhs;
}

// matrix add/sub
// *this += B;
template<typename T>
matrix_coo<T>& matrix_coo<T>::operator+=(const matrix_coo<T>& B) {
   //unsigned this_ind = 0, B_ind = 0;;

   // check sizes
   if ( m != B.m || n != B.n ) {
      cerr << "error: matrix_coo:+=: matrix sizes do not match - ("
           << m << "," << n << ") vs. ("
           << B.m << "," << B.n << ")" << endl;
      exit(-1);
   }

   // One (stupid!) option is to just append the values of B to *this and 
   // then call this->clean() to sort and add.  Don't know how fast
   // this will be though.  Shouldn't really matter, since CRS will be 
   // the main type.
 
   row_ind.insert(row_ind.begin(), B.row_ind.begin(), B.row_ind.end());
   col_ind.insert(col_ind.begin(), B.col_ind.begin(), B.col_ind.end());
   val.insert(val.begin(), B.val.begin(), B.val.end());
   this->clean();

   return *this;

   // {{{
   //// *this is empty
   //if ( val.size() == 0 ) {
   //   row_ind = B.row_ind;
   //   col_ind = B.col_ind;
   //   val = B.val;
   //   
   //   return *this;
   //}
  
   //// loop through entries of B
   //while ( B_ind < B.val.size() ) {

   //   // need to append entries to *this once we're past the last one in *this
   //   if ( this_ind == val.size() ) {
   //      row_ind.push_back(B.row_ind[B_ind]);
   //      col_ind.push_back(B.col_ind[B_ind]);
   //      val.push_back(B.val[B_ind]);
   //      
   //      ++B_ind;
   //      continue;
   //   }

   //   // entry in row before *this has anything
   //   if ( B.row_ind[B_ind] < row_ind[this_ind] ) {
   //      row_ind.insert(row_ind.begin() + this_ind, B.row_ind[B_ind]);
   //      col_ind.insert(col_ind.begin() + this_ind, B.col_ind[B_ind]);
   //      val.insert(val.begin() + this_ind, B.val[B_ind]);

   //      ++B_ind;
   //      continue;
   //   }

   //   // entry in B in same row as *this or next row
   //   else {
   //      if ( B.row_ind[B_ind] == row_ind[this_ind] ) {

   //         // make this_ind point to the next entry in *this past the
   //         // entry in B pointed to by B_ind
   //         while ( this_ind < val.size() && 
   //                 row_ind[this_ind] == B.row_ind[B_ind] && 
   //                 col_ind[this_ind] < B.col_ind[B_ind] )
   //            ++this_ind;

   //         // this_ind points to the next row
   //         // repeat the outer iteration
   //         if ( this_ind == val.size() || 
   //              row_ind[this_ind] > B.row_ind[B_ind] ) {
   //            continue;
   //         }

   //         // combine entries
   //         else if ( B.col_ind[B_ind] == col_ind[this_ind] ) {
   //            val[this_ind] += B.val[B_ind];
   //
   //            // if val below _ELEMENT_ZERO_TOL_, remove it
   //            if ( abs(val[this_ind]) < _ELEMENT_ZERO_TOL_ ) {
   //               row_ind.erase(row_ind.begin()+ this_ind);
   //               col_ind.erase(col_ind.begin()+ this_ind);
   //               val.erase(val.begin()+ this_ind);
   //            }

   //            ++B_ind;
   //            continue;
   //         }

   //         // insert entry of B into *this
   //         else if ( B.col_ind[B_ind] < col_ind[this_ind] ) {
   //            row_ind.insert(row_ind.begin() + this_ind, B.row_ind[B_ind]);
   //            col_ind.insert(col_ind.begin() + this_ind, B.col_ind[B_ind]);
   //            val.insert(val.begin() + this_ind, B.val[B_ind]);

   //            ++B_ind;
   //            continue;
   //         }
   //        
   //         // if we reach here, we must be out of elements in *this
   //         // that is, this_ind == val.size()
   //         // continue loop to reach section that handles this
   //         else { 
   //            row_ind.push_back(B.row_ind[B_ind]);
   //            col_ind.push_back(B.col_ind[B_ind]);
   //            val.push_back(B.val[B_ind]);

   //            ++B_ind;
   //            continue;
   //         }
   //      }

   //      // B points past the row this_ind points to; move this_ind and 
   //      // repeat the loop or just append entry of B
   //      else if ( B.row_ind[B_ind] > row_ind[this_ind] ) {

   //         while ( this_ind < val.size() && 
   //                 row_ind[this_ind] < B.row_ind[B_ind] )
   //            ++this_ind;

   //         // past last element of *this; append entry of B
   //         if ( this_ind == val.size() ) {
   //            row_ind.push_back(B.row_ind[B_ind]);
   //            col_ind.push_back(B.col_ind[B_ind]);
   //            val.push_back(B.val[B_ind]);

   //            ++B_ind;
   //            continue;
   //         }

   //         // now try to repeat the loop
   //      }
   //   }
   //}
   // 
   //return *this;
   // }}}
}

template<typename T>
matrix_coo<T> operator+(const matrix_coo<T>& lhs, const matrix_coo<T>& rhs) {
   return matrix_coo<T>(lhs) += rhs;
}

template<typename T>
matrix_coo<T>& matrix_coo<T>::operator-=(const matrix_coo<T>& B) {
   return (*this += -1.0*B);
}

template<typename T>
matrix_coo<T> operator-(const matrix_coo<T>& lhs, const matrix_coo<T>& rhs) {
   return matrix_coo<T>(lhs) -= rhs;
}


/////////////////////
// Type conversion //
/////////////////////
template<typename T>
matrix_crs<T> matrix_coo<T>::to_crs(void) {
   return matrix_crs<T>(row_ind, col_ind, val, m, n);
}


////////////
// Output //
////////////
template<typename T>
void matrix_coo<T>::print_full(void) {
   this->to_crs().print_full();
}


// Some template functions are defined in the header
// I couldn't get things to work any other way

// Force instantiation for specific types
// double
template class matrix_coo<double>;
template matrix_coo<double> eye_coo<double>(unsigned, unsigned);

template matrix_coo<double> operator*(const matrix_coo<double>&,const double&);
template matrix_coo<double> operator*(const double&,const matrix_coo<double>&);
template matrix_coo<double> operator/(const matrix_coo<double>&,const double&);

template matrix_coo<double> operator+(const matrix_coo<double>&,
                                      const matrix_coo<double>&);
template matrix_coo<double> operator-(const matrix_coo<double>&,
                                      const matrix_coo<double>&);

