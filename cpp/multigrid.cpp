// multigrid.cpp
//
// James Folberth
// Spring 2015

#include "multigrid.hpp"

using namespace std;

///////////
// Level //
///////////
// {{{
template<typename T>
vector<level<T>> build_levels_1d(function<matrix_crs<T>(unsigned)> build_A,
      valarray<T> f, function<matrix_crs<T>(unsigned)> build_P,
      function<matrix_crs<T>(unsigned)> build_R,
      function<void(const matrix_crs<T>&, const valarray<T>&,
         valarray<T>&, unsigned)> smoother_ip,
      unsigned L, unsigned Lx, const valarray<T>& v0, unsigned v0_level) {
   // Set up A,P,R,smoother, etc. for levels.  This is just a generic setup
   // for multigrid with A coarsened by PDE discretization.
   //
   // build_A, build_P, and build_R take exactly one argument: Lx
   // Use a lambda function to define build_A if you're using something
   // like model_problem_1d(Lx, sigma)

   // build and populate vector of levels
   vector<level<T>> levels(L);

   // A,P,R
   for (unsigned l = 0; l < L; ++l) {
      levels[l].A = build_A(Lx-l);
      levels[l].f = valarray<T>(0.,levels[l].A.n);
      levels[l].v = valarray<T>(0.,levels[l].A.n);
      levels[l].P = build_P(Lx-l);
      levels[l].R = build_R(Lx-l);
      levels[l].smoother_ip = smoother_ip;
   }

   // initial RHS and v
   levels[v0_level].f = f;
   levels[v0_level].v = v0;

   //for (auto it = levels.begin(); it != levels.end(); ++it) {
   //   //it->A.print_full();
   //   //it->P.print_full();
   //   it->R.print_full();
   //}

   return levels;
}

template<typename T>
vector<level<T>> build_levels_1d(function<matrix_crs<T>(unsigned)> build_A,
      function<matrix_crs<T>(unsigned)> build_P,
      function<matrix_crs<T>(unsigned)> build_R,
      function<void(const matrix_crs<T>&, const valarray<T>&,
         valarray<T>&, unsigned)> smoother_ip,
      unsigned L, unsigned Lx) {
   // Set up A,P,R,smoother, etc. for levels.  This is just a generic setup
   // for multigrid with A coarsened by PDE discretization.  
   //
   // Set f and v to zero vectors on all grids

   // build and populate vector of levels
   vector<level<T>> levels(L);

   // A,P,R
   for (unsigned l = 0; l < L; ++l) {
      levels[l].A = build_A(Lx-l);
      levels[l].f = valarray<T>(0.,levels[l].A.n);
      levels[l].v = valarray<T>(0.,levels[l].A.n);
      levels[l].P = build_P(Lx-l);
      levels[l].R = build_R(Lx-l);
      levels[l].smoother_ip = smoother_ip;
   }

   return levels;
}


template<typename T>
vector<level<T>> build_levels_2d(function<matrix_crs<T>(unsigned,unsigned)> build_A,
      valarray<T> f, function<matrix_crs<T>(unsigned,unsigned)> build_P,
      function<matrix_crs<T>(unsigned,unsigned)> build_R,
      function<void(const matrix_crs<T>&, const valarray<T>&,
         valarray<T>&, unsigned)> smoother_ip,
      unsigned L, unsigned Lx, unsigned Ly, const valarray<T>& v0, 
      unsigned v0_level) {
   // Set up A,P,R,smoother, etc. for levels.  This is just a generic setup
   // for multigrid with A coarsened by PDE discretization.
   //
   // build_A, build_P, and build_R take exactly one argument: Lx
   // Use a lambda function to define build_A if you're using something
   // like model_problem_1d(Lx, sigma)

   // build and populate vector of levels
   vector<level<T>> levels(L);

   // build/copy A,P,R, etc.
   for (unsigned l = 0; l < L; ++l) {
      levels[l].A = build_A(Lx-l, Ly-l);
      levels[l].f = valarray<T>(0.,levels[l].A.n);
      levels[l].v = valarray<T>(0.,levels[l].A.n);
      levels[l].P = build_P(Lx-l, Ly-l);
      levels[l].R = build_R(Lx-l, Ly-l);
      levels[l].smoother_ip = smoother_ip;
   }

   // initial RHS and v
   levels[v0_level].f = f;
   levels[v0_level].v = v0;

   return levels;
}

template<typename T>
vector<level<T>> build_levels_2d(function<matrix_crs<T>(unsigned,unsigned)> build_A,
      function<matrix_crs<T>(unsigned,unsigned)> build_P,
      function<matrix_crs<T>(unsigned,unsigned)> build_R,
      function<void(const matrix_crs<T>&, const valarray<T>&,
         valarray<T>&, unsigned)> smoother_ip,
      unsigned L, unsigned Lx, unsigned Ly) {
   // Set up A,P,R,smoother, etc. for levels.  This is just a generic setup
   // for multigrid with A coarsened by PDE discretization.
   //
   // Set f and v to zero vectors on all levels

   // build and populate vector of levels
   vector<level<T>> levels(L);

   // build/copy A,P,R, etc.
   for (unsigned l = 0; l < L; ++l) {
      levels[l].A = build_A(Lx-l, Ly-l);
      levels[l].f = valarray<T>(0.,levels[l].A.n);
      levels[l].v = valarray<T>(0.,levels[l].A.n);
      levels[l].P = build_P(Lx-l, Ly-l);
      levels[l].R = build_R(Lx-l, Ly-l);
      levels[l].smoother_ip = smoother_ip;
   }

   return levels;
}
// }}}


////////////
// Cycles //
////////////
// {{{
template<typename T>
void mucycle(vector<level<T>>& levels, typename vector<level<T>>::iterator it,
      unsigned nu1, unsigned nu2, unsigned mu) {
   // Perform a (nu1,nu2,mu) Mu-cycle on the levels                        
   // the iterator it points to our current position in the levels vector.

   valarray<T> resid(0.,it->f.size());

   // pre-smooth
   it->smoother_ip(it->A, it->f, it->v, nu1);

   if ( it == levels.end()-1 ) { // if we're on the coarsest grid
      // smooth the dick out of it
      // TODO direct solve
      it->smoother_ip(it->A, it->f, it->v, 1000);
      //resid = it->f - (it->A)*(it->v);
      //cout << "Direct solve h-resid (" << it->v.size() << ") = " << sqrt(1./(it->A.n+1))*norm(resid,2) << endl;
   } 

   else { // we're not on the coarsest grid

      // prepare to coarsen
      resid = it->f;
      resid -= (it->A)*(it->v);

      //print_vector((it->f) - (it->A)*(it->v)); // TODO why doesn't this work?
      next(it)->f = (it->R)*resid;
      next(it)->v *= 0.;

      // move to coarser grid
      for (unsigned m = 0; m < mu; ++m) {
         mucycle(levels, next(it), nu1, nu2, mu);
      }

      // correct this grid using coarse grid
      it->v += (next(it)->P)*(next(it)->v);
   }

   // post-smooth
   it->smoother_ip(it->A, it->f, it->v, nu2);
   return;
}

template<typename T>
void vcycle(vector<level<T>>& levels, typename vector<level<T>>::iterator it,
      unsigned nu1, unsigned nu2) {
   // Perform a (nu1,nu2) V-cycle on the levels
   // A V cycle is just a Mu-cycle with mu=1
   mucycle(levels, it, nu1, nu2, 1);
}

template<typename T>
void fmg(vector<level<T>>& levels, unsigned nu1, unsigned nu2,
      unsigned num_vcycles) {
   // Perform (nu1,nu2) full multigrid with num_vcycles vcycles at each 
   // interior stage

   auto it = levels.end()-1; // this is the last element in the vector

   // smooth the dick out of it
   // TODO direct solve
   it->smoother_ip(it->A, it->f, it->v, 1000);

   // Iterate through up through grids
   // We will end up going to the top grid and doing num_vcycles vcycles there
   for (; it != levels.begin(); --it) {
      
      // Prolong up to upper level
      (it-1)->v = (it->P)*(it->v);

      // do num_vcycles vcycles
      for (unsigned nu = 0; nu < num_vcycles; ++nu) {
         vcycle(levels, it-1, nu1, nu2);
      }
   }
}



// }}}


/////////////////////////
// Intergrid Operators //
/////////////////////////
// {{{
// 1D
template<typename T>
matrix_crs<T> operator_1d_interp_lin(unsigned l) {
// form the linear operator P that performs linear interpolation from grid
// l to grid l+1.  Number of points on the initial grid is 2^l-1

   unsigned m = pow(2,l+1)-1, n = pow(2,l)-1;

   // Build vectors for the i,j,value tuples
   // There are three entries per column; there are n columns
   vector<unsigned> row_ind(3*n), col_ind(3*n);
   vector<T> val(3*n);

   for (unsigned col=0; col < n; ++col) {
      row_ind[3*col]   = 2*col;
      row_ind[3*col+1] = 2*col+1;
      row_ind[3*col+2] = 2*col+2;

      col_ind[3*col]   = col;
      col_ind[3*col+1] = col;
      col_ind[3*col+2] = col;

      val[3*col]   = static_cast<T>(0.5);
      val[3*col+1] = static_cast<T>(1);
      val[3*col+2] = static_cast<T>(0.5);
   }
 
   // build CRS matrix and return
   return matrix_crs<T>(row_ind, col_ind, val, m, n);
}

template<typename T>
matrix_crs<T> operator_1d_restrict_inj(unsigned l) {
// form the linear operator R that restricts from grid l to grid l-1 using
// injection.  Number of grid points on the initial grid is 2^l-1

   unsigned m = pow(2,l-1)-1, n = pow(2,l)-1;

   // Build vectors for the i,j,value tuples
   // There is only one element per row of the matrix
   vector<unsigned> row_ind(m), col_ind(m);
   vector<T> val(m);

   for (unsigned row=0; row < m; ++row) {
      row_ind[row] = row;
      col_ind[row] = 2*row+1;
      val[row] = static_cast<T>(1.);
   }

   // build CRS matrix and return
   return matrix_crs<T>(row_ind, col_ind, val, m, n);
}

template<typename T>
matrix_crs<T> operator_1d_restrict_full(unsigned l) {
// form the linear operator R that restricts from grid l to grid l-1 using
// full weighting.  Number of grid points on the initial grid is 2^l-1

   unsigned m = pow(2,l-1)-1, n = pow(2,l)-1;

   // Build vectors for the i,j,value tuples
   // There are three elements per row of the matrix
   vector<unsigned> row_ind(3*m), col_ind(3*m);
   vector<T> val(3*m);

   for (unsigned row=0; row < m; ++row) {
      row_ind[3*row] = row;
      row_ind[3*row+1] = row;
      row_ind[3*row+2] = row;

      col_ind[3*row] = 2*row;
      col_ind[3*row+1] = 2*row+1;
      col_ind[3*row+2] = 2*row+2;

      val[3*row] = static_cast<T>(0.25);
      val[3*row+1] = static_cast<T>(0.5);
      val[3*row+2] = static_cast<T>(0.25);
   }

   // build CRS matrix and return
   return matrix_crs<T>(row_ind, col_ind, val, m, n);
}

// 2D
template<typename T>
matrix_crs<T> operator_2d_interp_lin(unsigned lx, unsigned ly) {
// form the linear operator P that performs linear interpolation from grid
// (lx,ly) to grid (lx+1,ly+1).  Number of points on the initial grid is 
// (2^lx-1)*(2^ly-1)

   unsigned m = (pow(2,lx+1)-1)*(pow(2,ly+1)-1), // size of new mat
            n = (pow(2,lx)-1)*(pow(2,ly)-1),
            nx = pow(2,lx)-1, ny = pow(2,ly)-1, // size of old grid
            nxp = pow(2,lx+1)-1;// nyp = pow(2,ly+1)-1; // size of new grid

   // Build vectors for the i,j,value tuples
   // There are nine entries per column (9 pt stencil); there are n columns
   vector<unsigned> row_ind, col_ind;
   vector<T> val;
   row_ind.reserve(9*n);
   col_ind.reserve(9*n);
   val.reserve(9*n);

   // row and col are the row and col of the grid, not the matrix
   for (unsigned row=0; row < ny; ++row) {
      for (unsigned col=0; col < nx; ++col) {
         // NW
         col_ind.push_back(nx*row+col);
         row_ind.push_back(nxp*(2*row)+(2*col));
         val.push_back(static_cast<T>(0.25));

         // N
         col_ind.push_back(nx*row+col);
         row_ind.push_back(nxp*(2*row)+(2*col+1));
         val.push_back(static_cast<T>(0.5));

         // NE
         col_ind.push_back(nx*row+col);
         row_ind.push_back(nxp*(2*row)+(2*col+2));
         val.push_back(static_cast<T>(0.25));


         // W
         col_ind.push_back(nx*row+col);
         row_ind.push_back(nxp*(2*row+1)+(2*col));
         val.push_back(static_cast<T>(0.5));

         // Center
         col_ind.push_back(nx*row+col);
         row_ind.push_back(nxp*(2*row+1)+(2*col+1));
         val.push_back(static_cast<T>(1.));

         // E
         col_ind.push_back(nx*row+col);
         row_ind.push_back(nxp*(2*row+1)+(2*col+2));
         val.push_back(static_cast<T>(0.5));


         // SW
         col_ind.push_back(nx*row+col);
         row_ind.push_back(nxp*(2*row+2)+(2*col));
         val.push_back(static_cast<T>(0.25));

         // S
         col_ind.push_back(nx*row+col);
         row_ind.push_back(nxp*(2*row+2)+(2*col+1));
         val.push_back(static_cast<T>(0.5));

         // SE
         col_ind.push_back(nx*row+col);
         row_ind.push_back(nxp*(2*row+2)+(2*col+2));
         val.push_back(static_cast<T>(0.25));

      }
   }
 
   // build CRS matrix and return
   return matrix_crs<T>(row_ind, col_ind, val, m, n);
}


template<typename T>
matrix_crs<T> operator_2d_restrict_inj(unsigned lx, unsigned ly) {
// form the linear operator R that restricts from grid (lx,ly) to grid 
// (lx-1,ly-1) using injection.  Number of grid points on the initial grid is
// (2^lx-1)*(2^ly-1)

   unsigned m = (pow(2,lx-1)-1)*(pow(2,ly-1)-1), // size of new mat
            n = (pow(2,lx)-1)*(pow(2,ly)-1),
            nx = pow(2,lx)-1, ny = pow(2,ly)-1, // size of old grid
            nxp = pow(2,lx-1)-1; //nyp = pow(2,ly-1)-1; // size of new grid


   // Build vectors for the i,j,value tuples
   // There is only one element per row of the matrix
   vector<unsigned> row_ind, col_ind;
   vector<T> val;
   row_ind.reserve(m);
   col_ind.reserve(m);
   val.reserve(m);

   // (row,col) are the position on the original grid, not the matrix
   for (unsigned row = 1; row < ny; row+=2) {
      for (unsigned col = 1; col < nx; col+=2) {
         //cout << "(i,j) = (" << row << "," << col << ")" << endl;
         row_ind.push_back(nxp*(row-1)/2 + (col-1)/2);
         col_ind.push_back(nx*row+col);
         val.push_back(static_cast<T>(1.));
      }
   }

   // build CRS matrix and return
   return matrix_crs<T>(row_ind, col_ind, val, m, n);
}

template<typename T>
matrix_crs<T> operator_2d_restrict_full(unsigned lx, unsigned ly) {
// form the linear operator R that restricts from grid (lx,ly) to grid
// (lx-1,ly-1) using full weighting.  Number of grid points on the initial
// grid is (2^lx-1)*(2^ly-1)

   unsigned m = (pow(2,lx-1)-1)*(pow(2,ly-1)-1), // size of new mat
            n = (pow(2,lx)-1)*(pow(2,ly)-1),
            nx = pow(2,lx)-1, ny = pow(2,ly)-1, // size of old grid
            nxp = pow(2,lx-1)-1;// nyp = pow(2,ly-1)-1; // size of new grid

   // Build vectors for the i,j,value tuples
   // There are nine entries per point on the new grid; there are m points on
   // the new grid
   vector<unsigned> row_ind, col_ind;
   vector<T> val;
   row_ind.reserve(9*m);
   col_ind.reserve(9*m);
   val.reserve(9*m);

   // row and col are the row and col of the grid, not the matrix
   unsigned new_center_row, new_center_col;
   for (unsigned row=1; row < ny; row+=2) {
      for (unsigned col=1; col < nx; col+=2) {
         new_center_row = (row-1)/2;
         new_center_col = (col-1)/2;

         // {{{ Debugging
         //cout << "(i,j) = (" << row << "," << col << ")" << endl;
         //cout << "new (i,j) = (" << new_center_row << "," 
         //     << new_center_col << ")" << endl;
         //
         //cout << "using (i,j) = (" << row-1 << "," << col-1 << ")" << endl;
         //cout << "using (i,j) = (" << row-1 << "," << col-0 << ")" << endl;
         //cout << "using (i,j) = (" << row-1 << "," << col+1 << ")" << endl;
         //cout << "using (i,j) = (" << row-0 << "," << col-1 << ")" << endl;
         //cout << "using (i,j) = (" << row-0 << "," << col-0 << ")" << endl;
         //cout << "using (i,j) = (" << row-0 << "," << col+1 << ")" << endl;
         //cout << "using (i,j) = (" << row+1 << "," << col-1 << ")" << endl;
         //cout << "using (i,j) = (" << row+1 << "," << col-0 << ")" << endl;
         //cout << "using (i,j) = (" << row+1 << "," << col+1 << ")" << endl;
         // }}}

         // NW
         row_ind.push_back(nxp*new_center_row+new_center_col);
         col_ind.push_back(nx*(row-1)+(col-1));
         val.push_back(static_cast<T>(1./16.));

         // N
         row_ind.push_back(nxp*new_center_row+new_center_col);
         col_ind.push_back(nx*(row-1)+(col));
         val.push_back(static_cast<T>(1./8.));

         // NE
         row_ind.push_back(nxp*new_center_row+new_center_col);
         col_ind.push_back(nx*(row-1)+(col+1));
         val.push_back(static_cast<T>(1./16.));

         // W
         row_ind.push_back(nxp*new_center_row+new_center_col);
         col_ind.push_back(nx*(row)+(col-1));
         val.push_back(static_cast<T>(1./8.));

         // Center
         row_ind.push_back(nxp*new_center_row+new_center_col);
         col_ind.push_back(nx*(row)+(col));
         val.push_back(static_cast<T>(1./4.));

         // E
         row_ind.push_back(nxp*new_center_row+new_center_col);
         col_ind.push_back(nx*(row)+(col+1));
         val.push_back(static_cast<T>(1./8.));

         // SW
         row_ind.push_back(nxp*new_center_row+new_center_col);
         col_ind.push_back(nx*(row+1)+(col-1));
         val.push_back(static_cast<T>(1./16.));

         // S
         row_ind.push_back(nxp*new_center_row+new_center_col);
         col_ind.push_back(nx*(row+1)+(col));
         val.push_back(static_cast<T>(1./8.));

         // SE
         row_ind.push_back(nxp*new_center_row+new_center_col);
         col_ind.push_back(nx*(row+1)+(col+1));
         val.push_back(static_cast<T>(1./16.));

      }
   }
 
   // build CRS matrix and return
   return matrix_crs<T>(row_ind, col_ind, val, m, n);
}


// }}}


// Force instatiation
// double
// Level
template vector<level<double>> build_levels_1d(
      function<matrix_crs<double>(unsigned)> build_A, valarray<double> f,
      function<matrix_crs<double>(unsigned)> build_P,
      function<matrix_crs<double>(unsigned)> build_R,
      function<void(const matrix_crs<double>&, 
         const valarray<double>&, valarray<double>&, unsigned)>,
      unsigned L, unsigned Lx, const valarray<double>& v0, unsigned v0_level);

template vector<level<double>> build_levels_1d(
      function<matrix_crs<double>(unsigned)> build_A,
      function<matrix_crs<double>(unsigned)> build_P,
      function<matrix_crs<double>(unsigned)> build_R,
      function<void(const matrix_crs<double>&, 
         const valarray<double>&, valarray<double>&, unsigned)>,
      unsigned L, unsigned Lx);

template vector<level<double>> build_levels_2d(
      function<matrix_crs<double>(unsigned,unsigned)> build_A,
      valarray<double> f,
      function<matrix_crs<double>(unsigned,unsigned)> build_P,
      function<matrix_crs<double>(unsigned,unsigned)> build_R,
      function<void(const matrix_crs<double>&, 
         const valarray<double>&, valarray<double>&, unsigned)>,
      unsigned L, unsigned Lx, unsigned Ly, const valarray<double>& v0,
      unsigned v0_level);

template vector<level<double>> build_levels_2d(
      function<matrix_crs<double>(unsigned,unsigned)> build_A,
      function<matrix_crs<double>(unsigned,unsigned)> build_P,
      function<matrix_crs<double>(unsigned,unsigned)> build_R,
      function<void(const matrix_crs<double>&, 
         const valarray<double>&, valarray<double>&, unsigned)>,
      unsigned L, unsigned Lx, unsigned Ly);

// Cycles
template void mucycle(vector<level<double>>& levels, 
      vector<level<double>>::iterator, unsigned nu1, unsigned nu2,
      unsigned mu);

template void vcycle(vector<level<double>>& levels, 
      vector<level<double>>::iterator, unsigned nu1, unsigned nu2);

template void fmg(vector<level<double>>& levels, unsigned nu1, unsigned nu2,
      unsigned num_vcycles);

// Operators
template matrix_crs<double> operator_1d_interp_lin<double>(unsigned);
template matrix_crs<double> operator_1d_restrict_inj<double>(unsigned);
template matrix_crs<double> operator_1d_restrict_full<double>(unsigned);

template matrix_crs<double> operator_2d_interp_lin<double>(unsigned,
      unsigned);
template matrix_crs<double> operator_2d_restrict_inj<double>(unsigned,
      unsigned);
template matrix_crs<double> operator_2d_restrict_full<double>(unsigned,
      unsigned);


