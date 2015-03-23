// model_problems.cpp
//
// James Folberth
// Spring 2015

#include "model_problems.hpp"

template<typename T>
matrix_crs<T> model_problem_1d(unsigned L, T sigma) {

   unsigned m = pow(2,L)-1, n = m;
   T h2 = pow(1./static_cast<T>(m+1),2);

   vector<unsigned> rind, cind;
   vector<T> val;

   if (L < 2) {
      cerr << "error: model_problems.cpp:model_problem_1d: Must have L >= 2."
           << endl;
      exit(-1);
   }

   // preallocate enough space for building tridiagonal mxm 
   rind.resize(3*m-2); cind.resize(3*m-2); val.resize(3*m-2);

   // first row
   rind[0] = 0; rind[1] = 0;
   cind[0] = 0; cind[1] = 1;
   val[0] = 2./h2+sigma; val[1] = -1./h2;

   // interior rows
   for (unsigned i=1; i < m-1; ++i) {
      rind[3*i-1] = i;
      rind[3*i] = i;
      rind[3*i+1] = i;

      cind[3*i-1] = i-1;
      cind[3*i] = i;
      cind[3*i+1] = i+1;

      val[3*i-1] = -1./h2;
      val[3*i] = 2./h2+sigma;
      val[3*i+1] = -1./h2;
   }

   // last row
   rind[3*m-4] = m-1; rind[3*m-3] = m-1;
   cind[3*m-4] = m-2; cind[3*m-3] = m-1;
   val[3*m-4] = -1./h2; val[3*m-3] = 2./h2+sigma;

   return matrix_crs<T>(rind, cind, val, m, n);
}


template<typename T>
matrix_crs<T> model_problem_2d(unsigned Lx, unsigned Ly, T sigma) {
   // Form 2D system using 1D x and y system and Kronecker product

   unsigned m = pow(2,Lx)-1, n = pow(2,Ly)-1;
   matrix_crs<T> Ax = model_problem_1d<T>(Lx,static_cast<T>(0.0));
   matrix_crs<T> Ix = eye_crs<T>(m,m);
   matrix_crs<T> Ay = model_problem_1d<T>(Ly,static_cast<T>(0.0));
   matrix_crs<T> Iy = eye_crs<T>(n,n);

   // row-wise ordering
   return kron(Ax,Iy)+kron(Ix,Ay)+sigma*eye_crs<T>(m*n,m*n);
   // column-wise ordering
   //return kron(Iy,Ax)+kron(Ay,Ix)+sigma*eye_crs<T>(m*n,m*n);
}


// Force instantiation for specific types
template matrix_crs<double> model_problem_1d<double>(unsigned,double);
template matrix_crs<double> model_problem_2d<double>(unsigned,unsigned,double);
