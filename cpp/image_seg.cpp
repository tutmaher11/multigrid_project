// image_seg.cpp
//
// Spring 2015
// James Folberth, Nathan Heavner, Rachel Tutmaher

#include "image_seg.hpp"

////////////////////////
// Image Segmentation //
////////////////////////
// {{{

// First Level
//////////////
void build_first_level(const valarray<double>& I, const seg_params& params,
      vector<image_level>::iterator it) {

   // connectivity
   it->A = build_A1(I, params);

   // variance (init to zeros)
   it->S = zeros<double>(I.size(), I.size());
 
   // weighted boundary length
   build_L(it, params);

   // area 
   it->V = it->A;
   fill(it->V.val.begin(), it->V.val.end(), 1.);

   // boundary length
   build_G(it, params);

   // Build initial saliency vector
   it->Gamma = diag<double>(it->L);
   it->Gamma /= diag<double>(it->G);

}


matrix_crs<double> build_A1(const valarray<double>& I,
      const seg_params& params) {
// {{{
// Build A by rows, since we store in CRS
// The ordering of I for a 5x5 image (e.g. square.png) is
//
// 0   1   2   3   4
// 5   6   7   8   9
// 10  11  12  13  14
// 15  16  17  18  19
// 20  21  22  23  24
//

   unsigned n = params.n;
   double a = params.alpha;

   vector<unsigned> row_ptr, col_ind;
   vector<double> val;
   row_ptr.resize(I.size()+1);
   col_ind.reserve(4*2+4*(n-2)*3+(n-2)*(n-2)*4);
   val.reserve(4*2+4*(n-2)*3+(n-2)*(n-2)*4);

   if ( _DEBUG_ >= 2) {
      cout << "debug 2: building A1 (nnz = " 
           << 4*2+4*(n-2)*3+(n-2)*(n-2)*4 << ")" << endl;
   }

   row_ptr[0] = 0;
   for (unsigned ind = 0; ind < I.size(); ++ind) {
      row_ptr[ind+1] = row_ptr[ind];

      if ( 0 <= ind && ind < n ) { // on top
         if ( ind % n == 0 ) { // on left side
            // we must have ind == 0, top-left
            if ( _DEBUG_ >= 2) {
               cout << "debug 2: ind = " << ind 
                    << " top row, left side" << endl;
            }

            ++row_ptr[ind+1];
            col_ind.push_back(ind+1);
            val.push_back(exp(-a*abs(I[ind]-I[ind+1])));

            ++row_ptr[ind+1];
            col_ind.push_back(ind+n);
            val.push_back(exp(-a*abs(I[ind]-I[ind+n])));
         }

         else if ( (ind + 1) % n == 0 ) { // on right side
            // we must have ind == n-1, top-right
            if ( _DEBUG_ >= 2) {
               cout << "debug 2: ind = " << ind 
                    << " top row, right side" << endl;
            }
 
            ++row_ptr[ind+1];
            col_ind.push_back(ind-1);
            val.push_back(exp(-a*abs(I[ind]-I[ind-1])));

            ++row_ptr[ind+1];
            col_ind.push_back(ind+n);
            val.push_back(exp(-a*abs(I[ind]-I[ind+n])));
         }

         else { // in the middle
            ++row_ptr[ind+1];
            if ( _DEBUG_ >= 2) {
               cout << "debug 2: ind = " << ind 
                    << " top row, middle side" << endl;
            }
 
            col_ind.push_back(ind-1);
            val.push_back(exp(-a*abs(I[ind]-I[ind-1])));

            ++row_ptr[ind+1];
            col_ind.push_back(ind+1);
            val.push_back(exp(-a*abs(I[ind]-I[ind+1])));

            ++row_ptr[ind+1];
            col_ind.push_back(ind+n);
            val.push_back(exp(-a*abs(I[ind]-I[ind+n])));
         }
      }

      else if ( (n <= ind) && (ind < (n-1)*n) ) { // in the middle rows
         if ( ind % n == 0 ) { // on left side
            if ( _DEBUG_ >= 2) {
               cout << "debug 2: ind = " << ind 
                    << " middle row, left side" << endl;
            }
 
            ++row_ptr[ind+1];
            col_ind.push_back(ind-n);
            val.push_back(exp(-a*abs(I[ind]-I[ind-n])));

            ++row_ptr[ind+1];
            col_ind.push_back(ind+1);
            val.push_back(exp(-a*abs(I[ind]-I[ind+1])));

            ++row_ptr[ind+1];
            col_ind.push_back(ind+n);
            val.push_back(exp(-a*abs(I[ind]-I[ind+n])));
         }

         else if ( (ind + 1) % n == 0 ) { // on right side
            if ( _DEBUG_ >= 2) {
               cout << "debug 2: ind = " << ind 
                    << " middle row, right side" << endl;
            }
 
            ++row_ptr[ind+1];
            col_ind.push_back(ind-n);
            val.push_back(exp(-a*abs(I[ind]-I[ind-n])));

            ++row_ptr[ind+1];
            col_ind.push_back(ind-1);
            val.push_back(exp(-a*abs(I[ind]-I[ind-1])));

            ++row_ptr[ind+1];
            col_ind.push_back(ind+n);
            val.push_back(exp(-a*abs(I[ind]-I[ind+n])));
         }

         else { // in the middle
            if ( _DEBUG_ >= 2) {
               cout << "debug 2: ind = " << ind 
                    << " middle row, middle side" << endl;
            }
 
            ++row_ptr[ind+1];
            col_ind.push_back(ind-n);
            val.push_back(exp(-a*abs(I[ind]-I[ind-n])));

            ++row_ptr[ind+1];
            col_ind.push_back(ind-1);
            val.push_back(exp(-a*abs(I[ind]-I[ind-1])));

            ++row_ptr[ind+1];
            col_ind.push_back(ind+1);
            val.push_back(exp(-a*abs(I[ind]-I[ind+1])));

            ++row_ptr[ind+1];
            col_ind.push_back(ind+n);
            val.push_back(exp(-a*abs(I[ind]-I[ind+n])));
         }
      }

      else { // in bottom row
         if ( ind % n == 0 ) { // on left side
            if ( _DEBUG_ >= 2) {
               cout << "debug 2: ind = " << ind 
                    << " bottom row, left side" << endl;
            }
 
            ++row_ptr[ind+1];
            col_ind.push_back(ind-n);
            val.push_back(exp(-a*abs(I[ind]-I[ind-n])));

            ++row_ptr[ind+1];
            col_ind.push_back(ind+1);
            val.push_back(exp(-a*abs(I[ind]-I[ind+1])));
         }

         else if ( (ind + 1) % n == 0 ) { // on right side
            if ( _DEBUG_ >= 2) {
               cout << "debug 2: ind = " << ind 
                    << " bottom row, right side" << endl;
            }
 
            ++row_ptr[ind+1];
            col_ind.push_back(ind-n);
            val.push_back(exp(-a*abs(I[ind]-I[ind-n])));

            ++row_ptr[ind+1];
            col_ind.push_back(ind-1);
            val.push_back(exp(-a*abs(I[ind]-I[ind-1])));
         }

         else { // in the middle
            if ( _DEBUG_ >= 2) {
               cout << "debug 2: ind = " << ind 
                    << " bottom row, middle side" << endl;
            }
 
            ++row_ptr[ind+1];
            col_ind.push_back(ind-n);
            val.push_back(exp(-a*abs(I[ind]-I[ind-n])));

            ++row_ptr[ind+1];
            col_ind.push_back(ind-1);
            val.push_back(exp(-a*abs(I[ind]-I[ind-1])));

            ++row_ptr[ind+1];
            col_ind.push_back(ind+1);
            val.push_back(exp(-a*abs(I[ind]-I[ind+1])));
         }
      }
   }

   return matrix_crs<double>(row_ptr, col_ind, val, I.size(), I.size(), 1);
// }}}
}

void build_L(vector<image_level>::iterator it, 
      const seg_params& params) {
// {{{
//     {  -A_{ij}                  i != j
// L = {
//     {   \sum_{k\neq i} A_{ik}   i == j
//
// Equivalently,
// 1) L \gets -A
// 2) L_{ii} \gets -\sum_{k\neq i} L_{ik}
//
// Should be 
//    L = -A; for i in 1:size(L,1); L[i,i] -= sum(L[i,:]); end
// in julia.
//

   it->L = -1.*it->A;

   for (unsigned i = 0; i < it->L.m; ++i) {
      
      // compute full row sum
      double sum = 0.;
      for (unsigned jp = it->L.row_ptr[i]; jp < it->L.row_ptr[i+1]; ++jp) {
         sum += it->L.val[jp];
      }

      // Insert into matrix
      for (unsigned jp = it->L.row_ptr[i]; jp < it->L.row_ptr[i+1]; ++jp) {
         unsigned j = it->L.col_ind[jp];

         if ( j < i && i != it->L.m-1 ) {
            continue;
         }

         else if ( j == i ) {
            // this spot already exists in the matrix,
            // so we don't need to insert

            it->L.val[jp] -= sum;
            break;
         }

         else { // j > i or i == m-1 (last row)
            // the spot j == i doesn't already exist in the matrix, 
            // so we need to insert a new element and adjust row_ptr
            
            if ( i != it->L.m-1 ) {
               it->L.col_ind.insert(it->L.col_ind.begin() + jp, i);
               it->L.val.insert(it->L.val.begin() + jp, -sum);
            }
            else {
               it->L.col_ind.insert(it->L.col_ind.end(), i);
               it->L.val.insert(it->L.val.end(), -sum);
            }

            // adjust row pointers
            for ( unsigned row = i+1; row <= it->L.m; ++row) {
               it->L.row_ptr[row] += 1;
            }

            break;
         }
      }
   }

// }}}
}

void build_G(vector<image_level>::iterator it, 
      const seg_params& params) {
// {{{
//     {  -V_{ij}                  i != j
// G = {
//     {   \sum_{k\neq i} V_{ik}   i == j
//
// Equivalently,
// 1) G \gets -V
// 2) G_{ii} \gets -\sum_{k\neq i} G_{ik}
//
// Should be 
//    G = -V; for i in 1:size(G,1); G[i,i] -= sum(G[i,:]); end
// in julia.
//

   it->G = -1.*it->V;

   for (unsigned i = 0; i < it->G.m; ++i) {
      
      // compute full row sum
      double sum = 0.;
      for (unsigned jp = it->G.row_ptr[i]; jp < it->G.row_ptr[i+1]; ++jp) {
         sum += it->G.val[jp];
      }

      // Insert into matrix
      for (unsigned jp = it->G.row_ptr[i]; jp < it->G.row_ptr[i+1]; ++jp) {
         unsigned j = it->G.col_ind[jp];

         if ( j < i && i != it->G.m-1 ) {
            continue;
         }

         else if ( j == i ) {
            // this spot already exists in the matrix,
            // so we don't need to insert

            it->G.val[jp] -= sum;
            break;
         }

         else { // j > i or i == m-1 (last row)
            // the spot j == i doesn't already exist in the matrix, 
            // so we need to insert a new element and adjust row_ptr
            
            if ( i != it->G.m-1 ) {
               it->G.col_ind.insert(it->G.col_ind.begin() + jp, i);
               it->G.val.insert(it->G.val.begin() + jp, -sum);
            }
            else {
               it->G.col_ind.insert(it->G.col_ind.end(), i);
               it->G.val.insert(it->G.val.end(), -sum);
            }

            // adjust row pointers
            for ( unsigned row = i+1; row <= it->G.m; ++row) {
               it->G.row_ptr[row] += 1;
            }

            break;
         }
      }
   }

// }}}
}


// }}}


///////////////////////
// Image Seg Drivers //
///////////////////////
// {{{

// main driver for image segmentation
void image_seg(const string& img_filename, seg_params& params) {

   cv::Mat image;

   image = load_image(img_filename);
   valarray<double> I = image_to_intensity(image, params);
  

   // First Level stuff
   vector<image_level> levels(1);
   build_first_level(I, params, levels.begin());

}

// Overloaded drivers
/////////////////////

// use default segmentation parameters
void image_seg(const string& img_filename) {

   // (not necessarily good...) default parameters 
   // eqn 21 of Inglis et al. 2010
   seg_params params(100.,100.,100.,0.1,0.1,0.15,5.,1.);

   image_seg(img_filename, params);
}

// }}}


///////////////////
// OpenCV things //
///////////////////
// {{{

// Load image from file specified by filename into a CV Mat with uchar entries
// The image is blended to a grayscale image and sacled so the pixels
// have values between 0 and 1
// Matrices must be size n x n for the image segmentation method
cv::Mat load_image(const string& img_filename) {

   cv::Mat img;
   img = cv::imread(img_filename, 0); // load image to grayscale; uchar vals

   // Check to see if we read the image properly
   if ( img.data == NULL ) {
      cerr << "error: image_seg.cpp:load_image: OpenCV imread loaded empty "
           << "matrix.  filename = " << img_filename << endl;
      exit(-1);
   }

   // Check if the image is square
   if ( img.rows != img.cols && img.rows > 0 ) {
      cerr << "error: image_seg.cpp:load_image: the image \"" << img_filename 
           << "\" is not square." << endl;
      exit(-1);
   }

   return img;
}


// Compute the intensity vector from a matrix read in via `load_image`
// The intensity vector has values scaled from 0 to 1 (0->0, 255->1)
// The entries of the vector are the entries of the image matrix sorted
// row-wise
valarray<double> image_to_intensity(cv::Mat img, seg_params& params) {

   valarray<double> I(0., img.rows*img.cols);
   
   // Set the size of the image in the parameters class
   // We'll use it later to build A on the first level
   params.n = static_cast<unsigned>(img.rows);

   for (int i=0; i < img.rows; ++i) {

      const uchar* img_i = img.ptr<uchar>(i); // pointer to row i
      for (int j=0; j < img.cols; ++j) {
         // scale to [0,1] from range [0, 255] (uchar)
         I[img.rows*j + i] = static_cast<double>( img_i[j] )/ 255.;
      }
   }
  
   return I;
}

// }}}


int main(void) {

   image_seg("square.png");
   //image_seg("peppers.jpg");

   return 0;
}
