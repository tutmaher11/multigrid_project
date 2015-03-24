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
   //print_vector(I);

   matrix_crs<double> A1 = build_A1(I, params);
   //A1.print_full();
   //cout << A1 << endl;

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
