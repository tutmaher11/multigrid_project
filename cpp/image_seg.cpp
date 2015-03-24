// image_seg.cpp
//
// Spring 2015
// James Folberth, Nathan Heavner, Rachel Tutmaher

#include "image_seg.hpp"

////////////////////////
// Image Segmentation //
////////////////////////
// {{{
// }}}


///////////////////////
// Image Seg Drivers //
///////////////////////
// {{{

// main driver for image segmentation
void image_seg(const string& img_filename, const seg_params& params) {

   cv::Mat image;

   image = load_image(img_filename);

   valarray<double> I = image_to_intensity(image);
   //print_vector(I);


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
valarray<double> image_to_intensity(cv::Mat img) {

   valarray<double> I(0., img.rows*img.cols);

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
