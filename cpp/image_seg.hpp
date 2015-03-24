// image_seg.hpp
//
// Spring 2015
// James Folberth, Nathan Heavner, Rachel Tutmaher

#ifndef _IMAGE_SEG_HPP_
#define _IMAGE_SEG_HPP_ 

#include <iostream>
#include <iomanip>
#include <valarray>
#include <vector>
#include <string>
#include <opencv2/opencv.hpp> // We use OpenCV for image I/O and Mat 

#include "matrix_coo.hpp"
#include "matrix_crs.hpp"
#include "classical_solvers.hpp"
#include "multigrid.hpp"
#include "utils.hpp"

using namespace std;

////////////////////////
// Image Segmentation //
////////////////////////
// {{{

// Class to store all image segmentation variables for a specific level
// The idea here is to create an STL vector of levels
class image_level {
   public:
      matrix_crs<double> A; // coupling matrix
      matrix_crs<double> S; // variance matrix
      matrix_crs<double> L; // weighted boundary length matrix
      matrix_crs<double> V; // area matrix
      matrix_crs<double> G; // boundary length matrix

      valarray<double> Gamma; // saliency vector

      // construction
      image_level() = default;

      // destruction, move, copy
      image_level(const image_level& ) = default;
      image_level& operator=(const image_level& ) = default;
      image_level(image_level&& ) = default;
      image_level& operator=(image_level&& ) = default;
      ~image_level() = default;

};


// Class to store the image segmentation parameters
// These parameters are global, in the sense that the shouldn't vary across 
// the levels
class seg_params {
   public:
      double alpha;     // top-level intensity scaling factor
      double alpha_til; // coarse-level intensity rescaling factor
      double beta;      // coarse-level variance rescaling factor
      double theta;     // coarsening strength threshold
      double gamma;     // saliency threshold
      double d1;        // sharpening threshold
      double sigma;     // segment detection threshold level
      double rho;       // variance rescaling detection threshold level

      // construction
      seg_params() = default;
      seg_params(double init_alpha, double init_alpha_til, double init_beta,
            double init_theta, double init_gamma, double init_d1,
            double init_sigma, double init_rho) {
         alpha = init_alpha;
         alpha_til = init_alpha_til;
         beta = init_beta;
         theta = init_theta;
         gamma = init_gamma;
         d1 = init_d1;
         sigma = init_sigma;
         rho = init_rho;
      }

      // destruction, move, copy
      seg_params(const seg_params& ) = default;
      seg_params& operator=(const seg_params& ) = default;
      seg_params(seg_params&& ) = default;
      seg_params& operator=(seg_params&& ) = default;
      ~seg_params() = default;

};

// }}}


///////////////////////
// Image Seg Drivers //
///////////////////////
// {{{

void image_seg(const string& img_filename, const seg_params& params);

void image_seg(const string& img_filename);


// }}}



///////////////////
// OpenCV things //
///////////////////
// {{{

cv::Mat load_image(const string& img_filename);

valarray<double> image_to_intensity(cv::Mat img);

// }}}



#endif
