// model_problem.hpp
//
// James Folberth
// Spring 2015

#ifndef _MODEL_PROBLEMS_HPP_
#define _MODEL_PROBLEMS_HPP_

#include <vector>

#include "matrix_crs.hpp"

using namespace std;

template<typename T>
matrix_crs<T> model_problem_1d(unsigned L, T sigma = 0.);

template<typename T>
matrix_crs<T> model_problem_2d(unsigned Lx, unsigned Ly, T sigma = 0.);

#endif
