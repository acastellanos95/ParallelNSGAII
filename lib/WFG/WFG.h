//
// Created by https://git.cse.msu.edu/talukde1/cse890/tree/4ca3b3323f010d6fbf2638d3420808ec3af78ce2/onsga2r/src/wfg-suite
//

#ifndef NSGAII_LIB_WFG_WFG_H_
#define NSGAII_LIB_WFG_WFG_H_

#include <vector>
using std::vector;

class WFG {
  const double PI_value = 3.1415926535897932384626433832795;
 public:
  bool ArgsOK(const vector<double> &z, const int k, const int M);
  //** True if all elements of "x" are in [0,1], and m is in [1, x.size()]. ***
  bool shape_args_ok(const vector<double> &x, const int m);
  bool vector_in_01(const vector<double> &x);
  double correct_to_01(const double &a, const double &epsilon);
  vector<double> WFG_normalise_z(const vector<double> &z);
  //** Construct a vector with the elements v[head], ..., v[tail-1]. **********
  vector<double> subvector(const vector<double> &v, const int head, const int tail);
  void evaluate_wfg1(vector<double> &x, vector<double> &obj);
  void evaluate_wfg2(vector<double> &x, vector<double> &obj);
  vector<double> WFG1(const std::vector<double> &z, const int k, const int M);
  vector<double> WFG2(const std::vector<double> &z, const int k, const int M);
  //** t1 from WFG1. **********************************************************
  std::vector<double> WFG1_t1(const std::vector<double> &y, const int k);
  //** t2 from WFG1. **********************************************************
  std::vector<double> WFG1_t2(const std::vector<double> &y, const int k);
  //** t3 from WFG1. **********************************************************
  std::vector<double> WFG1_t3(const std::vector<double> &y);
  //** t4 from WFG1. **********************************************************
  std::vector<double> WFG1_t4(const std::vector<double> &y, const int k, const int M);
  //** t2 from WFG2. **********************************************************
  std::vector<double> WFG2_t2(const std::vector<double> &y, const int k);
  //** t3 from WFG2. Effectively as per WFG4, t2. *****************************
  std::vector<double> WFG2_t3(const std::vector<double> &y, const int k, const int M);
  //** The polynomial bias transformation function. ***************************
  double b_poly(const double &y, const double &alpha);
  //** The flat region bias transformation function. **************************
  double b_flat(const double &y, const double &A, const double &B, const double &C);
  //** The linear shift transformation function. ******************************
  double s_linear(const double &y, const double &A);
  //** The weighted sum reduction transformation function. ********************
  double r_sum(const std::vector<double> &y, const std::vector<double> &w);
  //** The non-separable reduction transformation function. *******************
  double r_nonsep(const std::vector<double> &y, const int A);
  //** Given the last transition vector, get the fitness values for WFG1. *****
  std::vector<double> WFG1_shape(const std::vector<double> &t_p);
  //** Given the last transition vector, get the fitness values for WFG2. *****
  std::vector<double> WFG2_shape(const std::vector<double> &t_p);
  //** Construct a vector of length M-1, with values "1,0,0,..." if ***********
  //** "degenerate" is true, otherwise with values "1,1,1,..." if   ***********
  //** "degenerate" is false.                                       ***********
  vector<short> WFG_create_A(const int M, const bool degenerate);
  //** Given the vector "x" (the last value of which is the sole distance ****
  //** parameter), and the shape function results in "h", calculate the   ****
  //** scaled fitness values for a WFG problem.                           ****
  vector<double> WFG_calculate_f(const vector<double> &x, const vector<double> &h);
  vector<double> calculate_f(const double &D,
                             const vector<double> &x,
                             const vector<double> &h,
                             const vector<double> &S);
  vector<double> calculate_x(const vector<double> &t_p, const vector<short> &A);
  //** The convex shape function. (m is indexed from 1.) **********************
  double convex(const std::vector<double> &x, const int m);
  //** The mixed convex/concave shape function. *******************************
  double mixed(const std::vector<double> &x, const int A, const double &alpha);
  //** The disconnected shape function. ***************************************
  double disc(const std::vector<double> &x, const int A, const double &alpha, const double &beta);
};

#endif //NSGAII_LIB_WFG_WFG_H_
