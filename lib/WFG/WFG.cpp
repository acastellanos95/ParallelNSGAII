//
// https://git.cse.msu.edu/talukde1/cse890/tree/4ca3b3323f010d6fbf2638d3420808ec3af78ce2/onsga2r/src/wfg-suite
//

#include "WFG.h"

#include <vector>
#include <cassert>
#include <cmath>
using std::vector;

//** True if "k" in [1,z.size()), "M" >= 2, and "k" mod ("M"-1) == 0. *******
bool WFG::ArgsOK(const vector<double> &z, const int k, const int M) {
  const int n = static_cast< int >( z.size());

  return k >= 1 && k < n && M >= 2 && k % (M - 1) == 0;
}

//** True if all elements of "x" are in [0,1], and m is in [1, x.size()]. ***
bool WFG::shape_args_ok(const vector<double> &x, const int m) {
  const int M = static_cast< int >( x.size());

  return vector_in_01(x) && m >= 1 && m <= M;
}

bool WFG::vector_in_01(const vector<double> &x) {
  for (int i = 0; i < static_cast< int >( x.size()); i++) {
    if (x[i] < 0.0 || x[i] > 1.0) {
      return false;
    }
  }

  return true;
}

double WFG::correct_to_01(const double &a, const double &epsilon = 1.0e-10) {
  assert(epsilon >= 0.0);

  const double min = 0.0;
  const double max = 1.0;

  const double min_epsilon = min - epsilon;
  const double max_epsilon = max + epsilon;

  if (a <= min && a >= min_epsilon) {
    return min;
  } else if (a >= max && a <= max_epsilon) {
    return max;
  } else {
    return a;
  }
}

//** Reduces each paramer in "z" to the domain [0,1]. ***********************
vector<double> WFG::WFG_normalise_z(const vector<double> &z) {
  vector<double> result;

  for (int i = 0; i < static_cast<int>(z.size()); i++) {
    const double bound = 2.0 * (i + 1);

    assert(z[i] >= 0.0);
    assert(z[i] <= bound);

    result.push_back(z[i] / bound);
  }

  return result;
}

//** Construct a vector with the elements v[head], ..., v[tail-1]. **********
vector<double> WFG::subvector(const vector<double> &v, const int head, const int tail) {
  assert(head >= 0);
  assert(head < tail);
  assert(tail <= static_cast< int >( v.size()));

  vector<double> result;

  for (int i = head; i < tail; i++) {
    result.push_back(v[i]);
  }

  return result;
}

vector<double> WFG::WFG1(const vector<double> &z, const int k, const int M) {
  assert(ArgsOK(z, k, M));

  vector<double> y = WFG_normalise_z(z);

  y = WFG1_t1(y, k);
  y = WFG1_t2(y, k);
  y = WFG1_t3(y);
  y = WFG1_t4(y, k, M);

  return WFG1_shape(y);
}

vector<double> WFG::WFG2(const vector<double> &z, const int k, const int M) {
  assert(ArgsOK(z, k, M));
  assert((static_cast< int >( z.size()) - k) % 2 == 0);

  vector<double> y = WFG_normalise_z(z);

  y = WFG1_t1(y, k);
  y = WFG2_t2(y, k);
  y = WFG2_t3(y, k, M);

  return WFG2_shape(y);
}

void WFG::evaluate_wfg1(vector<double> &x, vector<double> &obj) {
  int k = 4, M = obj.size();
  vector<double> f = WFG1(x, k, M);
  for (unsigned int i = 0; i < f.size(); i++)
    obj[i] = f[i];
}

void WFG::evaluate_wfg2(vector<double> &x, vector<double> &obj) {
  int k = 4, M = obj.size();
  vector<double> f = WFG2(x, k, M);
  for (unsigned int i = 0; i < f.size(); i++)
    obj[i] = f[i];
}

vector<double> WFG::WFG1_t1(const vector<double> &y, const int k) {
  const int n = static_cast< int >( y.size());

  assert(vector_in_01(y));
  assert(k >= 1);
  assert(k < n);

  vector<double> t;

  for (int i = 0; i < k; i++) {
    t.push_back(y[i]);
  }

  for (int i = k; i < n; i++) {
    t.push_back(s_linear(y[i], 0.35));
  }

  return t;
}

vector<double> WFG::WFG1_t2(const vector<double> &y, const int k) {
  const int n = static_cast< int >( y.size());

  assert(vector_in_01(y));
  assert(k >= 1);
  assert(k < n);

  vector<double> t;

  for (int i = 0; i < k; i++) {
    t.push_back(y[i]);
  }

  for (int i = k; i < n; i++) {
    t.push_back(b_flat(y[i], 0.8, 0.75, 0.85));
  }

  return t;
}

vector<double> WFG::WFG1_t3(const vector<double> &y) {
  const int n = static_cast< int >( y.size());

  assert(vector_in_01(y));

  vector<double> t;

  for (int i = 0; i < n; i++) {
    t.push_back(b_poly(y[i], 0.02));
  }

  return t;
}

vector<double> WFG::WFG1_t4(const vector<double> &y, const int k, const int M) {
  const int n = static_cast< int >( y.size());

  assert(vector_in_01(y));
  assert(k >= 1);
  assert(k < n);
  assert(M >= 2);
  assert(k % (M - 1) == 0);

  vector<double> w;

  for (int i = 1; i <= n; i++) {
    w.push_back(2.0 * i);
  }

  vector<double> t;

  for (int i = 1; i <= M - 1; i++) {
    const int head = (i - 1) * k / (M - 1);
    const int tail = i * k / (M - 1);

    const vector<double> &y_sub = subvector(y, head, tail);
    const vector<double> &w_sub = subvector(w, head, tail);

    t.push_back(r_sum(y_sub, w_sub));
  }

  const vector<double> &y_sub = subvector(y, k, n);
  const vector<double> &w_sub = subvector(w, k, n);

  t.push_back(r_sum(y_sub, w_sub));

  return t;
}

vector<double> WFG::WFG2_t2(const vector<double> &y, const int k) {
  const int n = static_cast< int >( y.size());
  const int l = n - k;

  assert(vector_in_01(y));
  assert(k >= 1);
  assert(k < n);
  assert(l % 2 == 0);

  vector<double> t;

  for (int i = 0; i < k; i++) {
    t.push_back(y[i]);
  }

  for (int i = k + 1; i <= k + l / 2; i++) {
    const int head = k + 2 * (i - k) - 2;
    const int tail = k + 2 * (i - k);

    t.push_back(r_nonsep(subvector(y, head, tail), 2));
  }

  return t;
}

vector<double> WFG::WFG2_t3(const vector<double> &y, const int k, const int M) {
  const int n = static_cast< int >( y.size());

  assert(vector_in_01(y));
  assert(k >= 1);
  assert(k < n);
  assert(M >= 2);
  assert(k % (M - 1) == 0);

  const vector<double> w(n, 1.0);

  vector<double> t;

  for (int i = 1; i <= M - 1; i++) {
    const int head = (i - 1) * k / (M - 1);
    const int tail = i * k / (M - 1);

    const vector<double> &y_sub = subvector(y, head, tail);
    const vector<double> &w_sub = subvector(w, head, tail);

    t.push_back(r_sum(y_sub, w_sub));
  }

  const vector<double> &y_sub = subvector(y, k, n);
  const vector<double> &w_sub = subvector(w, k, n);

  t.push_back(r_sum(y_sub, w_sub));

  return t;
}

double WFG::b_poly(const double &y, const double &alpha) {
  assert(y >= 0.0);
  assert(y <= 1.0);
  assert(alpha > 0.0);
  assert(alpha != 1.0);

  return correct_to_01(std::pow(y, alpha));
}

double WFG::b_flat(const double &y, const double &A, const double &B, const double &C) {
  assert(y >= 0.0);
  assert(y <= 1.0);
  assert(A >= 0.0);
  assert(A <= 1.0);
  assert(B >= 0.0);
  assert(B <= 1.0);
  assert(C >= 0.0);
  assert(C <= 1.0);
  assert(B < C);
  assert(B != 0.0 || A == 0.0);
  assert(B != 0.0 || C != 1.0);
  assert(C != 1.0 || A == 1.0);
  assert(C != 1.0 || B != 0.0);

  const double tmp1 = std::min(0.0, floor(y - B)) * A * (B - y) / B;
  const double tmp2 = std::min(0.0, floor(C - y)) * (1.0 - A) * (y - C) / (1.0 - C);

  return correct_to_01(A + tmp1 - tmp2);
}

double WFG::s_linear(const double &y, const double &A) {
  assert(y >= 0.0);
  assert(y <= 1.0);
  assert(A > 0.0);
  assert(A < 1.0);

  return correct_to_01(fabs(y - A) / fabs(floor(A - y) + A));
}

double WFG::r_sum(const vector<double> &y, const vector<double> &w) {
  assert(y.size() != 0);
  assert(w.size() == y.size());
  assert(vector_in_01(y));

  double numerator = 0.0;
  double denominator = 0.0;

  for (int i = 0; i < static_cast< int >( y.size()); i++) {
    assert(w[i] > 0.0);

    numerator += w[i] * y[i];
    denominator += w[i];
  }

  return correct_to_01(numerator / denominator);
}

double WFG::r_nonsep(const std::vector<double> &y, const int A) {
  const int y_len = static_cast< int >( y.size());

  assert(y_len != 0);
  assert(vector_in_01(y));
  assert(A >= 1);
  assert(A <= y_len);
  assert(y.size() % A == 0);

  double numerator = 0.0;

  for (int j = 0; j < y_len; j++) {
    numerator += y[j];

    for (int k = 0; k <= A - 2; k++) {
      numerator += fabs(y[j] - y[(j + k + 1) % y_len]);
    }
  }

  const double tmp = ceil(A / 2.0);
  const double denominator = y_len * tmp * (1.0 + 2.0 * A - 2.0 * tmp) / A;

  return correct_to_01(numerator / denominator);
}

vector<double> WFG::WFG1_shape(const vector<double> &t_p) {
  assert(vector_in_01(t_p));
  assert(t_p.size() >= 2);

  const int M = static_cast< int >( t_p.size());

  const vector<short> &A = WFG_create_A(M, false);
  const vector<double> &x = calculate_x(t_p, A);

  vector<double> h;

  for (int m = 1; m <= M - 1; m++) {
    h.push_back(convex(x, m));
  }
  h.push_back(mixed(x, 5, 1.0));

  return WFG_calculate_f(x, h);
}

vector<double> WFG::WFG2_shape(const vector<double> &t_p) {
  assert(vector_in_01(t_p));
  assert(t_p.size() >= 2);

  const int M = static_cast< int >( t_p.size());

  const vector<short> &A = WFG_create_A(M, false);
  const vector<double> &x = calculate_x(t_p, A);

  vector<double> h;

  for (int m = 1; m <= M - 1; m++) {
    h.push_back(convex(x, m));
  }
  h.push_back(disc(x, 5, 1.0, 1.0));

  return WFG_calculate_f(x, h);
}

//** Construct a vector of length M-1, with values "1,0,0,..." if ***********
//** "degenerate" is true, otherwise with values "1,1,1,..." if   ***********
//** "degenerate" is false.                                       ***********
vector<short> WFG::WFG_create_A(const int M, const bool degenerate) {
  assert(M >= 2);

  if (degenerate) {
    vector<short> A(M - 1, 0);
    A[0] = 1;

    return A;
  } else {
    return vector<short>(M - 1, 1);
  }
}

vector<double> WFG::WFG_calculate_f(const vector<double> &x, const vector<double> &h) {
  assert(vector_in_01(x));
  assert(vector_in_01(h));
  assert(x.size() == h.size());

  const int M = static_cast< int >( h.size());

  vector<double> S;

  for (int m = 1; m <= M; m++) {
    S.push_back(m * 2.0);
  }

  return calculate_f(1.0, x, h, S);
}

vector<double> WFG::calculate_x(const vector<double> &t_p, const vector<short> &A) {
  assert(vector_in_01(t_p));
  assert(t_p.size() != 0);
  assert(A.size() == t_p.size() - 1);

  vector<double> result;

  for (int i = 0; i < static_cast< int >( t_p.size()) - 1; i++) {
    assert(A[i] == 0 || A[i] == 1);

    const double tmp1 = std::max<double>(t_p.back(), A[i]);
    result.push_back(tmp1 * (t_p[i] - 0.5) + 0.5);
  }

  result.push_back(t_p.back());

  return result;
}

vector<double> WFG::calculate_f(const double &D,
                                const vector<double> &x,
                                const vector<double> &h,
                                const vector<double> &S) {
  assert(D > 0.0);
  assert(vector_in_01(x));
  assert(vector_in_01(h));
  assert(x.size() == h.size());
  assert(h.size() == S.size());

  vector<double> result;

  for (int i = 0; i < static_cast< int >( h.size()); i++) {
    assert(S[i] > 0.0);

    result.push_back(D * x.back() + S[i] * h[i]);
  }

  return result;
}

double WFG::convex(const vector<double> &x, const int m) {
  assert(shape_args_ok(x, m));

  const int M = static_cast< int >( x.size());
  double result = 1.0;

  for (int i = 1; i <= M - m; i++) {
    result *= 1.0 - cos(x[i - 1] * PI_value / 2.0);
  }

  if (m != 1) {
    result *= 1.0 - sin(x[M - m] * PI_value / 2.0);
  }

  return correct_to_01(result);
}

double WFG::mixed(const vector<double> &x, const int A, const double &alpha) {
  assert(vector_in_01(x));
  assert(x.size() != 0);
  assert(A >= 1);
  assert(alpha > 0.0);

  const double tmp = 2.0 * A * PI_value;

  return correct_to_01(pow(1.0 - x[0] - cos(tmp * x[0] + PI_value / 2.0) / tmp, alpha));
}

double WFG::disc(const vector<double> &x, const int A, const double &alpha, const double &beta) {
  assert(vector_in_01(x));
  assert(x.size() != 0);
  assert(A >= 1);
  assert(alpha > 0.0);
  assert(beta > 0.0);

  const double tmp1 = A * pow(x[0], beta) * PI_value;
  return correct_to_01(1.0 - pow(x[0], alpha) * pow(cos(tmp1), 2.0));
}