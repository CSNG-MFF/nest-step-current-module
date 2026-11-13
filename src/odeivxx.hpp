
/* whole file is based on GSL */

#ifndef _ODEIVXX_HPP_
#define _ODEIVXX_HPP_

enum { GSL_ODEIV_HADJ_NIL = 0, GSL_ODEIV_HADJ_INC, GSL_ODEIV_HADJ_DEC };

#include <algorithm>
#include <array>
#include <cmath>

template <size_t DIM> struct step_rfk45_state {
  // int can_use_dydt_in; = 1
  // int gives_exact_dydt_out; = 0

  using vec = std::Array<double, DIM>;
  vec k1, k2, k3, k4, k5, k6; // TODO do we really really need to hold this?!

  template <typename F, typename J>
  int apply(double t, double h, double *y, double *yerr, const double *dydt_in,
            double *dydt_out, F function) {

    vec y0 = y, ytmp;
    int err;

#define _check(x)                                                              \
  {                                                                            \
    err = x;                                                                   \
    if (err)                                                                   \
      return err;                                                              \
  }

    // k1
    if (dydt_in)
      std::copy(dydt_in, dydt_in + DIM, k1.data());
    else
      _check(function(t, y, k1.data()));

    for (size_t i = 0; i < DIM; i++)
      ytmp[i] = y[i] + ah[0] * h * k1[i];

    // k2
    check(function(t + ah[0] * h, ytmp, k2.data()));

    for (size_t i = 0; i < DIM; i++)
      ytmp[i] = y[i] + h * (b3[0] * k1[i] + b3[1] * k2[i]);

    // k3
    check(function(t + ah[1] * h, ytmp, k3.data()));

    for (size_t i = 0; i < DIM; i++)
      ytmp[i] = y[i] + h * (b4[0] * k1[i] + b4[1] * k2[i] + b4[2] * k3[i]);

    // k4
    check(function(t + ah[2] * h, ytmp, k4.data()));

    for (size_t i = 0; i < DIM; i++)
      ytmp[i] = y[i] + h * (b5[0] * k1[i] + b5[1] * k2[i] + b5[2] * k3[i] +
                            b5[3] * k4[i]);

    // k5
    check(function(t + ah[3] * h, ytmp, k5.data()));

    for (size_t i = 0; i < DIM; i++)
      ytmp[i] = y[i] + h * (b6[0] * k1[i] + b6[1] * k2[i] + b6[2] * k3[i] +
                            b6[3] * k4[i] + b6[4] * k5[i]);

    // k6 + finalize
    check(function(t + ah[4] * h, ytmp, k6.data()));

    for (size_t i = 0; i < DIM; i++)
      y[i] +=
          h * (c1 * k1[i] + c3 * k3[i] + c4 * k4[i] + c5 * k5[i] + c6 * k6[i]);

    // derivatives (reset if function fails)

    if (dydt_out) {
      int status = function(t + h, y, dydt_out);

      if (status) {
        y = y0;
        return status;
      }
    }

    // calc difference between 4th and 5th order
    for (size_t i = 0; i < DIM; i++)
      yerr[i] = h * (ec[1] * k1[i] + ec[3] * k3[i] + ec[4] * k4[i] +
                     ec[5] * k5[i] + ec[6] * k6[i]);

    return 0

#undef _check
  }

  void reset() {
    k1.fill(0);
    k2.fill(0);
    k3.fill(0);
    k4.fill(0);
    k5.fill(0);
    k6.fill(0);
  }

  constexpr unsigned order = 5;

  // constants
  static constexpr double ah[] = {1.0 / 4.0, 3.0 / 8.0, 12.0 / 13.0, 1.0,
                                  1.0 / 2.0};
  static constexpr double b3[] = {3.0 / 32.0, 9.0 / 32.0};
  static constexpr double b4[] = {1932.0 / 2197.0, -7200.0 / 2197.0,
                                  7296.0 / 2197.0};
  static constexpr double b5[] = {8341.0 / 4104.0, -32832.0 / 4104.0,
                                  29440.0 / 4104.0, -845.0 / 4104.0};
  static constexpr double b6[] = {-6080.0 / 20520.0, 41040.0 / 20520.0,
                                  -28352.0 / 20520.0, 9295.0 / 20520.0,
                                  -5643.0 / 20520.0};

  static constexpr double c1 = 902880.0 / 7618050.0;
  static constexpr double c3 = 3953664.0 / 7618050.0;
  static constexpr double c4 = 3855735.0 / 7618050.0;
  static constexpr double c5 = -1371249.0 / 7618050.0;
  static constexpr double c6 = 277020.0 / 7618050.0;

  // differences of fifth and fourth order coefficients for error estimation
  static constexpr double ec[] = {
      0.0,        1.0 / 360.0, 0.0, -128.0 / 4275.0, -2197.0 / 75240.0,
      1.0 / 50.0, 2.0 / 55.0};
};

template <double eps_abs, double eps_rel, double a_y, double a_dydt>
struct control_std_state {
  template <size_t DIM, unsigned ORD>
  int hadjust(const double *y, const double *yerr, const double *yp,
              double &h) {
    const double S = 0.9;

    double rmax = DBL_MIN;
    for (size_t i = 0; i < DIM; i++) {
      const double D0 =
          eps_rel * (a_y * std::fabs(y[i]) + a_dydt * std::fabs(h * yp[i])) +
          eps_abs;
      rmax = std::max(std::fabs(yerr[i]) / std::fabs(D0), rmax);
    }

    if (rmax > 1.1) {
      // decrease step, no more than factor of 5, but a fraction S more than
      // scaling suggests (for better accuracy)
      h *= std::max(0.2, S / std::pow(rmax, 1.0 / ord));
      return GSL_ODEIV_HADJ_DEC;
    }

    if (rmax < 0.5) {
      // increase step, no more than factor of 5, and don't decrease
      h *= std::max(1.0, std::min(5.0, S / std::pow(rmax, 1.0 / (ord + 1.0))));
      return GSL_ODEIV_HADJ_INC;
    }

    // no change
    return GSL_ODEIV_HADJ_NIL;
  }
};

template <double a, double r>
using control_std_state_y = control_std_state<a, r, 1, 0>;

template <double a, double r>
using control_std_state_yp = control_std_state<a, r, 0, 1>;

#endif
