/*****************************************************************************
 *   Copyright (C) 2004-2018 The pykep development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *                                                                           *
 *   https://gitter.im/esa/pykep                                             *
 *   https://github.com/esa/pykep                                            *
 *                                                                           *
 *   act@esa.int                                                             *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.               *
 *****************************************************************************/

#include "lambert_original.hpp"

#include "const.hpp"
#include <iostream>

namespace astro_cpp {

/// Constructor
/** Constructs and solves a Lambert problem.
 *
 * \param[in] R1 first cartesian position
 * \param[in] R2 second cartesian position
 * \param[in] tof time of flight
 * \param[in] mu gravity parameter
 * \param[in] cw when 1 a retrograde orbit is assumed
 * \param[in] multi_revs maximum number of multirevolutions to compute
 */
lambert_original::lambert_original(const Vector3D &r1, const Vector3D &r2,
				   const double &tof, const double &mu,
				   const int &cw, const int &multi_revs)
    : m_r1(r1), m_r2(r2), m_tof(tof), m_mu(mu), m_has_converged(true),
      m_multi_revs(multi_revs) {
  // 1 - Getting lambda and T
  m_c = (r1 - r2).mag();
  double R1 = r1.mag();
  double R2 = r2.mag();
  m_s = (m_c + R1 + R2) / 2.0;

  Vector3D ir1 = r1.unit_vector(), ir2 = r2.unit_vector(),
	   ih = Vector3D::cross_product(ir1, ir2).unit_vector(),
	   it1 = Vector3D::cross_product(ih, ir1),
	   it2 = Vector3D::cross_product(ih, ir2);

  double lambda2 = 1.0 - m_c / m_s;
  m_lambda = sqrt(lambda2);

  if ((ih.z < 0.0) ^ cw) // Transfer angle is larger than 180 degrees as seen
			 // from abive the z axis
  {
    m_lambda = -m_lambda;
    it1.invert();
    it2.invert();
  }

  double lambda3 = m_lambda * lambda2;
  double T = sqrt(2.0 * m_mu / m_s / m_s / m_s) * m_tof;

  // 2 - We now have lambda, T and we will find all x
  // 2.1 - Let us first detect the maximum number of revolutions for which there
  // exists a solution
  m_Nmax = static_cast<int>(T / math::PI);
  double T00 = acos(m_lambda) + m_lambda * sqrt(1.0 - lambda2);
  double T0 = (T00 + m_Nmax * math::PI);
  double T1 = 2.0 / 3.0 * (1.0 - lambda3), DT = 0.0, DDT = 0.0, DDDT = 0.0;
  if ((m_Nmax <= multi_revs) && (m_Nmax > 0)) {
    if (T < T0) { // We use Halley iterations to find xM and TM
      // std::cout << "adjust revolutions Halley\n";
      int it = 0;
      double err = 1.0;
      double T_min = T0;
      double x_old = 0.0, x_new = 0.0;
      while (1) {
	dTdx(DT, DDT, DDDT, x_old, T_min);
	if (DT != 0.0) {
	  x_new = x_old - DT * DDT / (DDT * DDT - DT * DDDT / 2.0);
	}
	err = fabs(x_old - x_new);
	if ((err < 1e-13) || (it > 12)) {
	  break;
	}
	x2tof(T_min, x_new, m_Nmax);
	x_old = x_new;
	it++;
      }
      if (T_min > T) {
	m_Nmax -= 1;
      }
    }
  }
  // We exit this if clause with Mmax being the maximum number of revolutions
  // for which there exists a solution. We crop it to m_multi_revs
  m_Nmax = std::min(m_multi_revs, m_Nmax);

  // 2.2 We now allocate the memory for the output variables
  m_v1.resize(m_Nmax * 2 + 1);
  m_v2.resize(m_Nmax * 2 + 1);
  m_iters.resize(m_Nmax * 2 + 1);
  m_x.resize(m_Nmax * 2 + 1);

  // 3 - We may now find all solutions in x,y
  // 3.1 0 rev solution
  // 3.1.1 initial guess
  if (T >= T00) {
    m_x[0] = -(T - T00) / (T - T00 + 4);
  } else if (T <= T1) {
    m_x[0] = T1 * (T1 - T) / (2.0 / 5.0 * (1 - lambda2 * lambda3) * T) + 1;
  } else {
    m_x[0] = pow((T / T00), 0.69314718055994529 / log(T1 / T00)) - 1.0;
  }
  // 3.1.2 Householder iterations
  m_iters[0] = householder(T, m_x[0], 0.0, 1e-5, 15);
  // 3.2 multi rev solutions
  double tmp;
  for (int i = 1; i < m_Nmax + 1; ++i) {
    // 3.2.1 left Householder iterations
    tmp = pow((i * math::PI + math::PI) / (8.0 * T), 2.0 / 3.0);
    m_x[2 * i - 1] = (tmp - 1) / (tmp + 1);
    m_iters[2 * i - 1] = householder(T, m_x[2 * i - 1], i, 1e-8, 15);
    // 3.2.1 right Householder iterations
    tmp = pow((8.0 * T) / (i * math::PI), 2.0 / 3.0);
    m_x[2 * i] = (tmp - 1) / (tmp + 1);
    m_iters[2 * i] = householder(T, m_x[2 * i], i, 1e-8, 15);
  }

  // 4 - For each found x value we reconstruct the terminal velocities
  double gamma = sqrt(m_mu * m_s / 2.0);
  double rho = (R1 - R2) / m_c;
  double sigma = sqrt(1 - rho * rho);
  double vr1, vt1, vr2, vt2, y;
  for (size_t i = 0; i < m_x.size(); ++i) {
    y = sqrt(1.0 - lambda2 + lambda2 * m_x[i] * m_x[i]);
    vr1 =
	gamma * ((m_lambda * y - m_x[i]) - rho * (m_lambda * y + m_x[i])) / R1;
    vr2 =
	-gamma * ((m_lambda * y - m_x[i]) + rho * (m_lambda * y + m_x[i])) / R2;
    double vt = gamma * sigma * (y + m_lambda * m_x[i]);
    vt1 = vt / R1;
    vt2 = vt / R2;
    m_v1[i] = vr1 * ir1 + vt1 * it1;
    m_v2[i] = vr2 * ir2 + vt2 * it2;
  }
}

int lambert_original::householder(const double T, double &x0, const int N,
				  const double eps, const int iter_max) {
  int it = 0;
  double err = 1.0;
  double xnew = 0.0;
  double tof = 0.0, delta = 0.0, DT = 0.0, DDT = 0.0, DDDT = 0.0;
  while ((err > eps) && (it < iter_max)) {
    // std::cout << "x = " << x0 << '\n';
    x2tof(tof, x0, N);
    // std::cout << "T = " << tof << '\n';
    dTdx(DT, DDT, DDDT, x0, tof);
    // std::cout << "DT = " << DT << "; DDT = " << DDT << "; DDDT = " << DDDT <<
    // "\n\n";
    delta = tof - T;
    double DT2 = DT * DT;
    xnew = x0 - delta * (DT2 - delta * DDT / 2.0) /
		    (DT * (DT2 - delta * DDT) + DDDT * delta * delta / 6.0);
    err = fabs(x0 - xnew);
    x0 = xnew;
    it++;
  }
  return it;
}

void lambert_original::dTdx(double &DT, double &DDT, double &DDDT,
			    const double x, const double T) {
  double l2 = m_lambda * m_lambda;
  double l3 = l2 * m_lambda;
  double umx2 = 1.0 - x * x;
  double y = sqrt(1.0 - l2 * umx2);
  double y2 = y * y;
  double y3 = y2 * y;

  DT = 1.0 / umx2 * (3.0 * T * x - 2.0 + 2.0 * l3 * x / y);
  DDT = 1.0 / umx2 * (3.0 * T + 5.0 * x * DT + 2.0 * (1.0 - l2) * l3 / y3);
  DDDT = 1.0 / umx2 *
	 (7.0 * x * DDT + 8.0 * DT - 6.0 * (1.0 - l2) * l2 * l3 * x / y3 / y2);

  // std::cout << "DT = " << DT << "; DDT = " << DDT << "; DDDT = " << DDDT <<
  // "\n\n";
}

void lambert_original::x2tof2(double &tof, const double x, const int N) {
  double a = 1.0 / (1.0 - x * x);
  if (a > 0) // ellipse
  {
    double alfa = 2.0 * acos(x);
    double beta = 2.0 * asin(sqrt(m_lambda * m_lambda / a));
    if (m_lambda < 0.0)
      beta = -beta;
    tof = ((a * sqrt(a) *
	    ((alfa - sin(alfa)) - (beta - sin(beta)) + 2.0 * math::PI * N)) /
	   2.0);
  } else {
    double alfa = 2.0 * std::acosh(x);
    double beta = 2.0 * std::asinh(sqrt(-m_lambda * m_lambda / a));
    if (m_lambda < 0.0)
      beta = -beta;
    tof = (-a * sqrt(-a) * ((beta - sinh(beta)) - (alfa - sinh(alfa))) / 2.0);
  }
}

void lambert_original::x2tof(double &tof, const double x, const int N) {
  double battin = 0.01;
  double lagrange = 0.2;
  double dist = fabs(x - 1);
  if (dist < lagrange && dist > battin) { // We use Lagrange tof expression
    // std::cout << "Lagrange" << '\n';
    x2tof2(tof, x, N);
    return;
  }
  double K = m_lambda * m_lambda;
  double E = x * x - 1.0;
  double rho = fabs(E);
  double z = sqrt(1 + K * E);
  if (dist < battin) { // We use Battin series tof expression
    // std::cout << "Battin" << '\n';
    double eta = z - m_lambda * x;
    double S1 = 0.5 * (1.0 - m_lambda - x * eta);
    double Q = hypergeometricF(S1, 1e-11);
    Q = 4.0 / 3.0 * Q;
    tof = (eta * eta * eta * Q + 4.0 * m_lambda * eta) / 2.0 +
	  N * math::PI / pow(rho, 1.5);
    return;
  } else { // We use Lancaster tof expresion
    // std::cout << "Lancaster" << '\n';
    double y = sqrt(rho);
    double g = x * z - m_lambda * E;
    double d = 0.0;
    if (E < 0) {
      double l = acos(g);
      d = N * math::PI + l;
    } else {
      double f = y * (z - m_lambda * x);
      d = log(f + g);
    }
    tof = (x - m_lambda * z - d / y) / E;
    return;
  }
}

double lambert_original::hypergeometricF(double z, double tol) {
  double Sj = 1.0;
  double Cj = 1.0;
  double err = 1.0;
  double Cj1 = 0.0;
  double Sj1 = 0.0;
  int j = 0;
  while (err > tol) {
    Cj1 = Cj * (3.0 + j) * (1.0 + j) / (2.5 + j) * z / (j + 1);
    Sj1 = Sj + Cj1;
    err = fabs(Cj1);
    Sj = Sj1;
    Cj = Cj1;
    j = j + 1;
  }
  return Sj;
}

/// Gets velocity at r1
/**
 *
 * \return an std::vector containing 3-d arrays with the cartesian components of
 * the velocities at r1 for all 2N_max+1 solutions
 */
const std::vector<Vector3D> &lambert_original::get_v1() const { return m_v1; }

/// Gets velocity at r2
/**
 *
 * \return an std::vector containing 3-d arrays with the cartesian components of
 * the velocities at r2 for all 2N_max+1 solutions
 */
const std::vector<Vector3D> &lambert_original::get_v2() const { return m_v2; }

/// Gets r1
/**
 *
 * \return a 3-d array with the cartesian components of r1
 */
const Vector3D &lambert_original::get_r1() const { return m_r1; }

/// Gets r2
/**
 *
 * \return a 3-d array with the cartesian components of r2
 */
const Vector3D &lambert_original::get_r2() const { return m_r2; }

/// Gets the time of flight between r1 and r2
/**
 *
 * \return the time of flight
 */
const double &lambert_original::get_tof() const { return m_tof; }

/// Gets the x variable
/**
 * Gets the x variable for each solution found (0 revs, 1,1,2,2,3,3 .... N,N)
 *
 * \return the x variables in an std::vector
 */
const std::vector<double> &lambert_original::get_x() const { return m_x; }

/// Gets gravitational parameter
/**
 *
 * \return the gravitational parameter
 */
const double &lambert_original::get_mu() const { return m_mu; }

/// Gets number of iterations
/**
 *
 * \return an std::vector containing the iterations taken to compute each one of
 * the solutions
 */
const std::vector<int> &lambert_original::get_iters() const { return m_iters; }

/// Gets N_max
/**
 *
 * \return the maximum number of revolutions. The number of solutions to the
 * problem will be Nmax*2 +1
 */
int lambert_original::get_Nmax() const { return m_Nmax; }

} // namespace kep_toolbox
