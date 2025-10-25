#ifndef LAMBERT_ORIGINAL_HPP
#define LAMBERT_ORIGINAL_HPP

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

#include "vector3d.hpp"
#include <cmath>
#include <vector>

namespace astro_cpp {

/// Lambert Problem
/**
 * This class represent a Lambert's problem. When instantiated it assumes a
 * prograde orbit (unless otherwise stated) and evaluates all the solutions up
 * to a maximum number of multiple revolutions. After the object is instantiated
 * the solutions can be retreived using the appropriate getters. Note that the
 * number of solutions will be N_max*2 + 1, where N_max is the maximum number of
 * revolutions.
 *
 * NOTE: The class has been tested extensively via monte carlo runs checked with
 * numerical propagation. Compared to the previous Lambert Solver in the
 * keplerian_toolbox it is 1.7 times faster (on average as defined by
 * lambert_test.cpp). With respect to Gooding algorithm it is 1.3 - 1.5 times
 * faster (zero revs - multi revs). The algorithm is described in detail in the
 * publication below and its original with the author.
 *
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */

#define _USE_MATH_DEFINES

class lambert_original {
public:
  lambert_original(const Vector3D &r1, const Vector3D &r2, const double &tof,
		   const double &mu, const int &cw, const int &multi_revs);
  const std::vector<Vector3D> &get_v1() const;
  const std::vector<Vector3D> &get_v2() const;
  const Vector3D &get_r1() const;
  const Vector3D &get_r2() const;
  const double &get_tof() const;
  const double &get_mu() const;
  const std::vector<double> &get_x() const;
  const std::vector<int> &get_iters() const;
  int get_Nmax() const;

private:
  int householder(const double T, double &x0, const int N, const double eps,
		  const int itermax);
  void dTdx(double &DT, double &DDT, double &DDDT, const double x0,
	    const double tof);
  void x2tof(double &tof, const double x0, const int N);
  void x2tof2(double &tof, const double x0, const int N);
  double hypergeometricF(double z, double tol);

  const Vector3D m_r1, m_r2;
  const double m_tof;
  const double m_mu;
  std::vector<Vector3D> m_v1;
  std::vector<Vector3D> m_v2;
  std::vector<int> m_iters;
  std::vector<double> m_x;
  double m_s, m_c, m_lambda;
  int m_Nmax;
  bool m_has_converged;
  int m_multi_revs;
};
} // namespace kep_toolbox

#endif // LAMBERT_ORIGINAL_HPP
