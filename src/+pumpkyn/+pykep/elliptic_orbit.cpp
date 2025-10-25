#include "elliptic_orbit.hpp"

#include "const.hpp"
#include <cmath>
#include <iostream>

namespace astro_cpp {
EllipticOrbit::EllipticOrbit(double std_grav_param, double semi_major_axis,
			     double eccentricity,
			     double inclination,	// in radians
			     double longitude_of_AN,	// in radians
			     double arg_of_periapsis,	// in radians
			     double mean_anomaly_at_t0, // in radians
			     double t0)
    : mu{std_grav_param}, a{semi_major_axis}, e{eccentricity}, i{inclination},
      Om{longitude_of_AN}, om{arg_of_periapsis}, M0{mean_anomaly_at_t0},
      t0{t0} {
  const double inv_a = 1.0 / a;
  speed_root = std::sqrt(mu * inv_a);
  mean_motion = speed_root * inv_a;
  period_value = math::TAU / mean_motion;
  e_root = std::sqrt(1.0 - e * e);
  a_e_root = a * e_root;

  const double sin_Om = std::sin(Om), cos_Om = std::cos(Om),
	       sin_om = std::sin(om), cos_om = std::cos(om),
	       sin_i = std::sin(i), cos_i = std::cos(i),
	       cos_Om_cos_om = cos_Om * cos_om, sin_Om_sin_om = sin_Om * sin_om,
	       sin_Om_cos_om = sin_Om * cos_om, cos_Om_sin_om = cos_Om * sin_om;

  orb_to_ecl[0] = cos_Om_cos_om - sin_Om_sin_om * cos_i;
  orb_to_ecl[1] = -cos_Om_sin_om - sin_Om_cos_om * cos_i;
  orb_to_ecl[2] = sin_Om_cos_om + cos_Om_sin_om * cos_i;
  orb_to_ecl[3] = -sin_Om_sin_om + cos_Om_cos_om * cos_i;
  orb_to_ecl[4] = sin_om * sin_i;
  orb_to_ecl[5] = cos_om * sin_i;
}

// Elliptic Kepler's equation solver (about 40% faster than naive implementation
// with Newton's method) Based on the paper by F. Landis Markley "Kepler
// equation solver" (1995)
double EllipticOrbit::E_from_M(double M) const {
  constexpr double alpha_mult = 1.0 / (math::PI_SQR - 6.0);

  // Wrap M to [-PI; +PI]
  M = std::fmod(M + math::PI, math::TAU);
  if (M < 0.0) {
    M += math::TAU;
  }
  M -= math::PI;

  // Solving cubic equation based on Pade approximant for sin(E) (more info in
  // the paper)
  const double alpha = (3.0 * math::PI_SQR +
			1.6 * math::PI * (math::PI - std::abs(M)) / (1.0 + e)) *
		       alpha_mult,
	       d = 3.0 * (1.0 - e) + alpha * e, M_sqr = M * M,
	       alpha_d = alpha * d, q = 2.0 * alpha_d * (1.0 - e) - M_sqr,
	       q_sqr = q * q, r = 3.0 * alpha_d * (d - 1.0 + e) * M + M_sqr * M,
	       pre_w = std::abs(r) + std::sqrt(q_sqr * q + r * r),
	       w = std::cbrt(pre_w * pre_w);

  // Finally, get our initial guess
  double E = ((2.0 * r * w) / (w * w + w * q + q_sqr) + M) / d;

  // 3rd order Householder's method
  const double sin_E = std::sin(E), cos_E = std::cos(E), f = E - e * sin_E - M,
	       f_sqr = f * f, Df = 1.0 - e * cos_E, Df_sqr = Df * Df,
	       DDf = e * sin_E, DDDf = e * cos_E;
  E = E - (6.0 * f * Df_sqr - 3.0 * f_sqr * DDf) /
	      (6.0 * Df_sqr * Df - 6.0 * f * Df * DDf + f_sqr * DDDf);

  return E;
}

StateVectors EllipticOrbit::state_vectors_from_time(double t) const {
  double M = mean_motion * (t - t0) + M0, // Mean anomaly
      E = E_from_M(M);

  // Recompute sin(E) and cos(E) for final value of E
  const double sin_E = std::sin(E), cos_E = std::cos(E);

  // Radius vector in orbital basis
  const double r_orb_x = a * (cos_E - e), r_orb_y = a_e_root * sin_E;

  // Velocity in orbital basis
  const double v_mult = speed_root / (1 - e * cos_E), v_orb_x = -v_mult * sin_E,
	       v_orb_y = v_mult * e_root * cos_E;

  // Transform to the ecliptical basis
  return {{orb_to_ecl[0] * r_orb_x + orb_to_ecl[1] * r_orb_y,
	   orb_to_ecl[2] * r_orb_x + orb_to_ecl[3] * r_orb_y,
	   orb_to_ecl[4] * r_orb_x + orb_to_ecl[5] * r_orb_y},

	  {orb_to_ecl[0] * v_orb_x + orb_to_ecl[1] * v_orb_y,
	   orb_to_ecl[2] * v_orb_x + orb_to_ecl[3] * v_orb_y,
	   orb_to_ecl[4] * v_orb_x + orb_to_ecl[5] * v_orb_y}};
}

void EllipticOrbit::print() const {
  std::cout << "\nstd. grav. param. = " << mu << "\nsemi-major axis = " << a
	    << "\neccentricity = " << e
	    << "\ninclination = " << i * math::RAD_TO_DEG << " deg"
	    << "\nlongitude of the AN = " << Om * math::RAD_TO_DEG << " deg"
	    << "\narg. of periapsis = " << om * math::RAD_TO_DEG << " deg"
	    << "\nmean anomaly at t0 = " << M0 << " rad"
	    << "\nt0 = " << t0 << '\n';
}
} // namespace astro_cpp
