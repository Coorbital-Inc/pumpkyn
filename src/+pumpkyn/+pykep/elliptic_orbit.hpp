#ifndef ASTROCPP_ELLIPTIC_ORBIT_HPP
#define ASTROCPP_ELLIPTIC_ORBIT_HPP

#include "state_vectors.hpp"

namespace astro_cpp {
class EllipticOrbit {
public:
  EllipticOrbit(double std_grav_param, double semi_major_axis,
		double eccentricity,
		double inclination,	   // in radians
		double longitude_of_AN,	   // in radians
		double arg_of_periapsis,   // in radians
		double mean_anomaly_at_t0, // in radians
		double t0);
  void print() const;
  double period() const { return period_value; }
  StateVectors state_vectors_from_time(double t) const;

private:
  double mu, a, e, i, Om, om, M0, t0;

  // Some pre-computed quantities
  double speed_root, // sqrt(mu / a)
      mean_motion,   // sqrt(mu / a^3)
      period_value,  // orbital period
      e_root,	     // sqrt(1 - e^2)
      a_e_root,	     // a * e_root
      orb_to_ecl[6]; // transformation matrix
		     // from the orbital basis to
		     // the ecliptic one
  double E_from_M(double M) const;
};
} // namespace astro_cpp

#endif // ASTROCPP_ELLIPTIC_ORBIT_HPP
