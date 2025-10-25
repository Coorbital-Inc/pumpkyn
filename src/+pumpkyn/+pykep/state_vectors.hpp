#ifndef ASTROCPP_STATE_VECTORS_HPP
#define ASTROCPP_STATE_VECTORS_HPP

#include "vector3d.hpp"
#include <iostream>

namespace astro_cpp {
class StateVectors {
public:
  Vector3D r, v;
  void print() const;
  StateVectors propagated(double std_grav_param, double elapsed_time) const;

private:
  double delta_E_from_delta_M(double delta_M, double e_cos_E1,
			      double e_sin_E1) const;
};
} // namespace astro_cpp

#endif // ASTROCPP_STATE_VECTORS_HPP
