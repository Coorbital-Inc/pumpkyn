#ifndef ASTROCPP_MATH_CONST_HPP
#define ASTROCPP_MATH_CONST_HPP
#include <cmath>

//namespace astro_cpp::math {
namespace astro_cpp {
namespace math{
constexpr double PI = 3.14159265358979;
constexpr double INV_PI = 1 / PI;
constexpr double TAU = 2 * PI;
constexpr double INV_TAU = 1 / TAU;
constexpr double PI_OVER_2 = PI / 2;
constexpr double PI_SQR = PI * PI;
constexpr double LN_2 = 0.6931471805599453;
constexpr double RAD_TO_DEG = 180 / PI;
constexpr double DEG_TO_RAD = PI / 180;

// Stolen from StackOverflow
template <typename T> inline int sgn(const T &val) {
  return (val > T(0)) - (val < T(0));
}

// Wrap angle to [-Pi; +Pi]
template <typename T> inline double wrap_angle(const T &x) {
  // [-Tau; +Tau]
  x -= static_cast<int>(x * INV_TAU) * INV_TAU;
  // Check if already wrapped
  const bool is_wrapped = (x > -PI) || (x < PI);
  // Avoid branching
  return x - sgn(x) * is_wrapped * TAU;
}

} // namespace math
} // namespace astro_cpp
#endif // ASTROCPP_MATH_CONST_HPP
