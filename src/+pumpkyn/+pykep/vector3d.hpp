#ifndef ASTROCPP_MATH_VECTOR3D_HPP
#define ASTROCPP_MATH_VECTOR3D_HPP

#include <cmath>
#include <iostream>

namespace astro_cpp {
class Vector3D {
public:
  double x, y, z;

  Vector3D(double x, double y, double z);
  Vector3D();

  double mag_sqr() const;
  double mag() const;
  void print() const;

  Vector3D &operator+();      // prefix +
  Vector3D operator-() const; // prefix -
  Vector3D &operator+=(const Vector3D &v2);
  Vector3D operator+(const Vector3D &v2) const;
  Vector3D &operator-=(const Vector3D &v2);
  Vector3D operator-(const Vector3D &v2) const;
  Vector3D &operator*=(double c);
  Vector3D &operator/=(double c);
  // operators / and * are outside the class definition

  Vector3D &invert();
  Vector3D unit_vector() const;
  static double dot_product(const Vector3D &v1, const Vector3D &v2);
  static Vector3D cross_product(const Vector3D &v1, const Vector3D &v2);
};

// Vector3D methods' implementation
inline Vector3D::Vector3D(double x, double y, double z) : x{x}, y{y}, z{z} {}

inline Vector3D::Vector3D() : x(0.0), y(0.0), z(0.0) {}

inline double Vector3D::mag() const { return std::sqrt(x * x + y * y + z * z); }

inline double Vector3D::mag_sqr() const { return x * x + y * y + z * z; }

inline void Vector3D::print() const {
  std::cout << "\nx = " << x << "\ny = " << y << "\nz = " << z
	    << "\nmagnitude = " << mag() << '\n';
}

inline Vector3D &Vector3D::operator+() // prefix +
{
  return *this;
}

inline Vector3D Vector3D::operator-() const // prefix -
{
  return Vector3D(-x, -y, -z);
}

inline Vector3D &Vector3D::operator+=(const Vector3D &v2) {
  x += v2.x;
  y += v2.y;
  z += v2.z;
  return *this;
}

inline Vector3D Vector3D::operator+(const Vector3D &v2) const {
  return Vector3D(x + v2.x, y + v2.y, z + v2.z);
}

inline Vector3D &Vector3D::operator-=(const Vector3D &v2) {
  x -= v2.x;
  y -= v2.y;
  z -= v2.z;
  return *this;
}

inline Vector3D Vector3D::operator-(const Vector3D &v2) const {
  return Vector3D(x - v2.x, y - v2.y, z - v2.z);
}

inline Vector3D &Vector3D::operator*=(double c) {
  x *= c;
  y *= c;
  z *= c;
  return *this;
}

inline Vector3D operator*(double c, Vector3D v) { return v *= c; }

inline Vector3D &Vector3D::operator/=(double c) {
  const double inv_c = 1.0 / c;
  x *= inv_c;
  y *= inv_c;
  z *= inv_c;
  return (*this);
}

inline Vector3D operator/(Vector3D v, double c) { return v /= c; }

inline Vector3D &Vector3D::invert() {
  x = -x;
  y = -y;
  z = -z;
  return *this;
}

inline Vector3D Vector3D::unit_vector() const {
  double c = 1.0 / mag();
  return c * (*this);
}

inline double Vector3D::dot_product(const Vector3D &v1, const Vector3D &v2) {
  return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

inline Vector3D Vector3D::cross_product(const Vector3D &v1,
					const Vector3D &v2) {
  return Vector3D(v1.y * v2.z - v2.y * v1.z, v1.z * v2.x - v2.z * v1.x,
		  v1.x * v2.y - v2.x * v1.y);
}
} // namespace astro_cpp

#endif // ASTROCPP_MATH_VECTOR3D_HPP
