#include "Eigen/Dense"
#define PI   3.141592653589793238462643383
using Eigen::Vector4d;
using Eigen::Matrix4d;

void initialize_piecewise_coeffs(double coeffs[4], double depth, double height,
                                 double Rp, double Xi, double well_shift) {
  // Solve following system of equations for a, b, c, d where Rp is particle
  // diameter and Rm is location of minimum, or Rp + well_shift:
  //
  //  -----------\              -
  //              \             | Height
  //               \   -------  -
  //                \ /         | Depth
  //                 -          -
  //
  //   a + b*Rp +   c*Rp^2 +   d*Rp^3 = height / 2
  // 0*a + b    + 2*c*Rp   + 3*d*Rp^2 = slope
  //   a + b*Rm +   c*Rm^2 +   d*Rm^3 = -depth
  // 0*a + b    + 2*c*Rm   + 3*d*Rm^2 = 0
  double Rp2 = Rp*Rp;
  double Rp3 = Rp*Rp2;
  double Rm = Rp + well_shift;
  double Rm2 = Rm*Rm;
  double Rm3 = Rm*Rm2;
  double slope = -height / std::sqrt(PI) / Xi;
  Matrix4d A;
  A << 1, Rp, Rp2,  Rp3,
       0, 1,  2*Rp, 3*Rp2,
       1, Rm, Rm2,  Rm3,
       0, 1,  2*Rm, 3*Rm2;
  Vector4d b;
  b << height / 2.0, slope, -depth, 0.0;
  Vector4d coeffs_v = A.fullPivLu().solve(b);
  for (int i = 0; i < 4; i++) {
    coeffs[i] = coeffs_v[i];
  }
}
