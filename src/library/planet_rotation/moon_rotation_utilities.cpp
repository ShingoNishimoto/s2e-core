/**
 * @file moon_rotation_utilities.cpp
 * @brief Functions to calculate the moon rotation frame conversion
 * @note Ref: A Standardized Lunar Coordinate System for the Lunar Reconnaissance Orbiter and Lunar Datasets
 *            https://lunar.gsfc.nasa.gov/library/LunCoordWhitePaper-10-08.pdf
 *            https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430_moon_coord.pdf
 */

#include "moon_rotation_utilities.hpp"

#include <library/math/constants.hpp>

libra::Matrix<3, 3> CalcDcmEciToPrincipalAxis(const libra::Vector<3> moon_position_eci_m, const libra::Vector<3> moon_velocity_eci_m_s) {
  libra::Matrix<3, 3> dcm_eci2me = CalcDcmEciToMeanEarth(moon_position_eci_m, moon_velocity_eci_m_s);
  libra::Matrix<3, 3> dcm_me2pa = CalcDcmMeanEarthToPrincipalAxis();

  return dcm_me2pa * dcm_eci2me;
}

libra::Matrix<3, 3> CalcDcmEciToMeanEarth(const libra::Vector<3> moon_position_eci_m, const libra::Vector<3> moon_velocity_eci_m_s) {
  libra::Vector<3> me_ex_eci = -1.0 * moon_position_eci_m.CalcNormalizedVector();

  libra::Vector<3> moon_orbit_norm = libra::OuterProduct(moon_position_eci_m, moon_velocity_eci_m_s);
  libra::Vector<3> me_ez_eci = moon_orbit_norm.CalcNormalizedVector();

  libra::Vector<3> me_ey_eci = libra::OuterProduct(me_ez_eci, me_ex_eci);

  libra::Matrix<3, 3> dcm_eci_to_me;
  for (size_t i = 0; i < 3; i++) {
    dcm_eci_to_me[0][i] = me_ex_eci[i];
    dcm_eci_to_me[1][i] = me_ey_eci[i];
    dcm_eci_to_me[2][i] = me_ez_eci[i];
  }

  return dcm_eci_to_me;
}

libra::Matrix<3, 3> CalcDcmMeanEarthToPrincipalAxis() {
  // The correction values between DE430 Principal Axis and Mean Earth frame
  const double theta_x_rad = 0.285 * libra::arcsec_to_rad;
  const double theta_y_rad = 78.580 * libra::arcsec_to_rad;
  const double theta_z_rad = 67.573 * libra::arcsec_to_rad;

  libra::Matrix<3, 3> dcm_me_pa =
      libra::MakeRotationMatrixZ(theta_z_rad) * libra::MakeRotationMatrixY(theta_y_rad) * libra::MakeRotationMatrixX(theta_x_rad);

  return dcm_me_pa;
}

libra::Matrix<6, 6> CalcDcmMciToMoonEarthSynodic(const libra::Vector<3> earth_position_mci_m, const libra::Vector<3> earth_velocity_mci_m_s) {
  libra::Vector<3> x_hat = earth_position_mci_m.CalcNormalizedVector();
  libra::Vector<3> z_hat = (libra::OuterProduct(earth_position_mci_m, earth_velocity_mci_m_s)).CalcNormalizedVector();
  libra::Vector<3> y_hat = libra::OuterProduct(z_hat, x_hat);

  libra::Matrix<3, 3> C;
  for (uint8_t i = 0; i < 3; i++) {
    C(0, i) = x_hat(i);
    C(1, i) = y_hat(i);
    C(2, i) = z_hat(i);
  }

  // Angular velocity
  libra::Vector<3> omega = 1.0 / libra::InnerProduct(earth_position_mci_m, earth_position_mci_m) * libra::OuterProduct(earth_position_mci_m, earth_velocity_mci_m_s);

  // Skew matrix
  libra::Matrix<3, 3> OmegaX;
  OmegaX(0, 0) = 0;         OmegaX(0, 1) = -omega(2); OmegaX(0, 2) = omega(1);
  OmegaX(1, 0) = omega(2);  OmegaX(1, 1) = 0;         OmegaX(1, 2) = -omega(0);
  OmegaX(2, 0) = -omega(1); OmegaX(2, 1) = omega(0);  OmegaX(2, 2) = 0;

  // 6x6 state transition matrix
  libra::Matrix<6, 6> T;
  T.FillUp(0);
  libra::Matrix<3, 3> C_OmegaX = C * OmegaX;
  for (uint8_t i = 0; i < 3; i++) {
    for (uint8_t j = 0; j < 3; j++) {
      T(i, j) = C(i, j);
      T(3 + i, 3 + j) = C(i, j);
      T(3 + i, j) = -C_OmegaX(i, j);
    }
  }

  return T;
}
