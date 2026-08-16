/**
 * @file moon_rotation_utilities.hpp
 * @brief Functions to calculate the moon rotation frame conversion
 * @note Ref: A Standardized Lunar Coordinate System for the Lunar Reconnaissance Orbiter and Lunar Datasets
 *            https://lunar.gsfc.nasa.gov/library/LunCoordWhitePaper-10-08.pdf
 *            https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430_moon_coord.pdf
 */

#ifndef S2E_LIBRARY_PLANET_ROTATION_MOON_MEAN_EARTH_PRINCIPAL_AXIS_FRAME_HPP_
#define S2E_LIBRARY_PLANET_ROTATION_MOON_MEAN_EARTH_PRINCIPAL_AXIS_FRAME_HPP_

#include "library/math/matrix.hpp"
#include "library/math/vector.hpp"

/**
 * @fn CalcDcmEciToPrincipalAxis
 * @brief Calculate DCM from ECI to PA (Principal Axis) moon fixed frame
 * @param[in] moon_position_eci_m: Moon position vector @ ECI frame [m]
 * @param[in] moon_velocity_eci_m_s: Moon velocity vector @ ECI frame [m/s]
 */
libra::Matrix<3, 3> CalcDcmEciToPrincipalAxis(const libra::Vector<3> moon_position_eci_m, const libra::Vector<3> moon_velocity_eci_m_s);

/**
 * @fn CalcDcmEciToMeanEarth
 * @brief Calculate DCM from ECI to ME (Mean Earth) moon fixed frame
 * @param[in] moon_position_eci_m: Moon position vector @ ECI frame [m]
 * @param[in] moon_velocity_eci_m_s: Moon velocity vector @ ECI frame [m/s]
 */
libra::Matrix<3, 3> CalcDcmEciToMeanEarth(const libra::Vector<3> moon_position_eci_m, const libra::Vector<3> moon_velocity_eci_m_s);

/**
 * @fn CalcDcmMeToPrincipalAxis
 * @brief Calculate DCM from ME (Mean Earth) moon fixed frame to PA (Principal Axis) moon fixed frame
 */
libra::Matrix<3, 3> CalcDcmMeanEarthToPrincipalAxis();

/**
 * @fn CalcDcmMciToMoonEarthSynodic
 * @brief Calculate DCM from MCI to Moon-Earth synodic frame (rotating)
 * @param[in] earth_position_mci_m: Earth position vector @ MCI frame [m]
 * @param[in] earth_velocity_mci_m_s: Earth velocity vector @ MCI frame [m/s]
 */
libra::Matrix<6, 6> CalcDcmMciToMoonEarthSynodic(const libra::Vector<3> earth_position_mci_m, const libra::Vector<3> earth_velocity_mci_m_s);

#endif  // S2E_LIBRARY_PLANET_ROTATION_MOON_MEAN_EARTH_PRINCIPAL_AXIS_FRAME_HPP_
