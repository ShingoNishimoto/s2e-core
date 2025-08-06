/**
 * @file orbit.cpp
 * @brief Base class of orbit propagation
 */
#include "orbit.hpp"

libra::Quaternion Orbit::CalcQuaternion_i2lvlh() const {
  libra::Vector<3> lvlh_x = spacecraft_position_i_m_;  // x-axis in LVLH frame is position vector direction from geocenter to satellite
  libra::Vector<3> lvlh_ex = lvlh_x.CalcNormalizedVector();
  libra::Vector<3> lvlh_z =
      OuterProduct(spacecraft_position_i_m_, spacecraft_velocity_i_m_s_);  // z-axis in LVLH frame is angular momentum vector direction of orbit
  libra::Vector<3> lvlh_ez = lvlh_z.CalcNormalizedVector();
  libra::Vector<3> lvlh_y = OuterProduct(lvlh_z, lvlh_x);
  libra::Vector<3> lvlh_ey = lvlh_y.CalcNormalizedVector();

  libra::Matrix<3, 3> dcm_i2lvlh;
  dcm_i2lvlh[0][0] = lvlh_ex[0];
  dcm_i2lvlh[0][1] = lvlh_ex[1];
  dcm_i2lvlh[0][2] = lvlh_ex[2];
  dcm_i2lvlh[1][0] = lvlh_ey[0];
  dcm_i2lvlh[1][1] = lvlh_ey[1];
  dcm_i2lvlh[1][2] = lvlh_ey[2];
  dcm_i2lvlh[2][0] = lvlh_ez[0];
  dcm_i2lvlh[2][1] = lvlh_ez[1];
  dcm_i2lvlh[2][2] = lvlh_ez[2];

  libra::Quaternion q_i2lvlh = libra::Quaternion::ConvertFromDcm(dcm_i2lvlh);
  return q_i2lvlh.Normalize();
}

// NOTE: Currently Earth and Moon rotation exists.
void Orbit::TransformIToXcxf(std::string center_body_name) {
  libra::Matrix<3, 3> dcm_i_to_xcxf;
  libra::Vector<3> celestial_body_angular_velocity_i_rad_s{0.0};
  if (center_body_name == "EARTH") {
    dcm_i_to_xcxf = celestial_information_->GetEarthRotation().GetDcmJ2000ToEcef();
    celestial_body_angular_velocity_i_rad_s[2] = environment::earth_mean_angular_velocity_rad_s;
  }
  else if (center_body_name == "MOON") {
    dcm_i_to_xcxf = celestial_information_->GetMoonRotation().GetDcmJ2000ToMcmf();
    celestial_body_angular_velocity_i_rad_s[2] = environment::lunar_mean_angular_velocity_rad_s;
  }
  else {
    // TODO: implement others, or force to Earth.
  }

  libra::Vector<3> position_i_m = spacecraft_position_i_m_;
  std::string orbit_center_body_name = celestial_information_->GetCenterBodyName();
  if (orbit_center_body_name != center_body_name) {
    position_i_m += celestial_information_->GetPositionFromSelectedBody_i_m(orbit_center_body_name.c_str(), center_body_name.c_str());
  }

  libra::Vector<3> position_xcxf_m = dcm_i_to_xcxf * position_i_m;

  // convert velocity vector in Inertial frame to the vector in XCXF
  libra::Vector<3> we_cross_r = OuterProduct(celestial_body_angular_velocity_i_rad_s, position_i_m);
  libra::Vector<3> velocity_w_cross_r = spacecraft_velocity_i_m_s_ - we_cross_r;

  libra::Vector<3> velocity_xcxf_m_s = dcm_i_to_xcxf * velocity_w_cross_r;

  if (center_body_name == "EARTH") {
    spacecraft_position_ecef_m_ = position_xcxf_m;
    spacecraft_velocity_ecef_m_s_ = velocity_xcxf_m_s;
  }
  else {
    spacecraft_position_xcxf_m_ = position_xcxf_m;
    spacecraft_velocity_xcxf_m_s_= velocity_xcxf_m_s;
  }
}

// FIXME: This is only for Earth, need to update for other than geodetic.
void Orbit::TransformEcefToGeodetic(void) {
  spacecraft_geodetic_position_.UpdateFromEcef(spacecraft_position_ecef_m_);
  // Check altitude
  if (spacecraft_geodetic_position_.GetAltitude_m() < 0.0) {
    std::cout << "[Error Orbit]: The spacecraft altitude is smaller than zero." << std::endl;
    std::cout << "               The orbit or disturbance setting may have something wrong." << std::endl;
    std::exit(1);
  }
}

OrbitInitializeMode SetOrbitInitializeMode(const std::string initialize_mode) {
  if (initialize_mode == "DEFAULT") {
    return OrbitInitializeMode::kDefault;
  } else if (initialize_mode == "POSITION_VELOCITY_I") {
    return OrbitInitializeMode::kInertialPositionAndVelocity;
  } else if (initialize_mode == "ORBITAL_ELEMENTS") {
    return OrbitInitializeMode::kOrbitalElements;
  } else {
    std::cerr << "WARNINGS: orbit initialize mode is not defined!" << std::endl;
    std::cerr << "The orbit is automatically initialized as default mode" << std::endl;
    return OrbitInitializeMode::kDefault;
  }
}

std::string Orbit::GetLogHeader() const {
  std::string str_tmp = "";

  str_tmp += WriteVector("spacecraft_position", "i", "m", 3);
  str_tmp += WriteVector("spacecraft_position", "ecef", "m", 3);
  str_tmp += WriteVector("spacecraft_velocity", "i", "m/s", 3);
  str_tmp += WriteVector("spacecraft_velocity", "b", "m/s", 3);
  str_tmp += WriteVector("spacecraft_velocity", "ecef", "m/s", 3);
  str_tmp += WriteVector("spacecraft_acceleration", "i", "m/s2", 3);
  str_tmp += WriteScalar("spacecraft_latitude", "rad");
  str_tmp += WriteScalar("spacecraft_longitude", "rad");
  str_tmp += WriteScalar("spacecraft_altitude", "m");

  return str_tmp;
}

std::string Orbit::GetLogValue() const {
  std::string str_tmp = "";

  str_tmp += WriteVector(spacecraft_position_i_m_, 16);
  str_tmp += WriteVector(spacecraft_position_ecef_m_, 16);
  str_tmp += WriteVector(spacecraft_velocity_i_m_s_, 10);
  str_tmp += WriteVector(spacecraft_velocity_b_m_s_, 10);
  str_tmp += WriteVector(spacecraft_velocity_ecef_m_s_, 10);
  str_tmp += WriteVector(spacecraft_acceleration_i_m_s2_, 10);
  str_tmp += WriteScalar(spacecraft_geodetic_position_.GetLatitude_rad());
  str_tmp += WriteScalar(spacecraft_geodetic_position_.GetLongitude_rad());
  str_tmp += WriteScalar(spacecraft_geodetic_position_.GetAltitude_m());

  return str_tmp;
}
