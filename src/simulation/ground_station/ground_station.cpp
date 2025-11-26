/**
 * @file ground_station.cpp
 * @brief Base class of ground station
 */

#include "ground_station.hpp"

#include <environment/global/physical_constants.hpp>
#include <library/initialize/initialize_file_access.hpp>
#include <library/logger/log_utility.hpp>
#include <library/logger/logger.hpp>
#include <library/math/constants.hpp>
#include <library/utilities/macros.hpp>
#include <string>

GroundStation::GroundStation(const SimulationConfiguration* configuration, const unsigned int ground_station_id)
    : ground_station_id_(ground_station_id) {
  Initialize(configuration, ground_station_id_);
  number_of_spacecraft_ = configuration->number_of_simulated_spacecraft_;
  for (unsigned int i = 0; i < number_of_spacecraft_; i++) {
    is_visible_[i] = false;
    range_m_[i] = 0.0;
    range_rate_m_s_[i] = 0.0;
  }
}

GroundStation::~GroundStation() {}

void GroundStation::Initialize(const SimulationConfiguration* configuration, const unsigned int ground_station_id) {
  std::string gs_ini_path = configuration->ground_station_file_list_[0];
  auto conf = IniAccess(gs_ini_path);

  const char* section_base = "GROUND_STATION_";
  const std::string section_tmp = section_base + std::to_string(static_cast<long long>(ground_station_id));
  const char* Section = section_tmp.data();

  double latitude_deg = conf.ReadDouble(Section, "latitude_deg");
  double longitude_deg = conf.ReadDouble(Section, "longitude_deg");
  double height_m = conf.ReadDouble(Section, "height_m");
  geodetic_position_ = GeodeticPosition(latitude_deg * libra::deg_to_rad, longitude_deg * libra::deg_to_rad, height_m);
  position_ecef_m_ = geodetic_position_.CalcEcefPosition();

  elevation_limit_angle_deg_ = conf.ReadDouble(Section, "elevation_limit_angle_deg");

  configuration->main_logger_->CopyFileToLogDirectory(gs_ini_path);
}

void GroundStation::LogSetup(Logger& logger) { logger.AddLogList(this); }

void GroundStation::Update(const EarthRotation& celestial_rotation, const Spacecraft& spacecraft) {
  libra::Matrix<3, 3> dcm_ecef2eci = celestial_rotation.GetDcmJ2000ToEcef().Transpose();
  position_i_m_ = dcm_ecef2eci * position_ecef_m_;

  is_visible_[spacecraft.GetSpacecraftId()] = CalcIsVisible(spacecraft);
  CalcRaRR(spacecraft);
}


bool GroundStation::CalcIsVisible(const Spacecraft& spacecraft) {
  libra::Quaternion q_ecef_to_ltc = geodetic_position_.GetQuaternionXcxfToLtc();

  libra::Vector<3> sc_pos_ltc = q_ecef_to_ltc.FrameConversion(spacecraft.GetDynamics().GetOrbit().GetPosition_ecef_m() - position_ecef_m_);  // Satellite position in LTC frame [m]
  sc_pos_ltc = sc_pos_ltc.CalcNormalizedVector();
  libra::Vector<3> dir_gs_to_zenith = libra::Vector<3>(0);
  dir_gs_to_zenith[2] = 1;

  // Judge the satellite position angle is over the minimum elevation

  if (dot(sc_pos_ltc, dir_gs_to_zenith) < sin(elevation_limit_angle_deg_ * libra::deg_to_rad)) {
    return false;
  }

  std::string center_body_name = spacecraft.GetLocalEnvironment().GetCelestialInformation().GetGlobalInformation().GetCenterBodyName();
  if (center_body_name == "EARTH")
    return true;

  // Check the occultation by center of body
  libra::Vector<3> pos_center_body_eci = - spacecraft.GetLocalEnvironment().GetCelestialInformation().GetGlobalInformation().GetPositionFromCenter_i_m("EARTH");
  libra::Vector<3> los_gs_to_center_body_i = pos_center_body_eci - position_i_m_;
  double distance_gs_to_center_body = los_gs_to_center_body_i.CalcNorm();
  double zenith_edge = asin(spacecraft.GetLocalEnvironment().GetCelestialInformation().GetGlobalInformation().GetMeanRadiusFromName_m(center_body_name.c_str()) / distance_gs_to_center_body);
  double distance_edge = distance_gs_to_center_body * cos(zenith_edge);

  libra::Vector<3> los_gs_to_sc_i = (spacecraft.GetDynamics().GetOrbit().GetPosition_i_m() + pos_center_body_eci) - position_i_m_;
  double distance_sc = los_gs_to_sc_i.CalcNorm();
  los_gs_to_center_body_i = los_gs_to_center_body_i.CalcNormalizedVector();
  los_gs_to_sc_i = los_gs_to_sc_i.CalcNormalizedVector();
  double dot_sc_center_body = libra::InnerProduct(los_gs_to_center_body_i, los_gs_to_sc_i);
  // Clip
  if (fabs(dot_sc_center_body) > 1)
    dot_sc_center_body /= fabs(dot_sc_center_body);
  double zenith_sc = acos(dot_sc_center_body);
  if (distance_sc > distance_edge && zenith_sc < zenith_edge)
    return false;

  return true;
}


void GroundStation::CalcRaRR(const Spacecraft& spacecraft) {
  libra::Vector<3> rel_sc_position = spacecraft.GetDynamics().GetOrbit().GetPosition_ecef_m() - position_ecef_m_;
  const unsigned int sc_id = spacecraft.GetSpacecraftId();
  // Instantaneous geometric distance
  range_m_.at(sc_id) = (rel_sc_position).CalcNorm();
  // Instantaneous geometric range-rate NOTE: no need to consider inertial correction, since it cancel out.
  range_rate_m_s_.at(sc_id) = libra::InnerProduct(spacecraft.GetDynamics().GetOrbit().GetVelocity_ecef_m_s(), rel_sc_position.CalcNormalizedVector());
}


std::string GroundStation::GetLogHeader() const {
  std::string str_tmp = "";

  std::string head = "ground_station" + std::to_string(ground_station_id_) + "_";
  for (unsigned int i = 0; i < number_of_spacecraft_; i++) {
    std::string legend_base = head + "sc" + std::to_string(i);
    str_tmp += WriteScalar(legend_base + "_visible_flag");
    str_tmp += WriteScalar(legend_base + "_range", "m");
    str_tmp += WriteScalar(legend_base + "_range_rate", "m/s");
  }
  str_tmp += WriteVector("ground_station_position", "eci", "m", 3);
  // str_tmp += WriteVector("ground_station_position", "ecef", "m", 3);
  return str_tmp;
}

std::string GroundStation::GetLogValue() const {
  std::string str_tmp = "";

  for (unsigned int i = 0; i < number_of_spacecraft_; i++) {
    str_tmp += WriteScalar(is_visible_.at(i));
    str_tmp += WriteScalar(range_m_.at(i), 16);
    str_tmp += WriteScalar(range_rate_m_s_.at(i), 16);
  }
  str_tmp += WriteVector(position_i_m_, 16);
  // str_tmp += WriteVector(position_ecef_m_);
  return str_tmp;
}
