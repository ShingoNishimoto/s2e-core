/**
 * @file ground_station.hpp
 * @brief Base class of ground station
 */

#ifndef S2E_SIMULATION_GROUND_STATION_GROUND_STATION_HPP_
#define S2E_SIMULATION_GROUND_STATION_GROUND_STATION_HPP_

#include <library/geodesy/geodetic_position.hpp>
#include <library/math/vector.hpp>
#include <simulation/spacecraft/spacecraft.hpp>

#include "../simulation_configuration.hpp"

// FIXME: it assumes the gs is on Earth.
/**
 * @class GroundStation
 * @brief Base class of ground station
 */
class GroundStation : public ILoggable {
 public:
  /**
   * @fn GroundStation
   * @brief Constructor
   */
  GroundStation(const SimulationConfiguration* configuration, const unsigned int ground_station_id_);
  /**
   * @fn ~GroundStation
   * @brief Destructor
   */
  virtual ~GroundStation();

  /**
   * @fn Initialize
   * @brief Virtual function to initialize the ground station
   */
  virtual void Initialize(const SimulationConfiguration* configuration, const unsigned int ground_station_id);
  /**
   * @fn LogSetup
   * @brief Virtual function to log output setting for ground station related components
   */
  virtual void LogSetup(Logger& logger);
  /**
   * @fn Update
   * @brief Virtual function of main routine
   */
  virtual void Update(const EarthRotation& celestial_rotation, const Spacecraft& spacecraft);

  // Override functions for ILoggable
  /**
   * @fn GetLogHeader
   * @brief Override function of log header setting
   */
  virtual std::string GetLogHeader() const;
  /**
   * @fn GetLogValue
   * @brief Override function of log value setting
   */
  virtual std::string GetLogValue() const;

  // Getters
  /**
   * @fn GetGroundStationId
   * @brief Return ground station ID
   */
  int GetGroundStationId() const { return ground_station_id_; }
  /**
   * @fn GetGeodeticPosition
   * @brief Return ground station position in the geodetic frame
   */
  GeodeticPosition GetGeodeticPosition() const { return geodetic_position_; }
  /**
   * @fn GetPosition_xcxf_m
   * @brief Return ground station position in the ECEF frame [m]
   */
  Vector<3> GetPosition_xcxf_m() const { return position_ecef_m_; }
  /**
   * @fn GetPosition_i_m
   * @brief Return ground station position in the inertial frame [m]
   */
  Vector<3> GetPosition_i_m() const { return position_i_m_; }
  /**
   * @fn GetElevationLimitAngle_deg
   * @brief Return ground station elevation limit angle [deg]
   */
  double GetElevationLimitAngle_deg() const { return elevation_limit_angle_deg_; }
  /**
   * @fn IsVisible
   * @brief Return visible flag for the target spacecraft
   * @param [in] spacecraft_id: target spacecraft ID
   */
  bool IsVisible(const unsigned int spacecraft_id) const { return is_visible_.at(spacecraft_id); }

 protected:
  unsigned int ground_station_id_;      //!< Ground station ID
  GeodeticPosition geodetic_position_;  //!< Ground Station Position in the geodetic frame
  Vector<3> position_ecef_m_{0.0};      //!< Ground Station Position in the ECEF frame [m]
  Vector<3> position_i_m_{0.0};         //!< Ground Station Position in the inertial frame [m]
  double elevation_limit_angle_deg_;    //!< Minimum elevation angle to work the ground station [deg]

  std::map<int, double> range_m_;         //!< Geometric range between GS and each spacecraft ID (not care antenna)
  std::map<int, double> range_rate_m_s_;  //!< Geometric range rate between GS and each spacecraft ID (not care antenna)
  std::map<int, bool> is_visible_;        //!< Visible flag for each spacecraft ID (not care antenna)
  unsigned int number_of_spacecraft_;     //!< Number of spacecraft in the simulation

  /**
   * @fn CalcIsVisible
   * @brief Calculate the visibility for the target spacecraft
   * @param [in] spacecraft: spacecraft
   * @return True when the satellite is visible from the ground station
   */
  bool CalcIsVisible(const Spacecraft& spacecraft);
  /**
   * @fn CalcIsRaRR
   * @brief Calculate the range and range rate for the target spacecraft
   * @param [in] spacecraft: spacecraft
   */
  void CalcRaRR(const Spacecraft& spacecraft);
};

#endif  // S2E_SIMULATION_GROUND_STATION_GROUND_STATION_HPP_
