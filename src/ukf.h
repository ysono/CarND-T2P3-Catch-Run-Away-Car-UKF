#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "tools.h"
#include "null_stream.h"
#include "Eigen/Dense"
#include <tuple>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:
  NullStream logger_;

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* time when the state is true, in us
  long time_us_;

  ///* Process noise standard deviation for longitudinal acceleration in m/s^2
  double std_process_a_;

  //* Initial state noise standard deviation for yaw rate in rad/s
  double std_process_yawd_;

  ///* Process noise standard deviation for yaw acceleration in rad/s^2
  double std_process_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laser_px_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laser_py_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radar_rho_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radar_phi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radar_rhodot_;

  MatrixXd Q_;
  MatrixXd R_laser_;
  MatrixXd R_radar_;

  MatrixXd H_laser_;
  MatrixXd H_laser_transpose_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* Weights of augmented sigma points
  ///* If we were to regenerate sigma points for measurement prediction,
  ///* then we'd need 11-element weights.
  VectorXd weights_15_;

  ///* NIS statistics
  int radar_nis_count_over_95_;
  int radar_nis_count_;

  Tools tools_;

  /**
   * Constructor
   */
  UKF(bool is_verbose);

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  std::tuple<MatrixXd, MatrixXd> Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLaser(const VectorXd &z);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(const VectorXd &z, const MatrixXd &X_sig, const MatrixXd &X_sig_diff);
};

#endif /* UKF_H */
