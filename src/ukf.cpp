#include "ukf.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::endl;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF(bool is_verbose) {
  if (is_verbose) {
    logger_.rdbuf(std::cout.rdbuf());
  }
  
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // state vector
  x_ = VectorXd(5);

  // process covariance matrix
  P_ = MatrixXd(5, 5);

  /*
   * Absolute value of acceleration should be the greatest when braking.
   * Assume our simulated bicycle can brake from (10kph == 2.8m/s) to stop in 1 sec
   *   => max absolute value of acceleration is 2.8 m/s^2
   * Use the rule of thumb
   *   => (2.8m/s^2 / 2 == 1.4m/s^2) is the std of acceleration.
   *
   * Assume our bicycle can turn at an equivalent of rotating 2pi degrees in 6 sec
   *   => max absolute value of yaw rate is pi/3 rad/s
   * Use the rule of thumb
   *   => (pi/3 rad/s / 2 == pi/6 rad/s) is the std of yaw rate.
   *
   * Assume our bicycle can switch (or wobble) from its tightest possible curve
   * to the other tightest possible curve in the opposite direction, over 2sec
   *   => max absolute value of yaw acceleration is (pi/3 rad/s * 2 / 2s == pi/3 rad/s^2)
   * Use the rule of thumb
   *   => (pi/3 rad/s^2 / 2 == pi/6 rad/s^2) is the std of yaw acceleration.
   */

  // Process noise standard deviation for longitudinal acceleration in m/s^2
  std_process_a_ = 1.4;

  // Initial state noise standard deviation for yaw rate in rad/s
  std_process_yawd_  = M_PI / 6;

  // Process noise standard deviation for yaw acceleration in rad/s^2
  std_process_yawdd_ = M_PI / 6;

  Q_ = MatrixXd(2, 2);
  Q_ << pow(std_process_a_, 2), 0,
        0, pow(std_process_yawdd_, 2);
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laser_px_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laser_py_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radar_rho_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radar_phi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radar_rhodot_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << pow(std_laser_px_, 2), 0,
              0, pow(std_laser_py_, 2);

  R_radar_ = MatrixXd(3, 3);
  R_radar_.fill(0.0);
  R_radar_(0, 0) = pow(std_radar_rho_, 2);
  R_radar_(1, 1) = pow(std_radar_phi_, 2);
  R_radar_(2, 2) = pow(std_radar_rhodot_, 2);

  H_laser_ = MatrixXd(2, 5);
  H_laser_.fill(0.0);
  H_laser_(0, 0) = 1.0;
  H_laser_(1, 1) = 1.0;

  H_laser_transpose_ = H_laser_.transpose();

  int n_aug = x_.size() + Q_.rows();
  lambda_ = 3 - n_aug;
  weights_15_ = VectorXd(1 + 2 * n_aug);
  weights_15_(0) = lambda_ / (lambda_ + n_aug);
  double weight_rest = 1 / (lambda_ + n_aug) / 2;
  for (int i = 1; i < weights_15_.size(); i++) {
    weights_15_(i) = weight_rest;
  }

  logger_ << "lambda is " << lambda_ << endl
    << "weights_15_" << endl << weights_15_ << endl;
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  if (is_initialized_) {
    double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
    MatrixXd X_pred, X_pred_diff;
    std::tie(X_pred, X_pred_diff) = Prediction(delta_t);

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
      // Use X_pred as the sigma points, without generating new correct sigma points based on x and P.
      UpdateRadar(meas_package.raw_measurements_, X_pred, X_pred_diff);
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
      UpdateLaser(meas_package.raw_measurements_);
    }
  } else {
    // Fill P_ with some default values, skipping those values that can be filled
    // by the initial measurement.
    P_.fill(0.0);
    
    // Max uncertainty about yaw angle puts it within a 2pi range => std of pi is more than large enough.
    P_(3, 3) = pow(M_PI, 2);
    
    P_(4, 4) = pow(std_process_yawd_, 2);

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      logger_ << "initial meas is radar" << endl
        << "default P is " << endl << P_ << endl;

      x_ = tools_.radar_z_to_x(meas_package.raw_measurements_);
      tools_.initial_radar_P(P_, R_radar_, Q_);
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      // A bicycle in the city should max 15mph == 7m/s. Use the rule of thumb.
      // Only lidar, not radar, needs a default variance of v.
      P_(2, 2) = pow(7 / 2, 2);

      logger_ << "initial meas is lidar" << endl
        << "default P is " << endl << P_ << endl;

      x_ = tools_.laser_z_to_x(meas_package.raw_measurements_);
      tools_.initial_laser_P(P_, R_laser_);
    }

    logger_ << "initial P is " << endl << P_ << endl
      << "first z is " << endl << meas_package.raw_measurements_ << endl
      << "initial x is " << endl << x_ << endl;

    is_initialized_ = true;
  }

  time_us_ = meas_package.timestamp_;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 * @return 1) The predicted X sigma points.
 *         2) The diff between each predicted sigma point and all predicted sigma points' weighted sum.
 *            This will be used for calculating T matrix during measurement update.
 *            If we wanted to regenerate sigma points for measurement update, we wouldn't need to return this here.
 */
std::tuple<MatrixXd, MatrixXd> UKF::Prediction(double delta_t) {
  MatrixXd X_sig_aug = tools_.augmented_sigma_points(x_, P_, Q_, lambda_);

  logger_ << "X_sig_aug" << endl << X_sig_aug << endl;

  MatrixXd X_pred = tools_.predict_x(X_sig_aug, delta_t);

  logger_ << "X_pred" << endl << X_pred << endl;

  x_.fill(0.0);
  for (int i = 0; i < X_pred.cols(); i++) {
    x_ += weights_15_(i) * X_pred.col(i);
  }
  // Deriving P below will need to diff each sigma x against the aggregated x_pred.
  // Hence the yaw angle does NOT need to be normalized yet.

  logger_ << "x_pred" << endl << x_ << endl;

  MatrixXd X_pred_diff(X_pred.rows(), X_pred.cols());

  //predicted state covariance matrix
  P_.fill(0.0);

  for (int i = 0; i < X_pred.cols(); i++) {
    VectorXd x_diff = X_pred.col(i) - x_;

    // Need to normalize the diff psi, regardless of whether the psi values in
    // both X_pred and x_pred had been normalized.
    // In practice, yaw rate does not need to be normalized.
    tools_.normalize_angle(x_diff(3));

    X_pred_diff.col(i) = x_diff;

    P_ += weights_15_(i) * x_diff * x_diff.transpose();
  }

  logger_ << "diff between X_pred and x_pred" << endl << X_pred_diff << endl
    << "P_pred" << endl << P_ << endl;

  // Normalize the yaw angle after P has been calculated.
  // In practice, yaw rate does not need to be normalized.
  tools_.normalize_angle(x_(3));

  return std::make_tuple(X_pred, X_pred_diff);
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLaser(const VectorXd &z) {
  // standard KF.

  MatrixXd S = H_laser_ * P_ * H_laser_transpose_ + R_laser_;
  MatrixXd K = P_ * H_laser_transpose_ * S.inverse();

  VectorXd z_pred = H_laser_ * x_;
  x_ += K * (z - z_pred);

  P_ -= K * H_laser_ * P_;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(const VectorXd &z, const MatrixXd &X_sig, const MatrixXd &X_sig_diff) {

  MatrixXd Z_pred = tools_.x_to_radar_z(X_sig);

  logger_ << "radar Z_pred" << endl << Z_pred << endl;

  VectorXd z_pred(3);
  z_pred.fill(0.0);
  for (int i = 0; i < Z_pred.cols(); i++) {
    z_pred += weights_15_(i) * Z_pred.col(i);
  }
  // The phis in Z_pred are normalized b/c of `atan2()`.
  // Hence the phi in z_pred is also normalized.

  logger_ << "radar z_pred" << endl << z_pred << endl;

  MatrixXd S(3, 3);
  S.topLeftCorner(3, 3) = R_radar_; // What's Eigen's proper way of initializing all elements as a copy of another matrix?

  MatrixXd T(5, 3);
  T.fill(0.0);

  for (int i = 0; i < Z_pred.cols(); i++) {
    VectorXd z_diff = Z_pred.col(i) - z_pred;

    // Need to normalize the diff phi, even though the phi values in
    // both Z_pred and z_pred are normalized.
    tools_.normalize_angle(z_diff(1));

    MatrixXd z_diff_transpose = z_diff.transpose();

    S += weights_15_(i) * z_diff * z_diff_transpose;

    // The yaw angle in x_diff is already normalized.
    T += weights_15_(i) * X_sig_diff.col(i) * z_diff_transpose;
  }

  MatrixXd S_inverse = S.inverse();
  MatrixXd K = T * S_inverse;

  VectorXd y = z - z_pred;
  x_ += K * y;
  tools_.normalize_angle(x_(3));

  P_ -= K * S * K.transpose();

  logger_ << "radar z" << endl << z << endl
    << "updated x" << endl << x_ << endl
    << "updated P" << endl << P_ << endl;

  double nis = y.transpose() * S_inverse * y;
  if (nis > 7.815) radar_nis_count_over_95_++;
  radar_nis_count_++;
}
