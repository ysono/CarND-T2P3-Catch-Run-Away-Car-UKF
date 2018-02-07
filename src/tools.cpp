#include <math.h>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;

void Tools::normalize_angle(double &angle) {
  while (angle >= M_PI) angle -= 2. * M_PI;
  while (angle < -M_PI) angle += 2. * M_PI;
}

VectorXd Tools::radar_z_to_x(const VectorXd &z) {
  double rho = z(0);
  double phi = z(1);
  double rho_dot = z(2);

  VectorXd x(5);
  
  // In absence of other info, conversion has to assume there is no lateral
  // movement, ie the foreign object is travelling in a straight line towards or
  // away from the observer.
  // x << rho * cos(phi),
  //      rho * sin(phi),
  //      rho_dot,
  //      phi,
  //      0;

  // Or alternatively, assume nothing about the foreign object's speed and direction of travel.
  x << rho * cos(phi),
       rho * sin(phi),
       0,
       0,
       0;

  return x;
}

void Tools::initial_radar_P(MatrixXd &P, const MatrixXd &R, const MatrixXd &Q) {
  double var_rho = R(0, 0);
  double var_rhodot = R(2, 2);

  double var_accel = Q(0, 0);

  // We want variance of `rho * cos(phi)` or its sin equivalent.
  // The max variance of rho should be >= the desired variance whose expression
  // has an additional variable of phi. Hence variance of rho should be substitutable.
  P(0, 0) = var_rho; 
  P(1, 1) = var_rho;

  // Process acceleration noise is an assumption used for any direction, whereas
  // radar's variance of rho_dot is valid for the secant direction only.
  // Choose whichever is greater.
  P(2, 2) = std::max(var_rhodot, var_accel);
}

VectorXd Tools::laser_z_to_x(const VectorXd &z) {
  double px = z(0);
  double py = z(1);

  VectorXd x(5);
  // It's not helpful to assume yaw angle as atan2(py, px); just use zero.
  x << px, py, 0, 0, 0;
  return x;
}

void Tools::initial_laser_P(MatrixXd &P, const MatrixXd &R) {
  double var_px = R(0, 0);
  double var_py = R(1, 1);

  P(0, 0) = var_px;
  P(1, 1) = var_py;
}

MatrixXd Tools::augmented_sigma_points(const VectorXd &x, const MatrixXd &P, const MatrixXd &Q, const double lambda) {
  int num_states = x.size();
  int num_aug = Q.rows();
  int num_aug_states = num_states + num_aug;

  //create augmented mean state
  VectorXd x_aug(num_aug_states);
  x_aug.head(num_states) = x;
  x_aug.tail(num_aug).fill(0.0);

  //create augmented covariance matrix
  MatrixXd P_aug(num_aug_states, num_aug_states);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(num_states, num_states) = P;
  P_aug.bottomRightCorner(num_aug, num_aug) = Q;

  return sigma_points(x_aug, P_aug, lambda);
}

// If we were regenerating sigma points for h(x), we would use this function directly.
MatrixXd Tools::sigma_points(const VectorXd &x, const MatrixXd &P, const double lambda) {
  int num_states = x.size();
  int num_sigmas = 1 + 2 * num_states;

  MatrixXd X_sig(num_states, num_sigmas);

  X_sig.col(0) = x;

  //create square root matrix
  MatrixXd L = P.llt().matrixL();
  L *= sqrt(lambda + num_states);

  for (int i = 0; i < num_states; i++) {
    X_sig.col(1 + i) =
      x + L.col(i);
    X_sig.col(1 + num_states + i) =
      x - L.col(i);
  }

  return X_sig;
}

MatrixXd Tools::predict_x(const MatrixXd &X_sig_aug, const double delta_t) {
  int num_sigmas = X_sig_aug.cols();

  MatrixXd X_pred(5, num_sigmas);

  for (int i = 0; i < num_sigmas; i++) {
    //extract values for better readability
    double p_x = X_sig_aug(0,i);
    double p_y = X_sig_aug(1,i);
    double v = X_sig_aug(2,i);
    double yaw = X_sig_aug(3,i);
    double yawd = X_sig_aug(4,i);
    double nu_a = X_sig_aug(5,i);
    double nu_yawdd = X_sig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
      px_p = p_x + v/yawd * ( sin(yaw + yawd*delta_t) - sin(yaw));
      py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
      px_p = p_x + v*delta_t*cos(yaw);
      py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    // Yaw angle must NOT be normalized here, b/c their weighted sum will be calculated.
    // Otherwise there might be values close to +pi and -pi, which would be summed to approx zero.
    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    X_pred(0,i) = px_p;
    X_pred(1,i) = py_p;
    X_pred(2,i) = v_p;
    X_pred(3,i) = yaw_p;
    X_pred(4,i) = yawd_p;
  }

  return X_pred;
}

MatrixXd Tools::x_to_radar_z(const MatrixXd &X_sig) {
  MatrixXd Z_pred(3, X_sig.cols());

  for (int i = 0; i < X_sig.cols(); i++) {
    VectorXd x_motion_pred = X_sig.col(i);
    double px = X_sig(0, i);
    double py = X_sig(1, i);
    double v = X_sig(2, i);
    double psi = X_sig(3, i);
    
    double rho = sqrt(pow(px, 2) + pow(py, 2));
    double phi = atan2(py, px);
    double rho_dot = (px * v * cos(psi) + py * v * sin(psi)) / rho;

    // Phi must NOT be normalized here, b/c their weighted sum will be calculated.
    // Otherwise there might be values close to +pi and -pi, which would be summed to approx zero.

    Z_pred(0, i) = rho;
    Z_pred(1, i) = phi;
    Z_pred(2, i) = rho_dot;
  }

  return Z_pred;
}

VectorXd Tools::CalculateRMSE(VectorXd &sse,
                              int &sse_count,
                              const VectorXd &estimation,
                              const VectorXd &ground_truth) {
  auto err = (estimation - ground_truth).array();
  VectorXd sq = err * err;
  sse += sq;
  sse_count ++;

  VectorXd mse = sse / sse_count;
  VectorXd rmse = mse.array().sqrt();
  return rmse;
}
