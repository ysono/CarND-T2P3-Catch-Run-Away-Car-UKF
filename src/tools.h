#ifndef TOOLS_H_
#define TOOLS_H_
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class Tools {
public:

  void normalize_angle(double &angle);

  VectorXd radar_z_to_x(const VectorXd &z);
  void initial_radar_P(MatrixXd &P, const MatrixXd &R, const MatrixXd &Q);

  VectorXd laser_z_to_x(const VectorXd &z);
  void initial_laser_P(MatrixXd &P, const MatrixXd &R);

  MatrixXd augmented_sigma_points(const VectorXd &x, const MatrixXd &P, const MatrixXd &Q, const double lambda);
  MatrixXd sigma_points(const VectorXd &x, const MatrixXd &P, const double lambda);

  MatrixXd predict_x(const MatrixXd &X_sig_aug, const double delta_t);

  MatrixXd x_to_radar_z(const MatrixXd &X_sig);

  VectorXd CalculateRMSE(VectorXd &sse,
                         int &sse_count,
                         const VectorXd &estimation,
                         const VectorXd &ground_truth);

};

#endif /* TOOLS_H_ */