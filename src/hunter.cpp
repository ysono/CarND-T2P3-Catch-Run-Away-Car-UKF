#include "Eigen/Dense"
#include "null_stream.h"
#include <iostream>

using Eigen::VectorXd;
using std::endl;

class Hunter {
private:

  NullStream logger_;

  VectorXd prev_target_;

public:

  Hunter(bool is_verbose) {
    if (is_verbose) {
      logger_.rdbuf(std::cout.rdbuf());
    }

    prev_target_ = VectorXd(2);
    prev_target_.fill(0.0);
  }

  VectorXd calc_target(const VectorXd &ukf_state, const double &hunter_x, const double hunter_y) {
    double px = ukf_state(0);
    double py = ukf_state(1);
    double v = ukf_state(2);
    double yaw = ukf_state(3);
    double yawdot = ukf_state(4);
    logger_ << "position is " << px << " " << py << endl;
    logger_ << "yaw and yawdot are " << yaw << " " << yawdot << endl;

    if (yawdot < 1e-3) {
      logger_ << "Straight line detected. Can't catch, so no-op this round." << endl;
      return prev_target_; // initially [0, 0].
    }

    double radius = v / yawdot;
    logger_ << "radius " << radius << endl;

    double center_x = px - radius * sin(yaw);
    double center_y = py + radius * cos(yaw);
    logger_ << "center is " << center_x << " " << center_y << endl;

    double angle_from_center_to_hunter = atan2(hunter_y - center_y, hunter_x - center_x);

    double distance_between_hunter_and_center = sqrt(pow(center_y - hunter_y, 2) + pow(center_x - hunter_x, 2));
    double distance_off_circumference_threshold = radius / 4;
    if (abs(distance_between_hunter_and_center - radius) < distance_off_circumference_threshold) {

      // more accurate to rely on directly measured px and py than on UKF-implied yaw.
      double angle_from_center_to_target = atan2(py - center_y, px - center_x);

      // yawdot as in radians per simulation event, not radians per sec.
      double hunter_yawdot = (angle_from_center_to_target - angle_from_center_to_hunter);
      // need to slow down as distance reduces, due to simulation's time resolution.
      hunter_yawdot /= 5;
      if (hunter_yawdot * yawdot > 0) {
        // if hunter would be chasing in the same direction, flip the direction.
        hunter_yawdot *= -1;
      }
      logger_ << "hunter_hawdot " << hunter_yawdot << endl;

      angle_from_center_to_hunter += hunter_yawdot;
    }

    double target_x = center_x + radius * cos(angle_from_center_to_hunter);
    double target_y = center_y + radius * sin(angle_from_center_to_hunter);

    prev_target_ << target_x, target_y;
    return prev_target_;
  }
};
