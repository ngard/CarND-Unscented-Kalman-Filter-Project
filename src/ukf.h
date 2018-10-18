#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
private:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

public:
  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

private:
  ///* state covariance matrix
  MatrixXd P_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* time when the state is true, in us
  long long time_us_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;
  double std_a2_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;
  double std_yawdd2_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;
  double std_laspx2_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;
  double std_laspy2_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;
  double std_radr2_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;
  double std_radphi2_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_;
  double std_radrd2_;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* Number of Sigma points
  int n_sigma_;

  ///* Q dimension
  int n_q_;

public:
  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

private:
  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);

  void GenerateAugmentedSigmaPoints(MatrixXd& Xsig_aug) const;
  void PredictionOnSigmaPoints(const MatrixXd& Xsig_aug, double dt);
  void CalcurateMeanAndCovariance();

  void PredictLidarMeasurement(int n_z, MatrixXd& S, MatrixXd& Zsig, VectorXd& z_pred);
  void PredictRadarMeasurement(int n_z, MatrixXd& S, MatrixXd& Zsig, VectorXd& z_pred);
};

#endif /* UKF_H */
