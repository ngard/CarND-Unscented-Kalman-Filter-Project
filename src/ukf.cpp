#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  n_x_ = 5;
  n_q_ = 2;
  n_aug_ = n_x_ + n_q_;
  lambda_ = 3 - n_aug_;
  n_sigma_ = 2*n_aug_ + 1;

  // initial state vector
  x_ = VectorXd(n_x_);
  x_.fill(0);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  P_.fill(0);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  //std_a_ = 30;
  std_a_ = 0.6;
  std_a2_ = std_a_ * std_a_;

  // Process noise standard deviation yaw acceleration in rad/s^2
  //std_yawdd_ = 30;
  std_yawdd_ = 0.6;
  std_yawdd2_ = std_yawdd_ * std_yawdd_;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;
  std_laspx2_ = std_laspx_*std_laspx_;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;
  std_laspy2_ = std_laspy_*std_laspy_;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;
  std_radr2_ = std_radr_*std_radr_;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;
  std_radphi2_ = std_radphi_ * std_radphi_;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  std_radrd2_ = std_radrd_ * std_radrd_;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  Xsig_pred_ = MatrixXd(n_x_, n_sigma_);
  Xsig_pred_.fill(0);
  weights_ = VectorXd(n_sigma_);
  for (int ii=0; ii<n_sigma_; ++ii)
    //set weights
    weights_(ii) = (ii==0) ? lambda_/(lambda_+n_aug_) : 0.5/(lambda_+n_aug_);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  //cerr << __func__ << endl;

  if (!is_initialized_) {

    if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      float rho = meas_package.raw_measurements_[0];
      float theta = meas_package.raw_measurements_[1];
      float rho_dot = meas_package.raw_measurements_[2];
      float px = rho*cos(theta), py = rho*sin(theta);
      float vel_abs = rho_dot;
      float yaw_angle = theta;
      float yaw_rate = 0;
      x_ << px, py, vel_abs, yaw_angle, yaw_rate;
      P_ <<
	1, 0, 0, 0, 0,
	0, 1, 0, 0, 0,
	0, 0, 1, 0, 0,
	0, 0, 0, 2, 0,
	0, 0, 0, 0, 5;
    }

    if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) {
      float px = meas_package.raw_measurements_[0];
      float py = meas_package.raw_measurements_[1];
      float vel_abs = 0;
      float yaw_angle = atan2(py,px);
      float yaw_rate = 0;
      x_ << px, py, vel_abs, yaw_angle, yaw_rate;
      P_ <<
	0.5, 0, 0, 0, 0,
	0, 0.5, 0, 0, 0,
	0, 0, 3, 0, 0,
	0, 0, 0, 1, 0,
	0, 0, 0, 0, 2;
    }
    
    is_initialized_ = true;
    time_us_ = meas_package.timestamp_;

    return;
  }

  float dt = (meas_package.timestamp_-time_us_)/1.e6;
  
  if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {

    Prediction(dt);

    UpdateRadar(meas_package);

    time_us_ = meas_package.timestamp_;
  }

  if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) {

    Prediction(dt);

    UpdateLidar(meas_package);

    time_us_ = meas_package.timestamp_;
  }

  cerr << "x_:" << endl;
  cerr << x_ << endl;
  cerr << "P_:" << endl;
  cerr << P_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double dt) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  //cerr << __func__ << endl;

  MatrixXd Xsig_aug(n_aug_,n_sigma_);
  GenerateAugmentedSigmaPoints(Xsig_aug);

  PredictionOnSigmaPoints(Xsig_aug, dt);

  CalcurateMeanAndCovariance();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  //cerr << __func__ << endl;

  int n_z = 2;
  MatrixXd S(n_z, n_z);
  S.fill(0);
  MatrixXd Zsig(n_z, n_sigma_);
  Zsig.fill(0);
  VectorXd z_pred(n_z);
  z_pred.fill(0);
  PredictLidarMeasurement(n_z, S, Zsig, z_pred);

  // UKF Update
  VectorXd z(n_z);
  z = meas_package.raw_measurements_;
  //calculate cross correlation matrix
  MatrixXd Tc(n_x_, meas_package.raw_measurements_.size());
  Tc.fill(0);
  for (int ii=0; ii < n_sigma_; ++ii) {
    //residual
    VectorXd z_diff = Zsig.col(ii) - z_pred;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(ii) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc += weights_(ii) * x_diff * z_diff.transpose();
  }
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //update state mean and covariance matrix
  x_ += K*z_diff;
  P_ -= K*S*K.transpose();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  //cerr << __func__ << endl;

  int n_z = 3;
  MatrixXd S(n_z, n_z);
  S.fill(0);
  MatrixXd Zsig(n_z, n_sigma_);
  Zsig.fill(0);
  VectorXd z_pred(n_z);
  z_pred.fill(0);
  PredictRadarMeasurement(n_z, S, Zsig, z_pred);

  // UKF Update
  VectorXd z(n_z);
  z = meas_package.raw_measurements_;
  //calculate cross correlation matrix
  MatrixXd Tc(n_x_, meas_package.raw_measurements_.size());
  Tc.fill(0);
  for (int ii=0; ii < n_sigma_; ++ii) {
    //residual
    VectorXd z_diff = Zsig.col(ii) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(ii) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc += weights_(ii) * x_diff * z_diff.transpose();
  }
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ += K*z_diff;
  P_ -= K*S*K.transpose();
}

void UKF::GenerateAugmentedSigmaPoints(MatrixXd& Xsig_aug) const
{
  //cerr << __func__ << endl;

  //create augmented mean state
  VectorXd x_aug(n_aug_);
  x_aug << x_, 0, 0;

  //create augmented sigma points
  MatrixXd P_aug(n_aug_,n_aug_);
  P_aug.fill(0);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  MatrixXd Q(n_q_,n_q_);
  Q << std_a2_, 0,
       0, std_yawdd2_;
  P_aug.block(n_x_,n_x_,n_q_,n_q_) = Q;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  Xsig_aug.fill(0);
  //create augmented sigma points
  for (int ii=0; ii< n_sigma_; ++ii)
    Xsig_aug.col(ii) = x_aug;
  double k = sqrt(lambda_+n_aug_);
  L *= k;

  Xsig_aug.block(0,1,n_aug_,n_aug_) += L;
  Xsig_aug.block(0,n_aug_+1,n_aug_,n_aug_) -= L;
}

void UKF::PredictionOnSigmaPoints(const MatrixXd& Xsig_aug, double dt)
{
  //cerr << __func__ << endl;

  double dt2 = dt*dt;
  for (int ii=0; ii<n_sigma_; ++ii) {
    VectorXd x = Xsig_aug.col(ii);
    double px = x(0), py = x(1), v = x(2), yaw = x(3), yaw_dot = x(4);
    double nu_a = x(5), nu_yaw_dot = x(6);

    VectorXd x_pred(n_x_);
    double px_pred = px, py_pred = py;
    double v_pred = v, yaw_pred = yaw, yaw_dot_pred = yaw_dot;
    if (fabs(yaw_dot) > 1.e-3) {
      //predict sigma points
      px_pred += v/yaw_dot*(sin(yaw+yaw_dot*dt)-sin(yaw)) + 0.5*dt2*cos(yaw)*nu_a;
      py_pred += v/yaw_dot*(-cos(yaw+yaw_dot*dt)+cos(yaw)) + 0.5*dt2*sin(yaw)*nu_a;
    } else { // yaw_dot==0
      //avoid division by zero
      px_pred += v*cos(yaw)*dt + 0.5*dt2*cos(yaw)*nu_a;
      py_pred += v*sin(yaw)*dt + 0.5*dt2*sin(yaw)*nu_a;
    }
    v_pred += dt*nu_a;
    yaw_pred += yaw_dot*dt + 0.5*dt2*nu_yaw_dot;
    yaw_dot_pred += dt*nu_yaw_dot;

    //write predicted sigma points into right column
    Xsig_pred_.col(ii) << px_pred, py_pred, v_pred, yaw_pred, yaw_dot_pred;
  }
}

void UKF::CalcurateMeanAndCovariance()
{
  //cerr << __func__ << endl;

  x_.fill(0.);
  for (int ii=0; ii<n_sigma_; ++ii) {
    //predict state mean
    x_ += weights_(ii)*Xsig_pred_.col(ii);
  }

  //predict state covariance matrix
  P_.fill(0.);
  for (int ii=0; ii<n_sigma_; ++ii) {
    // state difference
    VectorXd x_diff = Xsig_pred_.col(ii) - x_;

    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ += weights_(ii)*x_diff*x_diff.transpose();
  }
}

void UKF::PredictLidarMeasurement(int n_z, MatrixXd &S, MatrixXd &Zsig, VectorXd &z_pred)
{
  //cerr << __func__ << endl;

  //transform sigma points into measurement space
  for (int ii=0; ii < n_sigma_; ++ii) {
    VectorXd x = Xsig_pred_.col(ii);
    double px = x(0), py = x(1);

    VectorXd z(n_z);
    z << px, py;
    Zsig.col(ii) = z;

    //calculate mean predicted measurement
    z_pred += weights_(ii)*z;
  }

  //calculate innovation covariance matrix S
  for (int ii=0; ii < n_sigma_; ++ii) {
    VectorXd z_diff = Zsig.col(ii) - z_pred;
    S += weights_(ii)*z_diff*z_diff.transpose();
  }
  MatrixXd R(n_z,n_z);
  R << std_laspx2_, 0,
       0, std_laspy2_;
  S += R;
}

void UKF::PredictRadarMeasurement(int n_z, MatrixXd& S, MatrixXd& Zsig, VectorXd& z_pred)
{
  //cerr << __func__ << endl;

  //transform sigma points into measurement space
  for (int ii=0; ii < n_sigma_; ++ii) {
    VectorXd x = Xsig_pred_.col(ii);
    double px = x(0), py = x(1), v = x(2), yaw = x(3);

    double rho = sqrt(px*px+py*py);
    double phi = atan2(py,px);
    double rho_dot = (px*cos(yaw)+py*sin(yaw))*v/rho;
    VectorXd z(n_z);
    z << rho, phi, rho_dot;
    Zsig.col(ii) = z;

    //calculate mean predicted measurement
    z_pred += weights_(ii)*z;
  }

  //calculate innovation covariance matrix S
  for (int ii=0; ii < n_sigma_; ++ii) {
    VectorXd z_diff = Zsig.col(ii) - z_pred;
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    S += weights_(ii)*z_diff*z_diff.transpose();
  }
  MatrixXd R(n_z,n_z);
  R << std_radr2_, 0, 0,
       0, std_radphi2_, 0,
       0, 0, std_radrd2_;
  S += R;
}
