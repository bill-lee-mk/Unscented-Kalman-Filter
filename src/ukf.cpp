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

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  //P_.setIdentity(5, 5);

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // Number of sigma points
  n_sigma_points_ = 2 * n_aug_ + 1;

  //define sigma point spreading parameter
  lambda_ = 3 - n_aug_; 

  //create sigma point matrix
  MatrixXd Xsig_ = MatrixXd(n_x_, n_sigma_points_);

  //create augmented sigma point matrix
  MatrixXd Xsig_aug_ = MatrixXd(n_aug_, n_sigma_points_);

  //create predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, n_sigma_points_);

  // Weights of sigma points
  weights_ = VectorXd(n_sigma_points_);

  double weights_0 = lambda_ / (lambda_ + n_aug_);
  weights_(0) = weights_0;
  for (int i = 1; i < n_sigma_points_; i++) {  //2n+1 weights
	double weight = 0.5 / (n_aug_ + lambda_);
	weights_(i) = weight;
  }

  // Initialize Normalized Innovation Squared (NIS) value for both sensors
  NIS_laser_ = 0.;
  NIS_radar_ = 0.;

  // Sensor's measurement dimension
  n_z_radar_ = 3; // radar -> rho, phi, rho_dot
  n_z_lidar_ = 2; // lidar -> px, py

  // Measurement covariance matrices
  R_lidar_ = MatrixXd(n_z_lidar_, n_z_lidar_);
  R_radar_ = MatrixXd(n_z_radar_, n_z_radar_);


  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  /*
  
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_radr_ * std_radr_, 0, 0,
	  0, std_radphi_*std_radphi_, 0,
	  0, 0, std_radrd_*std_radrd_;
  
  */
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar measurements.
  */

    /*****************************************************************************
    *  Initialization
    ****************************************************************************/

	if (!is_initialized_) {

		if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {

			// Extract values from measurement
			float rho = meas_package.raw_measurements_(0);
			float phi = meas_package.raw_measurements_(1);
			float rho_dot = meas_package.raw_measurements_(2);

			// Convert from polar to cartesian coordinates
			float px = rho * cos(phi);
			float py = rho * sin(phi);

			// Initialize state
			x_ << px, py, rho_dot, 0.0, 0.0;
		}
		else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {

			// Extract values from measurement
			float px = meas_package.raw_measurements_(0);
			float py = meas_package.raw_measurements_(1);

			// Initialize state
			x_ << px, py, 0.0, 0.0, 0.0;
		}

		// Initialize state covariance matrix
		P_ = MatrixXd::Identity(n_x_, n_x_);

		// Update last measurement
		previous_timestamp_ = meas_package.timestamp_;

		// Done initializing, no need to predict or update
		is_initialized_ = true;
		return;
	}

	/*****************************************************************************
	*  Prediction
	****************************************************************************/

	// Compute elapsed time from last measurement
	float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;

	// Update last measurement
	previous_timestamp_ = meas_package.timestamp_;

	Prediction(dt);

	/*****************************************************************************
	*  Update
	****************************************************************************/

	if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
		// Radar updates 
		UpdateRadar(meas_package);
	else
		// Laser updates
		UpdateLidar(meas_package);
}

// Generate augmented sigma points matrix 
// Calculate its prediction using input variable dt
// Returm the predicted sigma points matrix
MatrixXd UKF::ComputeSigmaPoints(double delta_t) {

	// Augmented mean state
	VectorXd x_aug = VectorXd(n_aug_);

	x_aug.head(5) = x_;
	x_aug(5) = 0;
	x_aug(6) = 0;

	// Augmented state covariance
	MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

	P_aug.fill(0.0);
	P_aug.topLeftCorner(5, 5) = P_;
	P_aug(5, 5) = std_a_ * std_a_;
	P_aug(6, 6) = std_yawdd_ * std_yawdd_;

	// square root of P
	MatrixXd L = P_aug.llt().matrixL(); 

	// Compute amented sigma point matrix
	MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sigma_points_);
	Xsig_aug.col(0) = x_aug;
	for (int i = 0; i < n_aug_; i++) {
		Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
		Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
	}

	// Actual sigma points' prediction
	for (int i = 0; i < n_sigma_points_; i++) {

		// Auxiliary variables for readability
		double p_x      = Xsig_aug(0, i);
		double p_y      = Xsig_aug(1, i);
		double v        = Xsig_aug(2, i);
		double yaw      = Xsig_aug(3, i);
		double yawd     = Xsig_aug(4, i);
		double nu_a     = Xsig_aug(5, i);
		double nu_yawdd = Xsig_aug(6, i);

		// Sanity check
		if (fabs(p_x) < 0.001 && fabs(p_y) < 0.001) {
			p_x = 0.1;
			p_y = 0.1;
		}

		// Predicted state values
		double px_p, py_p;

		//avoid division by zero
		if (fabs(yawd) > 0.001) {
			px_p = p_x + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
			py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
		}
		else {
			px_p = p_x + v * delta_t * cos(yaw);
			py_p = p_y + v * delta_t * sin(yaw);
		}

		double v_p = v;
		double yaw_p = yaw + yawd * delta_t;
		double yawd_p = yawd;

		// Add noise
		px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
		py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
		v_p = v_p + nu_a * delta_t;

		yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
		yawd_p = yawd_p + nu_yawdd * delta_t;

		//write predicted sigma points into right column
		Xsig_pred_(0, i) = px_p;
		Xsig_pred_(1, i) = py_p;
		Xsig_pred_(2, i) = v_p;
		Xsig_pred_(3, i) = yaw_p;
		Xsig_pred_(4, i) = yawd_p;
	}

	return Xsig_pred_;
}


/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state vector, x_. 
  Predict sigma points, the state, and the state covariance matrix.
  */

    // Compute predicted sigma points matrix
	MatrixXd Xsig_pred_ = ComputeSigmaPoints(delta_t);

	// Predicted state mean (5 x 1) (vectorised form: x = Xsig_pred_ * weights_ could also work)
	VectorXd x = VectorXd(n_x_);
	x.fill(0.0);

	// iterate over sigma points
	for (int i = 0; i < n_sigma_points_; i++) {
		x = x + weights_(i) * Xsig_pred_.col(i);
	}

	//create predicted state covariance matrix (5 x 5)
	MatrixXd P = MatrixXd(n_x_, n_x_);
	P.fill(0.0);

	//iterate over sigma points
	for (int i = 0; i < n_sigma_points_; i++) {

		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x;

		// Normalize angle in range [-PI, PI]
		normalize_between_minus_PI_and_PI(x_diff(3));

		P = P + weights_(i) * x_diff * x_diff.transpose();
	}

	// Update state vector and covariance matrix
	x_ = x;
	P_ = P;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's position. 
  Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  /*****************************************************************************
  *  Prediction
  ****************************************************************************/

    // create sigma points matrix in lidar measurement space
	MatrixXd Zsig_lidar = MatrixXd(n_z_lidar_, n_sigma_points_);
	//MatrixXd Zsig_lidar = Xsig_pred_.block(0, 0, n_z_lidar_, n_sigma_points_);

	for (int i = 0; i < n_sigma_points_; i++) {

		// extract values for better readability
		double p_x = Xsig_pred_(0, i);
		double p_y = Xsig_pred_(1, i);

		Zsig_lidar(0, i) = p_x;
		Zsig_lidar(1, i) = p_y;
	}

	// Predicted measurement mean
	VectorXd z_pred_lidar = VectorXd(n_z_lidar_);
	z_pred_lidar.fill(0.0);
	for (int i = 0; i < n_sigma_points_; i++) {
		z_pred_lidar = z_pred_lidar + weights_(i) * Zsig_lidar.col(i);
	}

	// Predicted measurement covariance matrix
	MatrixXd S_lidar = MatrixXd(n_z_lidar_, n_z_lidar_);
	S_lidar.fill(0.0);
	for (int i = 0; i < n_sigma_points_; i++) { //2n+1 simga points
		
		//residual
		VectorXd z_diff_lidar = Zsig_lidar.col(i) - z_pred_lidar;

		// Normalize angle in range [-PI, PI]
		normalize_between_minus_PI_and_PI(z_diff_lidar(1));

		//while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
		//while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

		S_lidar = S_lidar + weights_(i) * z_diff_lidar * z_diff_lidar.transpose();
	}

	// Handle lidar measurement noise 
	R_lidar_ = MatrixXd(n_z_lidar_, n_z_lidar_);
	
	R_lidar_ << std_laspx_ * std_laspx_, 0, 
		        0, std_laspy_ * std_laspy_;

	S_lidar = S_lidar + R_lidar_;

	/*****************************************************************************
	*  Update
	****************************************************************************/

	// Parse laser measurement
	VectorXd z_lidar = VectorXd(n_z_lidar_);

	z_lidar << meas_package.raw_measurements_[0],
			   meas_package.raw_measurements_[1];


	// Ccompute cross correlation matrix
	MatrixXd Tc_lidar = MatrixXd(n_x_, n_z_lidar_); // 5 x 2
	Tc_lidar.fill(0.0);
	for (int i = 0; i < n_sigma_points_; i++) {

		// Residual
		VectorXd z_diff_lidar = Zsig_lidar.col(i) - z_pred_lidar;

		// Normalize angle
		normalize_between_minus_PI_and_PI(z_diff_lidar(1));

		// State difference
		VectorXd x_diff_lidar = Xsig_pred_.col(i) - x_;

		// Normalize angle
		normalize_between_minus_PI_and_PI(x_diff_lidar(3));

		Tc_lidar = Tc_lidar + weights_(i) * x_diff_lidar * z_diff_lidar.transpose();
	}

	// Compute Kalman gain;
	MatrixXd K_lidar = Tc_lidar * S_lidar.inverse();

	// Residual
	VectorXd z_diff_lidar = z_lidar - z_pred_lidar;

	// Normalize angle
	normalize_between_minus_PI_and_PI(z_diff_lidar(1));

	// Update state mean and covariance matrix
	x_ = x_ + K_lidar * z_diff_lidar;
	P_ = P_ - K_lidar * S_lidar * K_lidar.transpose();

	// Compute NIS for laser sensor
	NIS_laser_ = (meas_package.raw_measurements_ - z_pred_lidar).transpose() * S_lidar.inverse() *
		(meas_package.raw_measurements_ - z_pred_lidar);

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's position. 
  Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  /*****************************************************************************
  *  Prediction
  ****************************************************************************/

    // create sigma points matrix in radar measurement space
	MatrixXd Zsig_radar = MatrixXd(n_z_radar_, n_sigma_points_);

	//transform sigma points into radar measurement space
	for (int i = 0; i < n_sigma_points_; i++) {

		// extract values for better readability
		double p_x = Xsig_pred_(0, i);
		double p_y = Xsig_pred_(1, i);
		double v =   Xsig_pred_(2, i);
		double yaw = Xsig_pred_(3, i);

		double v1 = cos(yaw) * v;
		double v2 = sin(yaw) * v;

		// Measurement model

		//rho
		Zsig_radar(0, i) = sqrt(p_x * p_x + p_y * p_y);

		//phi
		if (fabs(p_y) > 0.001 && fabs(p_x) > 0.001)
			Zsig_radar(1, i) = atan2(p_y, p_x);
		else
			Zsig_radar(1, i) = 0.0;

		//rho_dot
		if (fabs(sqrt(p_x * p_x + p_y * p_y)) > 0.001)
			Zsig_radar(2, i) = (p_x * v1 + p_y * v2) / sqrt(p_x * p_x + p_y * p_y);
		else
			Zsig_radar(2, i) = 0.0;

	}

	// mean predicted measurement 
	VectorXd z_pred_radar = VectorXd(n_z_radar_);
	z_pred_radar.fill(0.0);
	for (int i = 0; i < n_sigma_points_; i++) {
		z_pred_radar = z_pred_radar + weights_(i) * Zsig_radar.col(i);
	}

	// Predicted measurement covariance matrix
	MatrixXd S_radar = MatrixXd(n_z_radar_, n_z_radar_);
	S_radar.fill(0.0);
	for (int i = 0; i < n_sigma_points_; i++) { //2n+1 simga points

		//residual
		VectorXd z_diff_radar = Zsig_radar.col(i) - z_pred_radar;

		// Normalize angle in range [-PI, PI]
		normalize_between_minus_PI_and_PI(z_diff_radar(1));

		S_radar = S_radar + weights_(i) * z_diff_radar * z_diff_radar.transpose();
	}

	// add radar measurement noise covariance matrix
	R_radar_ = MatrixXd(n_z_radar_, n_z_radar_);
	R_radar_ << std_radr_ * std_radr_, 0, 0,
		        0, std_radphi_ * std_radphi_, 0,
		        0, 0, std_radrd_ * std_radrd_;

	S_radar = S_radar + R_radar_;

	/*****************************************************************************
	*  Update
	****************************************************************************/

	// Create vector for incoming radar measurement
	VectorXd z_radar = VectorXd(n_z_radar_);
	z_radar << meas_package.raw_measurements_[0],
		       meas_package.raw_measurements_[1],
		       meas_package.raw_measurements_[2];

	// Compute cross correlation matrix
	MatrixXd Tc_radar = MatrixXd(n_x_, n_z_radar_); 	// 5 x 3
	Tc_radar.fill(0.0);
	for (int i = 0; i < n_sigma_points_; i++) {  // iterate over sigma points

		// Residual
		VectorXd z_diff_radar = Zsig_radar.col(i) - z_pred_radar;

		// Normalize angle in range [-PI, PI]
		normalize_between_minus_PI_and_PI(z_diff_radar(1));

		// State difference
		VectorXd x_diff_radar = Xsig_pred_.col(i) - x_;

		// Normalize angle
		normalize_between_minus_PI_and_PI(x_diff_radar(3));

		Tc_radar = Tc_radar + weights_(i) * x_diff_radar * z_diff_radar.transpose();
	}

	// Compute Kalman gain
	MatrixXd K_radar = Tc_radar * S_radar.inverse();

	// Residual
	VectorXd z_diff_radar = z_radar - z_pred_radar;

	// Normalize angle
	normalize_between_minus_PI_and_PI(z_diff_radar(1));

	// Update state mean and covariance matrix
	x_ = x_ + K_radar * z_diff_radar;
	P_ = P_ - K_radar * S_radar * K_radar.transpose();

	// Compute NIS for radar sensor
	NIS_radar_ = (meas_package.raw_measurements_ - z_pred_radar).transpose() * S_radar.inverse() *
		(meas_package.raw_measurements_ - z_pred_radar);
}


/**
* Helper function to normalize angles in range [-PI, PI]
* @param phi angle to be normalized
*/
void UKF::normalize_between_minus_PI_and_PI(double& phi)
{
	while (phi > M_PI) phi -= 2. * M_PI;
	while (phi <-M_PI) phi += 2. * M_PI;
}

