#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
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

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.8;

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

  /**
  TODO
  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  //State dimension.
  n_x_ = 5;
  //Augmented state.
  n_aug_ = 7;
  //Lambda constant for sigma point generation.
  lambda_ = 3;

  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_ + 1);
  weights_ = VectorXd(2*n_aug_ + 1);
  
  R_radar_ = MatrixXd(3, 3);
  R_lidar_ = MatrixXd(2, 2);

  P_ << std_laspx_*std_laspx_, 0, 0, 0, 0,
		0, std_laspy_*std_laspy_, 0, 0, 0,
		0, 0, 1, 0, 0,
		0, 0, 0, 1, 0,
		0, 0, 0, 0, 1;
      
  R_lidar_ << std_laspx_*std_laspx_, 0,
                 0, std_laspy_*std_laspy_;
  
  R_radar_ << std_radr_*std_radr_, 0, 0,
	0, std_radphi_*std_radphi_, 0,
	0, 0,std_radrd_*std_radrd_;
  
  NIS_laser_ = 0.0;
  NIS_radar_ = 0.0;
}

UKF::~UKF() {}


/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  if(!is_initialized_){
    double px, py, v, rhox, rhoy;
     if(meas_package.sensor_type_ == MeasurementPackage::LASER){
         //Handle lasar measurement.
         px = meas_package.raw_measurements_(0);
         py = meas_package.raw_measurements_(1);
         v =  4;
         rhox = 0.5;
         rhoy = 0;
         x_ << px, py, v, rhox, rhoy;
     }
     else{
         //Handle radar measurement.
         double rho = meas_package.raw_measurements_(0);
         double phi = meas_package.raw_measurements_(1);
         double rho_rate = meas_package.raw_measurements_(2);
         px = rho * cos(phi);
         py = rho * sin(phi);
         v = 4;
         rhox = rho_rate * cos(phi);
         rhoy = rho_rate * sin(phi);
         x_ << px, py, v, rhox, rhoy;
     }

	 time_us_ = meas_package.timestamp_;
     is_initialized_ = true;
     return;
  }

  double time_delta = (meas_package.timestamp_ - time_us_)/1000000.0;
  Prediction(time_delta);

  if(meas_package.sensor_type_ == MeasurementPackage::LASER){
      UpdateLidar(meas_package);
  }
  else{
      UpdateRadar(meas_package);
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
	int n_z = 2; //Lidar having measurement dimension two.
	int sigma_aug_size = 2*n_aug_ + 1;
	MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_ + 1);
	Zsig.fill(0.0);
	//Mean predicted measurement.
	VectorXd z_pred = MatrixXd(n_z, n_z);
	z_pred.fill(0.0);
	//Measurement covariance matrix.
	MatrixXd S = MatrixXd(n_z, n_z);
	S.fill(0.0);

    double px;
    double py;
	for(auto i = 0; i < sigma_aug_size; i++){
        px = Xsig_pred_(0, i);
        py = Xsig_pred_(1, i);

        Zsig(0, i) = px;
        Zsig(1, i) = py;
        z_pred += weights_(i) * Zsig.col(i);
	}

    //Measurement covariance matrix.
    for(auto i = 0; i < sigma_aug_size; i ++){
        VectorXd z_diff = Zsig.col(i) - z_pred;
        S+= weights_(i) * z_diff * z_diff.transpose();
    }

    // S + R : Add noise.
    S = S + R_lidar_;

    //Input vector.
    VectorXd z = VectorXd(n_z);
    z.fill(0.0);
    double input_px = meas_package.raw_measurements_(0);
    double input_py = meas_package.raw_measurements_(1);

    z << input_px, input_py;

    //Calculate cross correlation.
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    Tc.fill(0.0);
    for(auto i = 0; i < sigma_aug_size; i++){
        VectorXd x_diff = Xsig_pred_.col(i) - x_;

        if(x_diff(3) > M_PI){
            x_diff(3) -= 2.0 * M_PI;
        }
        else if(x_diff(3) < - M_PI){
            x_diff(3) += 2.0 * M_PI;
        }

        VectorXd z_diff = Zsig.col(i) - z_pred;
        Tc += weights_(i) * x_diff * z_diff.transpose();
    }

    // residual.
    VectorXd z_diff = z - z_pred;
    //Kalman gain.
    MatrixXd K = Tc * S.inverse();
    //Update mean and covariance.
    x_ = x_ + K * z_diff;
    P_ = P_ - K * S * K.transpose();

    //Calculate NIS.
    NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
	int n_z = 3; //Radar having measurement dimension two.
	int sigma_aug_size = 2*n_aug_ + 1;
	MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_ + 1);
	Zsig.fill(0.0);
	//Mean predicted measurement.
	VectorXd z_pred = MatrixXd(n_z, n_z);
	z_pred.fill(0.0);
	//Measurement covariance matrix.
	MatrixXd S = MatrixXd(n_z, n_z);
	S.fill(0.0);

    double px, py, v, yaw;
    double rho = 0.0;
    double phi = 0.0;
    double rho_rate = 0.0;
	for(auto i = 0; i < sigma_aug_size; i++){
        px  =  Xsig_pred_(0, i);
        py  =  Xsig_pred_(1, i);
        v   =  Xsig_pred_(2, i);
        yaw =  Xsig_pred_(3, i);

        rho = sqrt(px*px + py*py);
        phi = atan2(py, px);
        rho_rate = (px*cos(yaw)*v + py*sin(yaw)*v)/rho;

        Zsig(0, i) = rho;
        Zsig(1, i) = phi;
        Zsig(2, i) = rho_rate;
        z_pred += weights_(i) * Zsig.col(i);
	}

    // Measurement covariance matrix.
    for(auto i = 0; i < sigma_aug_size; i ++){
        VectorXd z_diff = Zsig.col(i) - z_pred;
        if(z_diff(1) > M_PI){
            z_diff(1) -= 2.0 * M_PI;
        }
        else{
            z_diff(1) += 2.0 * M_PI;
        }
        S+= weights_(i) * z_diff * z_diff.transpose();
    }
    // S = S + R: Radar update.
    S = S + R_radar_;
    //Input vector.
    VectorXd z = VectorXd(n_z);

    double input_rho = meas_package.raw_measurements_(0);
    double input_phi = meas_package.raw_measurements_(1);
    double input_rho_rate = meas_package.raw_measurements_(2);

    z << input_rho, input_phi, input_rho_rate;

    //Cross correlation Matrix.
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    Tc.fill(0.0);

    // Calculate cross correlation matrix.
    for(auto i = 0; i < sigma_aug_size; i++){
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        VectorXd z_diff = Zsig.col(i) - z_pred;

        // normalize angles.
		if(x_diff(3) > M_PI){
   			x_diff(3)-=2.*M_PI;
		}
		else if(x_diff(3) < -M_PI){
    		x_diff(3)+=2.*M_PI;
		}
		
		if(z_diff(1) > M_PI){
   			z_diff(1)-=2.*M_PI;
		}
		else if(z_diff(1) < -M_PI){
    		z_diff(1)+=2.*M_PI;
		}
        Tc += weights_(i) * x_diff * z_diff.transpose();
    }
    //residual.
    VectorXd z_diff = z - z_pred;

    if(z_diff(1) > M_PI){
        z_diff(1) -= 2.0 * M_PI;
    }
    else if(z_diff(1) < -M_PI){
        z_diff(1) += 2.0 * M_PI;
    }

    //Kalman gain.
    MatrixXd K = Tc * S.inverse();
    //Update mean and covariance.
    x_ = x_ + K * z_diff;
    P_ = P_ - K * S * K.transpose();

    //NIS Update.
    NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
}

/*Generate sigma points with consideration for noise.*/
void UKF::GenerateSigmaPoints(MatrixXd* XSig_out){

    double lambda_ = 3 - n_aug_;
    MatrixXd P = MatrixXd(5, 5);
    P.fill(0.0);

    VectorXd x = VectorXd(n_x_);
    VectorXd x_aug = VectorXd(n_aug_);
    MatrixXd P_aug = MatrixXd(7, 7);
    
    x_aug.head(5) = x;
    x_aug(5) = 0;
    x_aug(6) = 0;

    P_aug.fill(0.0);
    P_aug.topLeftCorner(5, 5) = P;
    P_aug(5, 5) = std_a_;
    P_aug(6, 6) = std_yawdd_;

    MatrixXd L = P_aug.llt().matrixL();
    MatrixXd X_sig = MatrixXd(n_aug_, 2*n_aug_ + 1);
    X_sig.col(0) = x_aug;

    for (int i = 0; i < n_aug_; i++){
        X_sig.col(i+1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i); 
        X_sig.col(i+1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i); 
    }

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t){
    MatrixXd XSig_aug = MatrixXd(n_aug_, 2*n_aug_ + 1);
    MatrixXd XSig_pred = MatrixXd(n_x_, 2*n_aug_ + 1);

    GenerateSigmaPoints(&XSig_aug);
    double px_p, py_p;

    for(int i; i < 2*n_aug_+1; i++){
        double p_x = XSig_aug(0, i);
        double p_y = XSig_aug(1, i);
        double v = XSig_aug(2, i);
        double yaw = XSig_aug(3, i);
        double yawd = XSig_aug(4, i);

        double nu_a = XSig_aug(5, i);
        double nu_yawdd = XSig_aug(6, i);

        if(fabs(yawd) > 0.001){
            px_p = p_x + v/yawd * (sin(yaw + yawd*delta_t) - sin(yaw));
            py_p = p_y + v/yawd * (cos(yaw) - cos(yaw + yawd*delta_t));
        }
        else{
            px_p = p_x + v*delta_t*cos(yaw);
            py_p = p_y + v*delta_t*sin(yaw);
        }

        double v_p = v;
        double yaw_p = yaw + yawd*delta_t;
        double yawd_p = yawd;

        px_p = px_p + 0.5*nu_a*delta_t*delta_t*cos(yaw);
        py_p = py_p + 0.5*nu_a*delta_t*delta_t*sin(yaw);
        v_p = v_p + nu_a*delta_t;

        yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
        yawd_p = yawd_p + nu_yawdd*delta_t;

        XSig_pred(0, i) = px_p;
        XSig_pred(1, i) = py_p;
        XSig_pred(2, i) = v_p;
        XSig_pred(3, i) = yaw_p;
        XSig_pred(4, i) = yawd_p;
    }
    Xsig_pred_ = XSig_pred;
	PredictMeanCovariance();
}

//Predict mean and covariance for sigma points.
void UKF::PredictMeanCovariance(void){

    MatrixXd XSig_pred = MatrixXd(n_x_, 2*n_aug_ + 1);
    VectorXd weights = VectorXd(2*n_aug_ + 1);
    
    VectorXd x_pred = VectorXd(n_x_);
    MatrixXd P_pred = MatrixXd(n_x_, n_x_);
    GenerateSigmaPoints(&XSig_pred);

    weights(0) = lambda_/(lambda_ + n_aug_);
    for(int i = 1; i < 2*n_aug_ + 1; i++){
        weights(i) = 0.5/(n_aug_ + lambda_);
    }
    
    x_pred.fill(0.0);
    P_pred.fill(0.0);

    for(int i = 0; i < 2*n_aug_ + 1; i++){
        x_pred = x_pred + weights(i) * XSig_pred.col(i);
    }
    
    for(int i = 0; i < 2*n_aug_ + 1; i++){

        VectorXd x_diff = XSig_pred.col(i) - x_pred;
        
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

        P_pred = P_pred + weights(i) * x_diff * x_diff.transpose(); 
    }

    x_ = x_pred;
    P_ = P_pred;
}
