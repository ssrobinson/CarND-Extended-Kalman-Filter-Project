#include <iostream>
#include "tools.h"

using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
    * Calculate the RMSE here.
  */
	VectorXd rmse(4);
	rmse << 0, 0, 0, 0;

	// check the measurement inputs
	if (estimations.size() == 0 || estimations.size() != ground_truth.size()){
		cout << " Error - estimation or ground_truth data invalid" << endl;
		return rmse;
	}

	// accumulate squared residuals
	for (int = 0; i < estimation.size(); i++){

		VectorXd residual = estimations[i] - ground_truth[i];

		// compute residuals
		residual = residual.array() * residual.array();
		rmse += residual;
	}
	// calculate mean rmse
	rmse = rmse / estimations.size();

	// calculate squared root of mean rmse
	rmse = rmse.array().sqrt();

	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
    * Calculate a Jacobian here.
  */
	MatrixXd Hj(3, 4);
	// state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	// pre-compute Hj terms
	float c1 = px * px + py * py;
	float c2 = sqrt(c1);
	float c3 = (c1 * c2);

	// check for division by zero
	if (fabs(c1) < 0.0001){
		cout << "Error - CalculateJacobian (Division by Zero)" << endl;
		return Hj;
	}

	//compute the Jacobian matrix
	Hj << (px / c2), (py / c2), 0, 0,
		-(py / c1), (px / c1), 0, 0,
		py * (vx * py - vy * px) / c3, px * (px * vy - py * vx) / c3, px / c2, py / c2;

	return Hj;
}
