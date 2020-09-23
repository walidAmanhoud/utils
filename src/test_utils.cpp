#include "Eigen/Eigen"
#include <iostream>
#include "Utils.h"


int main(int argc, char **argv)
{

	ros::init(argc, argv, "test_utils");


	Eigen::Matrix<float,7,1> joints;
	joints.setConstant(0.0f);
	joints << -M_PI/4.0f, M_PI/6.0f, M_PI/4.0f, -M_PI/2.0f, M_PI/3.0f, M_PI/2.0f, -M_PI/3.0f;
	std::cerr << Utils<float>::getGeometricJacobian(joints) << std::endl;
}

  
