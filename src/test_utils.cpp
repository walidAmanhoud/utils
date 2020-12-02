#include "Eigen/Eigen"
#include <iostream>
#include "Utils.h"


int main(int argc, char **argv)
{

	ros::init(argc, argv, "test_utils");


	Eigen::Matrix<float,7,1> joints;
	joints.setConstant(0.0f);
	joints << 0.052256,  0.595207, -0.142447,  -1.26305, 0.0819084,   1.52267,   1.52041;
	std::cerr << Utils<float>::getGeometricJacobian(joints) << std::endl;

	
	  Eigen::Matrix4f Hfk;
	  Hfk = Utils<float>::getForwardKinematics(joints);

	  Eigen::Vector3f xEE, xRCM, xTool, xRobotBasis, xTrocar;
	  Eigen::Matrix3f wRb;
	  xRobotBasis << 0.0f, 0.5f, 0.0f;
	  float toolOffset = 0.4f;
	  xTrocar << -0.447232, 0.0650098,  0.204723;



	  xEE = Hfk.block(0,3,3,1) + xRobotBasis;
	  wRb = Hfk.block(0,0,3,3);

	  xRCM = xEE+(xTrocar-xEE).dot(wRb.col(2))*wRb.col(2);
	  xTool = xEE+toolOffset*wRb.col(2);
	  std::cerr  << "bou: " << xRCM.transpose() << std::endl;
	  std::cerr  << "bou: " << xTool.transpose() << std::endl;

}

  
