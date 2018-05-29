#ifndef __UTILS_H__
#define __UTILS_H__

#include "ros/ros.h"
#include "Eigen/Eigen"

class Utils 
{
	public:

		// Class constructor
		Utils();

    static Eigen::Vector4f quaternionProduct(Eigen::Vector4f q1, Eigen::Vector4f q2);

    static Eigen::Matrix3f getSkewSymmetricMatrix(Eigen::Vector3f input);

    static Eigen::Vector4f rotationMatrixToQuaternion(Eigen::Matrix3f R);

  	static Eigen::Matrix3f quaternionToRotationMatrix(Eigen::Vector4f q);

		static void quaternionToAxisAngle(Eigen::Vector4f q, Eigen::Vector3f &axis, float &angle);

  	static Eigen::Vector4f slerpQuaternion(Eigen::Vector4f q1, Eigen::Vector4f q2, float t);

  	static float smoothRise(float x, float a, float b);

		static float smoothFall(float x, float a, float b);

		static float smoothRiseFall(float x, float a, float b, float c, float d);
};


#endif
