#ifndef __UTILS_H__
#define __UTILS_H__

#include "ros/ros.h"
#include "Eigen/Eigen"

template<typename T = float>
class Utils 
{
	public:

	// Class constructor
	Utils();

    static Eigen::Matrix<T,4,1> quaternionProduct(Eigen::Matrix<T,4,1> q1, Eigen::Matrix<T,4,1> q2);

    static Eigen::Matrix<T,3,3> getSkewSymmetricMatrix(Eigen::Matrix<T,3,1> input);

    static Eigen::Matrix<T,4,1> rotationMatrixToQuaternion(Eigen::Matrix<T,3,3> R);

  	static Eigen::Matrix<T,3,3> quaternionToRotationMatrix(Eigen::Matrix<T,4,1> q);

	static void quaternionToAxisAngle(Eigen::Matrix<T,4,1> q, Eigen::Matrix<T,3,1> &axis, T &angle);

  	static Eigen::Matrix<T,4,1> slerpQuaternion(Eigen::Matrix<T,4,1> q1, Eigen::Matrix<T,4,1> q2, T t);

  	static T smoothRise(T x, T a, T b);

	static T smoothFall(T x, T a, T b);

	static T smoothRiseFall(T x, T a, T b, T c, T d);

	static T deadZone(T x, T a, T b);		

	static T wrapToZero(T x, T a, T b);		
};

#endif
