#include "Utils.h"

template<typename T>
Utils<T>::Utils()
{

}


template<typename T>
Eigen::Matrix<T,4,1> Utils<T>::quaternionProduct(Eigen::Matrix<T,4,1> q1, Eigen::Matrix<T,4,1> q2)
{
  Eigen::Matrix<T,4,1> q;
  q(0) = q1(0)*q2(0)-(q1.segment(1,3)).dot(q2.segment(1,3));
  Eigen::Matrix<T,3,1> q1Im = (q1.segment(1,3));
  Eigen::Matrix<T,3,1> q2Im = (q2.segment(1,3));
  q.segment(1,3) = q1(0)*q2Im+q2(0)*q1Im+q1Im.cross(q2Im);

  return q;
}


template<typename T>
Eigen::Matrix<T,3,3> Utils<T>::getSkewSymmetricMatrix(Eigen::Matrix<T,3,1> input)
{
  Eigen::Matrix<T,3,3> output;

  output << 0.0f, -input(2), input(1),
            input(2), 0.0f, -input(0),
            -input(1), input(0), 0.0f;

  return output;
}


template<typename T>
Eigen::Matrix<T,4,1> Utils<T>::rotationMatrixToQuaternion(Eigen::Matrix<T,3,3> R)
{
  Eigen::Matrix<T,4,1> q;

  float r11 = R(0,0);
  float r12 = R(0,1);
  float r13 = R(0,2);
  float r21 = R(1,0);
  float r22 = R(1,1);
  float r23 = R(1,2);
  float r31 = R(2,0);
  float r32 = R(2,1);
  float r33 = R(2,2);

  float tr = r11+r22+r33;
  float tr1 = r11-r22-r33;
  float tr2 = -r11+r22-r33;
  float tr3 = -r11-r22+r33;

  if(tr>0)
  {  
    q(0) = sqrt(1.0f+tr)/2.0f;
    q(1) = (r32-r23)/(4.0f*q(0));
    q(2) = (r13-r31)/(4.0f*q(0));
    q(3) = (r21-r12)/(4.0f*q(0));
  }
  else if((tr1>tr2) && (tr1>tr3))
  {
    q(1) = sqrt(1.0f+tr1)/2.0f;
    q(0) = (r32-r23)/(4.0f*q(1));
    q(2) = (r21+r12)/(4.0f*q(1));
    q(3) = (r31+r13)/(4.0f*q(1));
  }     
  else if((tr2>tr1) && (tr2>tr3))
  {   
    q(2) = sqrt(1.0f+tr2)/2.0f;
    q(0) = (r13-r31)/(4.0f*q(2));
    q(1) = (r21+r12)/(4.0f*q(2));
    q(3) = (r32+r23)/(4.0f*q(2));
  }
  else
  {
    q(3) = sqrt(1.0f+tr3)/2.0f;
    q(0) = (r21-r12)/(4.0f*q(3));
    q(1) = (r31+r13)/(4.0f*q(3));
    q(2) = (r32+r23)/(4.0f*q(3));        
  }

  return q;
}


template<typename T>
Eigen::Matrix<T,3,3> Utils<T>::quaternionToRotationMatrix(Eigen::Matrix<T,4,1> q)
{
  Eigen::Matrix<T,3,3> R;

  T q0 = q(0);
  T q1 = q(1);
  T q2 = q(2);
  T q3 = q(3);

  R(0,0) = q0*q0+q1*q1-q2*q2-q3*q3;
  R(1,0) = 2.0f*(q1*q2+q0*q3);
  R(2,0) = 2.0f*(q1*q3-q0*q2);

  R(0,1) = 2.0f*(q1*q2-q0*q3);
  R(1,1) = q0*q0-q1*q1+q2*q2-q3*q3;
  R(2,1) = 2.0f*(q2*q3+q0*q1);

  R(0,2) = 2.0f*(q1*q3+q0*q2);
  R(1,2) = 2.0f*(q2*q3-q0*q1);
  R(2,2) = q0*q0-q1*q1-q2*q2+q3*q3;  

  return R;
}


template<typename T>
void Utils<T>::quaternionToAxisAngle(Eigen::Matrix<T,4,1> q, Eigen::Matrix<T,3,1> &axis, T &angle)
{
  if((q.segment(1,3)).norm() < 1e-3f)
  {
    axis = q.segment(1,3);
  }
  else
  {
    axis = q.segment(1,3)/(q.segment(1,3)).norm();
    
  }

  angle = 2*std::acos(q(0));
}

template<typename T>
Eigen::Matrix<T,4,1> Utils<T>::axisAngleToQuaterion(Eigen::Matrix<T,3,1> axis, T angle)
{
  Eigen::Matrix<T,4,1> q;
  q(0) = std::cos(angle/2);
  q(1) = axis(0)*std::sin(angle/2);
  q(2) = axis(1)*std::sin(angle/2);
  q(3) = axis(2)*std::sin(angle/2);
  return q;
}


template<typename T>
Eigen::Matrix<T,4,1> Utils<T>::slerpQuaternion(Eigen::Matrix<T,4,1> q1, Eigen::Matrix<T,4,1> q2, T t)
{

  Eigen::Matrix<T,4,1> q;

  // Change sign of q2 if dot product of the two quaterion is negative => allows to interpolate along the shortest path
  if(q1.dot(q2)<0.0f)
  {   
    q2 = -q2;
  }

  T dotProduct = q1.dot(q2);
  if(dotProduct > 1.0f)
  {
    dotProduct = 1.0f;
  }
  else if(dotProduct < -1.0f)
  {
    dotProduct = -1.0f;
  }

  T omega = acos(dotProduct);

  if(std::fabs(omega)<FLT_EPSILON)
  {
    q = q1.transpose()+t*(q2-q1).transpose();
  }
  else
  {
    q = (std::sin((1-t)*omega)*q1+std::sin(t*omega)*q2)/std::sin(omega);
  }

  return q;
}

template<typename T>
Eigen::Matrix<T,4,1> Utils<T>::slerpQuaternion(Eigen::Matrix<T,4,1>* q, Eigen::Matrix<T,Eigen::Dynamic,1> t, int size)
{

  if(size==1)
  {
    return q[0];
  }
  else
  {
    T sum = 0.0f;
    for(int k = 0; k <size; k++)
    {
      sum+=t(k);
    }
    if(sum<FLT_EPSILON)
    {
      return slerpQuaternion(slerpQuaternion(q,t,size-1),q[size-1],0.0f);
    }
    else
    {
      return slerpQuaternion(slerpQuaternion(q,t,size-1),q[size-1],t(size-1)/sum);
    }
  }
}


template<typename T>
Eigen::Matrix<T,3,3> Utils<T>::rodriguesRotation(Eigen::Matrix<T,3,1> v1, Eigen::Matrix<T,3,1> v2)
{
  // Compute rotation error between current orientation and plane orientation using Rodrigues' law
  v1.normalize();
  v2.normalize();

  Eigen::Matrix<T,3,1> w;
  w = v1.cross(v2);
  float c = v1.dot(v2);  
  float s = w.norm();
  w /= s;
  
  Eigen::Matrix<T,3,3> K;
  K << getSkewSymmetricMatrix(w);

  Eigen::Matrix<T,3,3> Re;
  if(fabs(s)< FLT_EPSILON)
  {
    Re = Eigen::Matrix<T,3,3>::Identity();
  }
  else
  {
    Re = Eigen::Matrix<T,3,3>::Identity()+s*K+(1-c)*K*K;
  }

  return Re;
}

template<typename T>
Eigen::Matrix<T,3,1> Utils<T>::quaternionToAngularVelocity(Eigen::Matrix<T,4,1> q1, Eigen::Matrix<T,4,1> q2, T gain)
{
  Eigen::Matrix<T,4,1> q1I, wq;
  q1I(0) = q1(0);
  q1I.segment(1,3) = -q1.segment(1,3);
  wq = 2.0f*gain*quaternionProduct(q2-q1,q1I);
  
  return wq.segment(1,3);
}


template<typename T>
Eigen::Matrix<T,3,3> Utils<T>::orthogonalProjector(Eigen::Matrix<T,3,1> v)
{ 
  return Eigen::Matrix<T,3,3>::Identity()-v*v.transpose();
}


template<typename T>
T Utils<T>::smoothRise(T x, T a, T b)
{
  T y; 
  if(x<a)
  {
    y = 0.0f;
  }
  else if(x>b)
  {
    y = 1.0f;
  }
  else
  {
    y = (1.0f+sin(M_PI*(x-a)/(b-a)-M_PI/2.0f))/2.0f;
  }

  return y;
}


template<typename T>
T Utils<T>::smoothFall(T x, T a, T b)
{
  return 1.0f-smoothRise(x,a,b);
}


template<typename T>
T Utils<T>::smoothRiseFall(T x, T a, T b, T c, T d)
{
  return smoothRise(x,a,b)*smoothFall(x,c,d);
}


template<typename T>
T Utils<T>::deadZone(T x, T a, T b)
{
  if(x < b && x > a)
  {
    return 0.0f;
  }
  else
  {
    return x;
  }
}

template<typename T>
T Utils<T>::wrapToZero(T x, T a, T b)
{
  if(x < b && x > a)
  {
    return x;
  }
  else
  {
    return 0.0f;
  }
}


template class Utils<float>;
template class Utils<double>;
