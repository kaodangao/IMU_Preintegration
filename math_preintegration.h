#ifndef _MATH_PREINTEGRATION_H
#define _MATH_PREINTEGRATION_H

#include <eigen3/Eigen/Eigen>

typedef Eigen::Matrix<double, 15, 1> Vector15d;
typedef Eigen::Matrix<double, 9, 1> Vector9d;
typedef Eigen::Matrix<double, 6, 1> Vector6d;
typedef Eigen::Matrix<double, 4, 1> Vector4d;
typedef Eigen::Matrix<double, 3, 1> Vector3d;
typedef Eigen::Matrix<double, 3, 3> Matrix3d;
typedef Eigen::Matrix<double, 6, 6> Matrix6d;
typedef Eigen::Matrix<double, 9, 9> Matrix9d;
typedef Eigen::Matrix<double, 6, 9> Matrix69;
typedef Eigen::Matrix<double, 9, 3> Matrix93;

Matrix3d hat( Vector3d v){
    Matrix3d V;
    V.fill(0.);
    V(0,1) = -v(2);
    V(0,2) = v(1);
    V(1,2) = -v(0);
    V(1,0) = v(2);
    V(2,0) = -v(1);
    V(2,1) = v(0);
    return V;
};

Matrix3d jacobian_right( Vector3d phi){
    Matrix3d Jr;
    float phi_abs = phi.norm();
    Matrix3d phi_hat = hat(phi);
    
    if(phi_abs < 0.0001){
        Jr = Matrix3d::Identity()-0.5*phi_hat+phi_hat*phi_hat/6;
    }
    else{
        Jr = Matrix3d::Identity()-(1-cos(phi_abs))/(phi_abs*phi_abs)*phi_hat+(phi_abs-sin(phi_abs))/(phi_abs*phi_abs*phi_abs)*phi_hat*phi_hat;
    }
    return Jr;
};

Matrix3d jacobian_right_inv( Vector3d phi){
    Matrix3d Jr_inv;
    float phi_abs = phi.norm();
    Matrix3d phi_hat = hat(phi);

    if(phi_abs < 0.0001){
        Jr_inv = Matrix3d::Identity()+0.5*phi_hat;
    }
    else{
        Jr_inv = Matrix3d::Identity()+0.5*phi_hat+(1/(phi_abs*phi_abs)-(1+cos(phi_abs))/(2*phi_abs*sin(phi_abs)))*phi_hat*phi_hat;
    }
    return Jr_inv;
};

Matrix3d Exp( Vector3d phi){

    float theta = phi.norm();
    
    if(theta == 0){
        return Matrix3d::Identity();
    }
    else{
        Vector3d v = sin(0.5*theta)/theta*phi;
        Eigen::Quaterniond Q(cos(0.5*theta),v(0),v(1),v(2));
        return Q.toRotationMatrix();
    }
};

Vector3d Log( Matrix3d R ){

    Eigen::Quaterniond Q(R);
    Vector3d v(Q.x(),Q.y(),Q.z());
    double n = v.norm();

    if(n==0 or Q.w()==0){
        return Vector3d::Zero();
    }
    else{
        return 2*atan(n/Q.w())/n*v;       
    }

}

#endif
