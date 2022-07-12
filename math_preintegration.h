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
    
    if(phi_abs < 0.001){
        Jr = Matrix3d::Identity()-0.5*phi_hat+phi_hat*phi_hat/6;
    }
    else{
        Jr = Matrix3d::Identity()-(1-cos(phi_abs))/(phi_abs*phi_abs)*phi_hat+(phi_abs-sin(phi_abs))/(phi_abs*phi_abs*phi_abs)*phi_hat*phi_hat;
    }
    // if(phi_abs < 0.001){
    //     Jr = Matrix3d::Identity()-0.5*phi_abs*phi_hat;
    // }
    // else{
    //     Jr = sin(phi_abs)/phi_abs*Matrix3d::Identity()+(1- sin(phi_abs)/phi_abs)*phi*phi.transpose()+((cos(phi_abs)-1)/phi_abs)*phi_hat;
    // }
    return Jr;
};

Matrix3d jacobian_right_inv( Vector3d phi){
    Matrix3d Jr_inv;
    float phi_abs = phi.norm();
    Matrix3d phi_hat = hat(phi);

    if(phi_abs < 0.001){
        Jr_inv = Matrix3d::Identity()+0.5*phi_hat;
    }
    else{
        Jr_inv = Matrix3d::Identity()+0.5*phi_hat+(1/(phi_abs*phi_abs)-(1+cos(phi_abs))/(2*phi_abs*sin(phi_abs)))*phi_hat*phi_hat;
    }
    return Jr_inv;
};

Matrix3d Exp( Vector3d phi){
    Matrix3d R;
    
    float phi_abs = phi.norm();
    if(phi_abs<0.001){

        R = Matrix3d::Identity()+hat(phi);
    }
    else{
        R = Matrix3d::Identity()+sin(phi_abs)/phi_abs*hat(phi)+(1-cos(phi_abs))/(phi_abs*phi_abs)*hat(phi)*hat(phi);

    }

    return R;
};

Vector3d Log( Matrix3d r){
    double d =  0.5*(r(0,0)+r(1,1)+r(2,2)-1);

    Vector3d omega;
    Vector3d dR;
    dR << r(2,1)-r(1,2),r(0,2)-r(2,0),r(1,0)-r(0,1);

    if (d>0.999 or d<-0.999){
        omega=0.5*dR;
    }
    else{
        double theta = acos(d);
        omega = theta/(2*sqrt(1-d*d))*dR;
    }
    return omega;
    // Eigen::AngleAxisd rotation_vector(r);
    // Vector3d omega;
    // omega = rotation_vector.angle()*rotation_vector.axis();
    // return omega;


}

#endif
