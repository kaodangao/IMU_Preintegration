#ifndef _MATH_PREINTEGRATION_H
#define _MATH_PREINTEGRATION_H

#include <Eigen/Eigen>

typedef Eigen::Matrix<double, 15, 1> Vector15d;
typedef Eigen::Matrix<double, 9, 1> Vector9d;
typedef Eigen::Matrix<double, 6, 1> Vector6d;
typedef Eigen::Matrix<double, 3, 1> Vector3d;
typedef Eigen::Matrix<double, 3, 3> Matrix3d;
typedef Eigen::Matrix<double, 6, 6> Matrix6d;
typedef Eigen::Matrix<double, 9, 9> Matrix9d;
typedef Eigen::Matrix<double, 6, 9> Matrix69;


Eigen::Matrix3d hat( Eigen::Vector3d v){
    Eigen::Matrix3d V;
    V.fill(0.);
    V(0,1) = -v(2);
    V(0,2) = v(1);
    V(1,2) = -v(0);
    V(1,0) = v(2);
    V(2,0) = -v(1);
    V(2,1) = v(0);
    return V;
};

Eigen::Matrix3d jacobian_right( Eigen::Vector3d phi){
    Eigen::Matrix3d Jr;
    float phi_abs = phi.norm();
    Eigen::Matrix3d phi_hat = hat(phi);
    if(phi_abs < 0.00001){
        return Eigen::Matrix3d::Identity()-0.5*phi_hat;
    }
    Jr = Eigen::Matrix3d::Identity()-(1-cos(phi_abs))/(phi_abs*phi_abs)*phi_hat+(phi_abs-sin(phi_abs))/(phi_abs*phi_abs*phi_abs)*phi_hat*phi_hat;
    return Jr;
};

Eigen::Matrix3d jacobian_right_inv( Eigen::Vector3d phi){
    Eigen::Matrix3d Jr_inv;
    float phi_abs = phi.norm();
    Eigen::Matrix3d phi_hat = hat(phi);
    if(phi_abs < 0.00001){
        return Eigen::Matrix3d::Identity()+0.5*phi_hat;
    }
    Jr_inv = Eigen::Matrix3d::Identity()+0.5*phi_hat+(1/(phi_abs*phi_abs)-(1+cos(phi_abs))/(2*phi_abs*sin(phi_abs)))*phi_hat*phi_hat;
    return Jr_inv;
};

Eigen::Matrix3d exp( Eigen::Vector3d phi){
    Eigen::Matrix3d R;
    float phi_abs = phi.norm();
    Eigen::Matrix3d phi_hat = hat(phi);
    if(phi_abs < 0.00001){
        return Eigen::Matrix3d::Identity()+phi_hat;
    }
    R = Eigen::Matrix3d::Identity()+sin(phi_abs)/phi_abs*phi_hat+(1-cos(phi_abs))/(phi_abs*phi_abs)*phi_hat*phi_hat;
    return R;
};

Eigen::Vector3d Log( Eigen::Matrix3d r){
    double d =  0.5*(r(0,0)+r(1,1)+r(2,2)-1);
    Eigen::Vector3d omega;
    Eigen::Vector3d dR;
    dR << r(2,1)-r(1,2),r(0,2)-r(2,0),r(1,0)-r(0,1);

    if (d>0.99999)
    {
        omega=0.5*dR;
    }
    else
    {
        double theta = acos(d);
        omega = theta/(2*sqrt(1-d*d))*dR;
    }

    return omega;
}

#endif
