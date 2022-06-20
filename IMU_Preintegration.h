#include <Eigen/Eigen>

#include "math_preintegration.h"

class IMU_Preintegration{
    public:

    IMU_Preintegration(Vector3d acc_bias, Vector3d gyro_bias, double acc_noise_sigma, double gyro_noise_sigma){
        delta_tij = 0;
        g << 0, 0, -9.8;

        delta_p = Vector3d::Zero();
        delta_v = Vector3d::Zero();
        delta_r = Matrix3d::Identity();

        ba = Vector3d::Zero(); 
        bg = Vector3d::Zero();

        pd_r_bg = Matrix3d::Zero();
        pd_p_bg = Matrix3d::Zero();
        pd_p_ba = Matrix3d::Zero();
        pd_v_bg = Matrix3d::Zero();
        pd_v_ba = Matrix3d::Zero();

        cov_acc = acc_noise_sigma*acc_noise_sigma*Matrix3d::Identity();
        cov_gyro = gyro_noise_sigma*gyro_noise_sigma*Matrix3d::Identity();
        cov_residual = Matrix9d::Identity();
    }

    void Calculate(Eigen::Vector3d acc, Eigen::Vector3d gyro, double delta_t){
        delta_tij += delta_t;
        Eigen::Matrix3d dR = exp((gyro-bg)*delta_t);

        Eigen::Matrix<double, 9, 9> A = Eigen::Matrix<double, 9, 9>::Identity();
        A.block(0,0,3,3) = dR.transpose();
        A.block(6,0,3,3) = -delta_r*hat(acc-ba)*delta_t;
        A.block(3,0,3,3) = -0.5*delta_r*hat(acc-ba)*delta_t*delta_t;
        A.block(3,6,3,3) = Eigen::Matrix<double, 3, 3>::Identity()*delta_t;

        Eigen::Matrix<double, 9, 3> B = Eigen::Matrix<double, 9, 3>::Zero();
        B.block(3,0,3,3) = 0.5*delta_r*delta_t*delta_t;
        B.block(6,0,3,3) = delta_r*delta_t;

        Eigen::Matrix<double, 9, 3> C = Eigen::Matrix<double, 9, 3>::Zero();
        C.block(0,0,3,3) = jacobian_right((gyro-bg)*delta_t)*delta_t;
        cov_residual = A*cov_residual*A.transpose() + B*cov_acc*B.transpose() + C*cov_gyro*C.transpose();

        pd_p_ba += pd_v_ba*delta_t - 0.5*delta_r*delta_t*delta_t;
        pd_p_bg += pd_v_bg*delta_t - 0.5*delta_r*hat(acc-ba)*pd_r_bg*delta_t*delta_t;
        pd_v_ba -= delta_r*delta_t;
        pd_v_bg -= delta_r*hat(acc-ba)*pd_r_bg*delta_t;
        pd_r_bg = dR.transpose()*pd_r_bg - jacobian_right((gyro-bg)*delta_t)*delta_t;

        
        delta_p += delta_v*delta_t + 0.5*delta_r*(acc-ba)*delta_t*delta_t;
        delta_v += delta_r*(acc-ba)*delta_t;
        delta_r = delta_r*dR;
    }

    double delta_tij_read(){return delta_tij;}
    Eigen::Vector3d delta_p_read(){return delta_p;}
    Eigen::Vector3d delta_v_read(){return delta_v;}
    Eigen::Matrix3d delta_r_read(){return delta_r;}
    Eigen::Vector3d ba_read(){return ba;}
    Eigen::Vector3d bg_read(){return bg;}
    Eigen::Matrix3d pd_r_bg_read(){return pd_r_bg;}
    Eigen::Matrix3d pd_p_bg_read(){return pd_p_bg;}
    Eigen::Matrix3d pd_p_ba_read(){return pd_p_ba;}
    Eigen::Matrix3d pd_v_bg_read(){return pd_v_bg;}
    Eigen::Matrix3d pd_v_ba_read(){return pd_v_ba;}
    Eigen::Matrix<double, 9, 9> cov_residual_read(){return cov_residual;};
    

    private:

    double delta_tij;
    Eigen::Vector3d g;

    Eigen::Vector3d delta_p;
    Eigen::Vector3d delta_v;
    Eigen::Matrix3d delta_r;

    Eigen::Vector3d ba;
    Eigen::Vector3d bg;

    Eigen::Matrix3d pd_r_bg;
    Eigen::Matrix3d pd_p_bg;
    Eigen::Matrix3d pd_p_ba;
    Eigen::Matrix3d pd_v_bg;
    Eigen::Matrix3d pd_v_ba;

    Eigen::Matrix3d cov_acc;
    Eigen::Matrix3d cov_gyro;
    Eigen::Matrix<double, 9, 9> cov_residual;

};