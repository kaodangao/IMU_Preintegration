#ifndef IMUSIM_PARAM_H
#define IMUSIM_PARAM_H

#include <eigen3/Eigen/Core>

class Param{

public:

    Param(){};

    // time
    int imu_frequency = 200;
    double imu_timestep = 1./imu_frequency;
    double t_start = 0.;
    double t_end = 2.2;  //  20 s

    // noise
    double gyro_bias_sigma = 1.0e-5;
    double acc_bias_sigma = 0.0001;

    double gyro_noise_sigma = 0.015;    // rad/s * 1/sqrt(hz)
    double acc_noise_sigma = 0.019;      //　m/(s^2) * 1/sqrt(hz)

};


#endif //IMUSIM_PARAM_H
