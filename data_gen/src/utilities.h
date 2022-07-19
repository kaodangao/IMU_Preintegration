#ifndef IMUSIMWITHPOINTLINE_UTILITIES_H
#define IMUSIMWITHPOINTLINE_UTILITIES_H

#include "imu.h"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <iostream>
#include <vector>
#include <fstream>

void LoadPose(std::string filename, std::vector<MotionData>& pose);
void save_Pose(std::string filename, std::vector<MotionData> pose);

#endif //IMUSIMWITHPOINTLINE_UTILITIES_H


