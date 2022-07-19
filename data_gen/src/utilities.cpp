#include "utilities.h"

void LoadPose(std::string filename, std::vector<MotionData>& pose)
{

    std::ifstream f;
    f.open(filename.c_str());

    if(!f.is_open())
    {
        std::cerr << " can't open LoadFeatures file "<<std::endl;
        return;
    }

    while (!f.eof()) {

        std::string s;
        std::getline(f,s);

        if(! s.empty())
        {
            std::stringstream ss;
            ss << s;

            MotionData data;
            double time;
            Eigen::Quaterniond q;
            Eigen::Vector3d t;
            Eigen::Vector3d gyro;
            Eigen::Vector3d acc;

            ss>>time;
            ss>>q.w();
            ss>>q.x();
            ss>>q.y();
            ss>>q.z();
            ss>>t(0);
            ss>>t(1);
            ss>>t(2);
            ss>>gyro(0);
            ss>>gyro(1);
            ss>>gyro(2);
            ss>>acc(0);
            ss>>acc(1);
            ss>>acc(2);


            data.timestamp = time;
            data.imu_gyro = gyro;
            data.imu_acc = acc;
            data.twb = t;
            data.Rwb = Eigen::Matrix3d(q);
            pose.push_back(data);

        }
    }

}

void save_Pose(std::string filename, std::vector<MotionData> pose)
{
    std::ofstream save_points;
    save_points.open(filename.c_str());

    for (int i = 0; i < pose.size(); ++i) {
        MotionData data = pose[i];
        double time = data.timestamp;
        Eigen::Quaterniond q(data.Rwb);
        Eigen::Vector3d t = data.twb;
        Eigen::Vector3d gyro = data.imu_gyro;
        Eigen::Vector3d acc = data.imu_acc;
        Eigen::Vector3d v = data.imu_velocity;
        Eigen::Vector3d ba = data.imu_acc_bias;
        Eigen::Vector3d bg = data.imu_gyro_bias;
        save_points<<time<<" "
                   <<q.w()<<" "
                   <<q.x()<<" "
                   <<q.y()<<" "
                   <<q.z()<<" "
                   <<t(0)<<" "
                   <<t(1)<<" "
                   <<t(2)<<" "
                   <<gyro(0)<<" "
                   <<gyro(1)<<" "
                   <<gyro(2)<<" "
                   <<acc(0)<<" "
                   <<acc(1)<<" "
                   <<acc(2)<<" "
                   <<v(0)<<" "
                   <<v(1)<<" "
                   <<v(2)<<" "
                   <<bg(0)<<" "
                   <<bg(1)<<" "
                   <<bg(2)<<" "
                   <<ba(0)<<" "
                   <<ba(1)<<" "
                   <<ba(2)<<" "
                   <<std::endl;
    }
}