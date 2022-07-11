#include <eigen3/Eigen/Eigen>
#include <fstream>
#include <iostream>
#include <g2o/core/block_solver.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/core/sparse_optimizer.h>
#include <g2o/solvers/eigen/linear_solver_eigen.h>
#include <g2o/solvers/dense/linear_solver_dense.h>

#include "math_preintegration.h"
#include "IMU_Preintegration.h"
#include "IMU_factor.h"

using namespace std;

#define FRAME_NUM 20 // The number of keyframes in one optimization
#define OPT_NUM 50 // The number of optimization interation

/**
 * read imu data                 
 * data order: timestamp(0), quaternion(1-4), translation(5-7), angular velocity(8-10), 
 *             acceleration(11-13), velocity(14-16), angular velocity bias(17-19),
 *             acceleration bias(21-23)
 * @param path txtfile path
 * @return a 2d matrix of imu data
 */
vector< vector<double> > read_data(string path){

    ifstream read_imu(path);
    string imu_data;
    vector< vector<double> > IMUDATA;
    vector<double> IMUDATA_LINE;
    stringstream ss;

    while(getline(read_imu, imu_data)){

        double data = 0;
        ss = stringstream();
        ss<<imu_data;
        int count = 0;

        while(ss >> data){

            IMUDATA_LINE.push_back(data);
        }
        IMUDATA.push_back(IMUDATA_LINE);
        IMUDATA_LINE.clear();
    }
    return IMUDATA;
}

int main(int argc, char *argv[]){

    Vector3d ba0(0,0,0);
    Vector3d bg0(0,0,0);
    IMU_Preintegration IMUP(ba0,bg0,0.019,0.015); //0.019 is noise parameter of angular velocity 0.015 is noise parameter of acceleration 
    
    vector< vector<double> > imudata;
    if(argv[1]){
        imudata = read_data(argv[1]);
    }
    else{
        imudata = read_data("../imu_pose_noise.txt");
    }

    g2o::OptimizationAlgorithmLevenberg* solver =
      new g2o::OptimizationAlgorithmLevenberg(g2o::make_unique<g2o::BlockSolverX>(
          g2o::make_unique<
              g2o::LinearSolverDense<g2o::BlockSolverX::PoseMatrixType>>()));
    g2o::SparseOptimizer optimizer;
    optimizer.setAlgorithm(solver);

    Vector3d g(0,0,-9.81);

    Eigen::Quaterniond qi;
    Vector3d ti;

    Eigen::Quaterniond qj;
    Vector3d tj;

    Vector3d velocity0(imudata[0][14],imudata[0][15],imudata[0][16]);
    Vector3d t0(imudata[0][5], imudata[0][6], imudata[0][7]);
    Eigen::Quaterniond q0(imudata[0][1], imudata[0][2], imudata[0][3], imudata[0][4]);

    qi = q0;
    ti = t0;

    double chi2_v_in = 0;
    double chi2_v_out = 0;

    double chi2_r_in = 0;
    double chi2_r_out = 0;

    double chi2_p_in = 0;
    double chi2_p_out = 0;

    double chi2_ba_in = 0;
    double chi2_ba_out = 0;

    double chi2_bg_in = 0;
    double chi2_bg_out = 0;

    Vector3d ba;
    Vector3d bg;

    std::srand(std::time(nullptr));

    std::ofstream save_vision_points;
    save_vision_points.open("../vision_pose.txt");


    for(int j = 0; j < FRAME_NUM; j++){

        Eigen::Quaterniond q1(imudata[j*20+20][1], imudata[j*20+20][2], imudata[j*20+20][3], imudata[j*20+20][4]);
        q1.normalize();
        qj = q1;

        tj << imudata[j*20+20][5], imudata[j*20+20][6], imudata[j*20+20][7];

        Vector3d vj(imudata[j*20+20][14],imudata[j*20+20][15],imudata[j*20+20][16]);

        bg << imudata[j*20+20][17],imudata[j*20+20][18],imudata[j*20+20][19];
        ba << imudata[j*20+20][20],imudata[j*20+20][21],imudata[j*20+20][22];

        for(int i = 0; i < 20; i++){

            Vector3d gyro;
            gyro << imudata[j*20+i][8], imudata[j*20+i][9], imudata[j*20+i][10];
            Vector3d acc;
            acc << imudata[j*20+i][11], imudata[j*20+i][12], imudata[j*20+i][13];

            bg << imudata[j*20+i][17],imudata[j*20+i][18],imudata[j*20+i][19];
            ba << imudata[j*20+i][20],imudata[j*20+i][21],imudata[j*20+i][22];

            chi2_bg_in += (bg0-bg).norm()*(bg0-bg).norm();
            chi2_ba_in += (ba0-ba).norm()*(ba0-ba).norm();



            IMUP.Calculate(acc, gyro, imudata[j*20+i+1][0]-imudata[j*20+i][0]);
        }
        // std::cout<<"delta_t = "<<IMUP.delta_tij_read()<<std::endl;
        // std::cout<<"imucov = "<<IMUP.cov_residual_read()<<std::endl;
        // std::cout<<"delta_p = "<<IMUP.delta_p_read()<<std::endl;
        // std::cout<<"delta_v = "<<IMUP.delta_v_read()<<std::endl;
        // std::cout<<"r/bg = "<<IMUP.pd_r_bg_read()<<std::endl;
        // std::cout<<"p/ba = "<<IMUP.pd_p_ba_read()<<std::endl;
        // std::cout<<"p/bg = "<<IMUP.pd_p_bg_read()<<std::endl;
        // std::cout<<"v/ba = "<<IMUP.pd_v_ba_read()<<std::endl;
        // std::cout<<"v/bg = "<<IMUP.pd_v_bg_read()<<std::endl;
        if( j == 0 ){
            VertexIMU* v0 = new VertexIMU();

            Vector15d state_i;
            state_i<<Log(q0.toRotationMatrix()),t0,velocity0,Vector6d::Zero();

            v0->setEstimate(state_i);
            v0->setId(j);
            v0->setFixed(true);
            optimizer.addVertex(v0); 

        }

        VertexIMU* v1 = new VertexIMU();

        double error = Vector3d::Random().normalized()(0)/100+1;
        Eigen::AngleAxisd r_error(0.01,Vector3d::Random().normalized());

        Vector15d state_j;
        std::cout<<"error:"<<error<<std::endl;
        state_j<<Log(r_error.toRotationMatrix()*qj.toRotationMatrix()), tj*error, vj, Vector6d::Zero();

        //state_j<<Vector15d::Zero();

        chi2_v_in += (state_j.segment(6,3)-vj).norm()*(state_j.segment(6,3)-vj).norm();
        chi2_r_in += (state_j.head(3)-Log(qj.toRotationMatrix())).norm()*(state_j.head(3)-Log(qj.toRotationMatrix())).norm();
        chi2_p_in += (state_j.segment(3,3)-tj).norm()*(state_j.segment(3,3)-tj).norm();

        qj = exp(state_j.head(3));
        tj = state_j.segment(3,3);
        

        save_vision_points
            <<qj.w()<<" "
            <<qj.x()<<" "
            <<qj.y()<<" "
            <<qj.z()<<" "
            <<state_j(3)<<" "
            <<state_j(4)<<" "
            <<state_j(5)<<" "
            <<std::endl;
        // std::cout<<"state_"<<j+1<<":"<<std::endl;
        // std::cout<<Log(qj.toRotationMatrix())<<std::endl;
        // std::cout<<tj<<std::endl;
        // std::cout<<imudata[j*20+20][14]<<std::endl;
        // std::cout<<imudata[j*20+20][15]<<std::endl;
        // std::cout<<imudata[j*20+20][16]<<std::endl;

        //velocity_i += g*IMUP.delta_tij_read()+qi.toRotationMatrix()*IMUP.delta_v_read();

        v1->setEstimate(state_j);
        v1->setId(j+1);
        v1->setFixed(false);
        optimizer.addVertex(v1);

        EdgeIMU* e = new EdgeIMU();

        Vector9d meas;
        meas << Log(IMUP.delta_r_read()),IMUP.delta_v_read(),IMUP.delta_p_read();
        
        e->setVertex(0,dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(j)));
        e->setVertex(1,dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(j+1)));
        e->setMeasurement(meas);
        e->setInformation(IMUP.cov_residual_read().inverse());
        e->parameter_input(IMUP.delta_tij_read(),IMUP.delta_p_read(),IMUP.delta_v_read(),IMUP.delta_r_read(),IMUP.pd_r_bg_read(),
                        IMUP.pd_p_bg_read(),IMUP.pd_p_ba_read(),IMUP.pd_v_bg_read(),IMUP.pd_v_ba_read());
        optimizer.addEdge(e);

        EdgeT* et = new EdgeT();

        Vector6d tji_meas;
        tji_meas<<Log(qj.toRotationMatrix()*qi.toRotationMatrix().transpose()),tj-qj.toRotationMatrix()*qi.toRotationMatrix().transpose()*ti;
        et->setVertex(0,dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(j)));
        et->setVertex(1,dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(j+1)));
        et->setMeasurement(tji_meas);
        et->setInformation(Matrix6d::Identity()*1000);
        optimizer.addEdge(et);

        IMUP.reset(ba0,bg0);
        qi = qj;
        ti << tj;
    }



    optimizer.initializeOptimization();
    optimizer.setVerbose( true );
    optimizer.optimize(OPT_NUM);

    VertexIMU* v_;

    std::ofstream save_points;
    save_points.open("../pre_imu_pose.txt");

    Eigen::Quaterniond Qwb(imudata[0][1], imudata[0][2], imudata[0][3], imudata[0][4]);
    Vector3d Pwb(imudata[0][5], imudata[0][6], imudata[0][7]);
    Vector3d Vw(imudata[0][14],imudata[0][15],imudata[0][16]);
    Vector3d acc(0,0,0);
    Vector3d gyro(0,0,0);
    std::cout<<"qwb:"<<Qwb.toRotationMatrix()<<std::endl;


    for(int i = 0; i < FRAME_NUM + 1; i++){

        Vector3d vj(imudata[i*20][14],imudata[i*20][15],imudata[i*20][16]);
        Eigen::Quaterniond qj(imudata[i*20][1], imudata[i*20][2], imudata[i*20][3], imudata[i*20][4]);

        tj<<imudata[i*20][5], imudata[i*20][6], imudata[i*20][7];


        
        v_ = static_cast<VertexIMU*>(optimizer.vertex(i));
        //std::cout<<"v"<<i<<": "<<std::endl<<v_->estimate()<<std::endl;
        // std::cout<<"r:"<<qj.toRotationMatrix()<<std::endl;
        // std::cout<<"r:"<<Log(qj.toRotationMatrix())<<std::endl;
        Vector3d pw = v_->estimate().segment(3,3);
        Vector3d vw = v_->estimate().segment(6,3);

        chi2_v_out += (vw-vj).norm()*(vw-vj).norm();
        chi2_r_out += (v_->estimate().segment(0,3)-Log(qj.toRotationMatrix())).norm()*(v_->estimate().segment(0,3)-Log(qj.toRotationMatrix())).norm();
        chi2_p_out += (pw-tj).norm()*(pw-tj).norm();
        chi2_bg_out += (bg-v_->estimate().segment(12,3)).norm()*(bg-v_->estimate().segment(12,3)).norm();
        chi2_ba_out += (ba-v_->estimate().segment(9,3)).norm()*(ba-v_->estimate().segment(9,3)).norm();


        if(i==FRAME_NUM){
            continue;
        }
        else{
            v_ = static_cast<VertexIMU*>(optimizer.vertex(i+1));
            for(int j = 0; j < 20; j++){
                
                bg<<imudata[i*20+j][17],imudata[i*20+j][18],imudata[i*20+j][19];

                ba<<imudata[i*20+j][20],imudata[i*20+j][21],imudata[i*20+j][22];   

                chi2_bg_out += (bg-v_->estimate().segment(12,3)-bg0).norm()*(bg-v_->estimate().segment(12,3)-bg0).norm();
                chi2_ba_out += (ba-v_->estimate().segment(9,3)-ba0).norm()*(ba-v_->estimate().segment(9,3)-ba0).norm(); 
            
                double dt = 0.005;
                gyro << imudata[i*20+j+1][8],imudata[i*20+j+1][9],imudata[i*20+j+1][10];
                acc << imudata[i*20+j+1][11],imudata[i*20+j+1][12],imudata[i*20+j+1][13];
                //std::cout<<"gyro:"<<gyro.transpose()<<std::endl;
                bg = v_->estimate().segment(12,3);
                ba = v_->estimate().segment(9,3);
                Vector3d v_p = v_->estimate().segment(6,3);

                // ba <<0,0,0;
                // bg <<0,0,0;

                // Pwb += Vw*dt+0.5*g*dt*dt+0.5*Qwb.toRotationMatrix()*(acc+ba)*dt*dt;
                // Vw += Qwb*(acc+ba)*dt+g*dt;
                // Qwb = Qwb.toRotationMatrix()*exp((gyro+bg)*dt);
                Eigen::Quaterniond dq;
                Eigen::Vector3d dtheta_half =  (gyro-bg) * dt /2.0;
                dq.w() = 1;
                dq.x() = dtheta_half.x();
                dq.y() = dtheta_half.y();
                dq.z() = dtheta_half.z();
                dq.normalize();
                
                /// imu 动力学模型 欧拉积分
                Eigen::Vector3d acc_w = Qwb * (acc-ba) + g;  // aw = Rwb * ( acc_body - acc_bias ) + gw
                Qwb = Qwb * dq;
                Pwb = Pwb + Vw * dt + 0.5 * dt * dt * acc_w;
                Vw = Vw + acc_w * dt;

                save_points
                    <<Qwb.w()<<" "
                    <<Qwb.x()<<" "
                    <<Qwb.y()<<" "
                    <<Qwb.z()<<" "
                    <<Pwb(0)<<" "
                    <<Pwb(1)<<" "
                    <<Pwb(2)<<" "
                    <<v_p(0)<<" "
                    <<v_p(1)<<" "
                    <<v_p(2)<<" "
                    <<bg(0)<<" "
                    <<bg(1)<<" "
                    <<bg(2)<<" "
                    <<ba(0)<<" "
                    <<ba(1)<<" "
                    <<ba(2)<<" "
                    <<imudata[i*20+j+1][0]<<" "
                    <<std::endl;

            }        
        }


    }
    std::cout<<"chi2_v_in : "<<chi2_v_in<<std::endl;
    std::cout<<"chi2_v_out : "<<chi2_v_out<<std::endl;
    std::cout<<"chi2_r_in : "<<chi2_r_in<<std::endl;
    std::cout<<"chi2_r_out : "<<chi2_r_out<<std::endl;
    std::cout<<"chi2_p_in : "<<chi2_p_in<<std::endl;
    std::cout<<"chi2_p_out : "<<chi2_p_out<<std::endl;
    std::cout<<"chi2_bg_in : "<<chi2_bg_in<<std::endl;
    std::cout<<"chi2_bg_out : "<<chi2_bg_out<<std::endl;
    std::cout<<"chi2_ba_in : "<<chi2_ba_in<<std::endl;
    std::cout<<"chi2_ba_out : "<<chi2_ba_out<<std::endl;

    return 0;
}