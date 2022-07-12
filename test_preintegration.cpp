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
#define OPT_NUM 100 // The number of optimization interation

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
    Vector3d vi;

    Eigen::Quaterniond qj;
    Vector3d tj;
    Vector3d vj;

    std::ofstream save_vision_points;
    save_vision_points.open("../vision_pose.txt");

    for(int j = 0; j < FRAME_NUM; j++){

        Eigen::Quaterniond q0(imudata[j*20][1], imudata[j*20][2], imudata[j*20][3], imudata[j*20][4]);
        qi = q0;
        ti << imudata[j*20][5], imudata[j*20][6], imudata[j*20][7];
        vi << imudata[j*20][14], imudata[j*20][15], imudata[j*20][16];

        Eigen::Quaterniond q1(imudata[j*20+20][1], imudata[j*20+20][2], imudata[j*20+20][3], imudata[j*20+20][4]);
        qj = q1;
        tj << imudata[j*20+20][5], imudata[j*20+20][6], imudata[j*20+20][7];
        vj << imudata[j*20+20][14], imudata[j*20+20][15], imudata[j*20+20][16];

        for(int i = 0; i < 20; i++){

            Vector3d gyro;
            gyro << imudata[j*20+i][8], imudata[j*20+i][9], imudata[j*20+i][10];
            Vector3d acc;
            acc << imudata[j*20+i][11], imudata[j*20+i][12], imudata[j*20+i][13];

            IMUP.Calculate(acc, gyro, imudata[j*20+i+1][0] - imudata[j*20+i][0]);
        }

        if( j == 0 ){

            VertexIMU* v0 = new VertexIMU();

            Vector15d state_i;
            state_i << Log(qi.toRotationMatrix()), ti, vi, Vector6d::Zero();

            v0->setEstimate(state_i);
            v0->setId(j);
            v0->setFixed(true);
            optimizer.addVertex(v0); 
        }

        VertexIMU* v1 = new VertexIMU();

        double error_coefficient = Vector3d::Random().normalized()(0)/100+1;
        Eigen::AngleAxisd r_error(0.01,Vector3d::Random().normalized());

        Vector15d state_j;
        //state_j << Log(r_error.toRotationMatrix()*qj.toRotationMatrix()), tj*error_coefficient, vj, Vector6d::Zero();
        state_j << Log(r_error.toRotationMatrix()*qj.toRotationMatrix()), tj*error_coefficient, Vector9d::Zero();

        Eigen::Quaterniond qv(Exp(state_j.head(3)));
        qj = qv;
        tj = state_j.segment(3,3);

        save_vision_points
            <<qv.w()<<" "
            <<qv.x()<<" "
            <<qv.y()<<" "
            <<qv.z()<<" "
            <<state_j(3)<<" "
            <<state_j(4)<<" "
            <<state_j(5)<<" "
            <<std::endl;

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
    }

    optimizer.initializeOptimization();
    optimizer.setVerbose( true );
    optimizer.optimize(OPT_NUM);

    VertexIMU* v_;

    std::ofstream save_points_pre;
    save_points_pre.open("../pre_imu_pose.txt");

    std::ofstream save_points_dir;
    save_points_dir.open("../dir_imu_pose.txt");

    Eigen::Quaterniond Qwb(imudata[0][1], imudata[0][2], imudata[0][3], imudata[0][4]);
    Vector3d Pwb(imudata[0][5], imudata[0][6], imudata[0][7]);
    Vector3d Vw(imudata[0][14],imudata[0][15],imudata[0][16]);
    Vector3d acc(0,0,0);
    Vector3d gyro(0,0,0);
    Vector3d ba;
    Vector3d bg;

    for(int i = 0; i < imudata.size()-1; i++){

        double dt = imudata[i+1][0]-imudata[i][0];
        gyro << imudata[i+1][8],imudata[i+1][9],imudata[i+1][10];
        acc << imudata[i+1][11],imudata[i+1][12],imudata[i+1][13];

        Eigen::Quaterniond dq;

        dq = Exp(gyro*dt);
        dq.normalize();
    
        Pwb += Vw*dt+0.5*g*dt*dt+0.5*dt*dt*(Qwb*acc);
        Vw += Qwb*acc*dt+g*dt;
        Qwb = Qwb*dq;

        save_points_dir
            <<Qwb.w()<<" "
            <<Qwb.x()<<" "
            <<Qwb.y()<<" "
            <<Qwb.z()<<" "
            <<Pwb(0)<<" "
            <<Pwb(1)<<" "
            <<Pwb(2)<<" "
            <<std::endl;
    }
    
    Eigen::Quaterniond Q0(imudata[0][1], imudata[0][2], imudata[0][3], imudata[0][4]);
    Qwb = Q0;
    Pwb << imudata[0][5], imudata[0][6], imudata[0][7];
    Vw << imudata[0][14],imudata[0][15],imudata[0][16];


    for(int i = 1; i < FRAME_NUM + 1; i++){

        v_ = static_cast<VertexIMU*>(optimizer.vertex(i));
        
        for(int j = 0; j < 20; j++){
            
            double dt = imudata[i*20+j-19][0]-imudata[i*20+j-20][0];
            gyro << imudata[i*20+j-19][8],imudata[i*20+j-19][9],imudata[i*20+j-19][10];
            acc << imudata[i*20+j-19][11],imudata[i*20+j-19][12],imudata[i*20+j-19][13];

            bg = v_->estimate().segment(12,3);
            ba = v_->estimate().segment(9,3);
            Vector3d v_p = v_->estimate().segment(6,3);

            Eigen::Quaterniond dq;
            dq = Exp((gyro-bg)*dt);
            dq.normalize();

            Pwb += Vw*dt+0.5*g*dt*dt+Qwb*(acc-ba)*0.5*dt*dt;
            Vw += Qwb*(acc-ba)*dt+g*dt;
            Qwb = Qwb*dq;

            save_points_pre
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
                <<std::endl;
            }        
    }

    return 0;
}