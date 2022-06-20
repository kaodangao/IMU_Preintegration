#include <Eigen/Eigen>
#include <fstream>
#include <iostream>
#include <g2o/core/block_solver.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/core/sparse_optimizer.h>
#include <g2o/solvers/eigen/linear_solver_eigen.h>
#include <g2o/solvers/dense/linear_solver_dense.h>
#include <g2o/core/base_vertex.h>
#include <g2o/core/base_binary_edge.h>
#include <g2o/core/base_multi_edge.h>
#include <g2o/types/slam3d/se3quat.h>

#include "math_preintegration.h"
#include "IMU_Preintegration.h"
#include "IMU_factor.h"

using namespace std;

vector< vector<double> > read_data(string path){
    ifstream read_imu(path);
    int flag = 0;
    string imu_data;
    vector< vector<double> > IMUDATA;
    vector<double> IMUDATA_LINE;
    stringstream ss;

    while(flag<21){
        flag++;
        getline(read_imu, imu_data);
        double data = 0;
        ss = stringstream();
        ss<<imu_data;
        int count = 0;
        while(count<14){
            ss >> data;
            if(count == 0){
                IMUDATA_LINE.push_back(data);

            }
            else if(count>7){
                IMUDATA_LINE.push_back(data);
            }
            count ++;

        }
        IMUDATA.push_back(IMUDATA_LINE);
        IMUDATA_LINE.clear();

    }
    return IMUDATA;
}

int main(int argc, char *argv[]){
    Vector3d ba(0,0,0);
    Vector3d bg(0,0,0);
    IMU_Preintegration IMUP(ba,bg,0.019,0.015);
    vector< vector<double> > imudata;
    imudata = read_data(argv[1]);
    for(int i = 0; i < imudata.size()-1; i++){
        Eigen::Vector3d gyro;
        gyro << imudata[i][1], imudata[i][2], imudata[i][3];
        Eigen::Vector3d acc;
        acc << imudata[i][4], imudata[i][5], imudata[i][6];

        IMUP.Calculate(acc, gyro, imudata[i+1][0]-imudata[i][0]);
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


    typedef g2o::BlockSolver<g2o::BlockSolverTraits<Eigen::Dynamic,Eigen::Dynamic>> Block;
    std::unique_ptr<Block::LinearSolverType> linearSolver (new g2o::LinearSolverDense<Block::PoseMatrixType>());
    std::unique_ptr<Block> solver_ptr (new Block(std::move(linearSolver)));
    g2o::OptimizationAlgorithmLevenberg* solver =   new g2o::OptimizationAlgorithmLevenberg(std::move(solver_ptr));

    g2o::SparseOptimizer optimizer;
    optimizer.setAlgorithm(solver);

    Eigen::Quaterniond q1(0.99875 ,0.0499792, 0, 0);
    Eigen::Vector3d p1(20 ,5 ,5);

    g2o::SE3Quat t1_(q1,p1);

    Eigen::Quaterniond q2(0.998598 ,0.0495645 ,0.0107507 ,0.0151907);
    Eigen::Vector3d p2(19.9926, 5.62822, 5.30902);

    g2o::SE3Quat t2_(q2,p2);
    std::cout<<"t2_"<<std::endl<<t2_<<std::endl;
    Eigen::Vector3d vel1(-0, 6.28319 ,3.14159);
    VertexIMU* v0 = new VertexIMU();
    Vector15d state_i = Vector15d::Zero();

    state_i<<Log(q1.matrix()),p1,vel1,Vector6d::Zero();
    v0->setEstimate(state_i);
    v0->setId(0);
    v0->setFixed(true);
    optimizer.addVertex(v0);

    VertexIMU* v1 = new VertexIMU();
    Vector15d state_j = Vector15d::Zero();
    Eigen::Vector3d k(0.1,0.1,0.1);
    state_j<<Log(q2.matrix())+k/10,p2+k,vel1+IMUP.delta_v_read(),Vector6d::Zero();
    v1->setEstimate(state_j);
    v1->setId(1);
    v1->setFixed(false);
    optimizer.addVertex(v1);

    EdgeIMU* e = new EdgeIMU();
    Vector9d meas;
    meas << Log(IMUP.delta_r_read()),IMUP.delta_v_read(),IMUP.delta_p_read();
    e->setVertex(0,dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(0)));
    e->setVertex(1,dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(1)));
    e->setMeasurement(meas);
    e->setInformation(IMUP.cov_residual_read().inverse());
    e->parameter_input(IMUP.delta_tij_read(),IMUP.delta_p_read(),IMUP.delta_v_read(),IMUP.delta_r_read(),IMUP.pd_r_bg_read(),
                       IMUP.pd_p_bg_read(),IMUP.pd_p_ba_read(),IMUP.pd_v_bg_read(),IMUP.pd_v_ba_read());

    optimizer.addEdge(e);

    optimizer.setVerbose( true );

    EdgeT* et = new EdgeT();
    g2o::SE3Quat t21 = t2_ * t1_.inverse();


    et->setVertex(0,dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(0)));
    et->setVertex(1,dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(1)));
    et->setMeasurement(t21);
    et->setInformation(Matrix6d::Identity());
    optimizer.addEdge(et);
    std::cout<<"nb"<<std::endl;
    optimizer.initializeOptimization();

    optimizer.optimize(30);

    VertexIMU* v0_ = static_cast<VertexIMU*>(optimizer.vertex(0));
    VertexIMU* v1_ = static_cast<VertexIMU*>(optimizer.vertex(1));
    Vector15d b = v1_->estimate();
    std::cout<<"vj:"<<b<<std::endl;
    return 0;
}