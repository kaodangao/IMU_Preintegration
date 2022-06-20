cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(preintegration)

LIST( APPEND CMAKE_MODULE_PATH /home/rixap/software/g2o-master/cmake_modules )
SET(G2O_ROOT /home/rixap/software/g2o-master)

find_package(Eigen3 REQUIRED)

find_package(G2O REQUIRED)

find_package(glog 0.6.0 REQUIRED)

include_directories(${G2O_INCLUDE_DIRS} "/usr/include/eigen3")
include_directories(${EIGEN3_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS "-std=c++14")
add_executable (preintegration test_preintegration.cpp)


set(G2O_LIBS g2o_cli g2o_ext_freeglut_minimal g2o_simulator g2o_solver_slam2d_linear g2o_types_icp g2o_types_slam2d g2o_core g2o_interface g2o_solver_csparse g2o_solver_structure_only g2o_types_sba g2o_types_slam3d g2o_csparse_extension g2o_opengl_helper g2o_solver_dense g2o_stuff g2o_types_sclam2d g2o_parser g2o_solver_pcg g2o_types_data g2o_types_sim3 cxsparse )

target_link_libraries(${PROJECT_NAME} ${G2O_LIBS} glog::glog)
