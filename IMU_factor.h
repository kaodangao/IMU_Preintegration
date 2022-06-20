#include <Eigen/Eigen>
#include <g2o/core/base_vertex.h>
#include <g2o/core/base_binary_edge.h>
#include <g2o/core/base_multi_edge.h>
#include "math_preintegration.h"

class VertexIMU : public g2o::BaseVertex <15, Vector15d>{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    virtual bool read(std::istream &in) {}

    virtual bool write(std::ostream &out) const {}

    virtual void setToOriginImpl() override {
        _estimate = Vector15d::Zero();
    }
    virtual void oplusImpl(const double* update_) override {
        Eigen::Map<const Vector15d> update(update_);

        Vector15d new_estimate;

        new_estimate.segment(0,3) = Log(exp(estimate().head(3))*exp(update.head(3)));
        new_estimate.segment(3,3) = exp(estimate().head(3))*update.segment(3,3)+_estimate.segment(3,3);
        new_estimate.tail(9) = estimate().tail(9)+ update.tail(9);

        setEstimate(new_estimate);
    }
};


class EdgeT : public g2o::BaseBinaryEdge <6, g2o::SE3Quat, VertexIMU, VertexIMU>{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    virtual bool read(std::istream &in) {}

    virtual bool write(std::ostream &out) const {}

    virtual void linearizeOplus() override {
        const VertexIMU* vi = static_cast<const VertexIMU*>(_vertices[0]);
        const VertexIMU* vj = static_cast<const VertexIMU*>(_vertices[1]);
        g2o::SE3Quat Ti(exp(vi->estimate().head(3)),vi->estimate().segment(3,3));

        g2o::SE3Quat Tj(exp(vj->estimate().head(3)),vj->estimate().segment(3,3));

        const g2o::SE3Quat& Tji = _measurement;
        g2o::SE3Quat e = Tj.inverse()*Tji*Ti;

        _jacobianOplusXi.block(0,0,3,3) = jacobian_right_inv(Log(e.rotation().matrix()));
        _jacobianOplusXi.block(3,3,3,3) = Tj.rotation().matrix().transpose()*Tji.rotation().matrix();

        _jacobianOplusXj.block(0,0,3,3) = -jacobian_right_inv(Log(e.rotation().matrix()));
        _jacobianOplusXj.block(3,0,3,3) = hat(e.translation());
        _jacobianOplusXj.block(3,3,3,3) = -Tj.rotation().matrix().transpose();
    }

    virtual void computeError() override {
        const VertexIMU* vi = static_cast<const VertexIMU*>(_vertices[0]);
        const VertexIMU* vj = static_cast<const VertexIMU*>(_vertices[1]);

        g2o::SE3Quat Ti(exp(vi->estimate().head(3)),vi->estimate().segment(3,3));
        
        g2o::SE3Quat Tj(exp(vj->estimate().head(3)),vj->estimate().segment(3,3));

        g2o::SE3Quat Tji(_measurement);


        g2o::SE3Quat error = Tj.inverse()*Tji*Ti;
        
        _error.head(3) =  Log(error.rotation().matrix());
        _error.segment(3,3) = error.translation();
    }
};

class EdgeIMU : public g2o::BaseBinaryEdge <9, Vector9d, VertexIMU, VertexIMU>{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    virtual bool read(std::istream &in) {}
    virtual bool write(std::ostream &out) const {}

    void parameter_input(   double deltatij,
                            Vector3d deltap,
                            Vector3d deltav,
                            Matrix3d deltar,
                            Matrix3d pd_r_gyro,
                            Matrix3d pd_p_gyro,
                            Matrix3d pd_p_acc,
                            Matrix3d pd_v_gyro,
                            Matrix3d pd_v_acc){
        delta_tij = deltatij;
        delta_p = deltap;
        delta_v = deltav;
        delta_r = deltar;
        pd_r_bg = pd_r_gyro;
        pd_p_bg = pd_p_gyro;
        pd_p_ba = pd_p_acc;
        pd_v_bg = pd_v_gyro;
        pd_v_ba = pd_v_acc;
        g << 0, 0, -9.8;

    }
    virtual void linearizeOplus() override {
        const VertexIMU* vi = static_cast<const VertexIMU*>(_vertices[0]);
        const VertexIMU* vj = static_cast<const VertexIMU*>(_vertices[1]);

        Matrix3d R_i = exp(vi->estimate().segment(0,3));
        Matrix3d R_j = exp(vj->estimate().segment(0,3));

        Vector3d P_i = vi->estimate().segment(3,3);
        Vector3d P_j = vj->estimate().segment(3,3);

        Vector3d V_i = vi->estimate().segment(6,3);
        Vector3d V_j = vj->estimate().segment(6,3);

        Vector3d delta_ba_i = vj->estimate().segment(9,3);
        Vector3d delta_bg_i = vj->estimate().segment(12,3);

        Matrix3d delta_r_meas = delta_r*exp(pd_r_bg*delta_bg_i);
        Vector3d residual_delta_rij = Log(delta_r_meas.transpose()*R_i.transpose()*R_j);
        Matrix3d jac_rij_inv = jacobian_right_inv(residual_delta_rij);

        //delta_r to ri 
        _jacobianOplusXi.block(0,0,3,3) = -jac_rij_inv*R_j.transpose()*R_i;
        //delta_v to ri
        _jacobianOplusXi.block(3,0,3,3) = hat(R_i.transpose()*(V_j-V_i-g*delta_tij));
        //delta_v to vi
        _jacobianOplusXi.block(3,6,3,3) = -R_i.transpose();
        //delta_p to ri
        _jacobianOplusXi.block(6,0,3,3) = hat(R_i.transpose()*(P_j-P_i-V_i*delta_tij-0.5*g*delta_tij*delta_tij));
        //delta_p to pi
        _jacobianOplusXi.block(6,3,3,3) = -Matrix3d::Identity();
        //delta_p to vi
        _jacobianOplusXi.block(6,6,3,3) = -R_i.transpose()*delta_tij;

        //delta_r to rj
        _jacobianOplusXj.block(0,0,3,3) = jac_rij_inv;
        //delta_v to vj
        _jacobianOplusXj.block(3,6,3,3) = R_i.transpose();
        //delta_p to pj
        _jacobianOplusXj.block(6,3,3,3) = R_i.transpose()*R_j; 
        //delta_r to delta_bgi
        _jacobianOplusXj.block(0,12,3,3) = -jac_rij_inv*exp(-residual_delta_rij)*jacobian_right(pd_r_bg*delta_bg_i)*pd_r_bg;
        //delta_v to delta_bai
        _jacobianOplusXj.block(3,9,3,3) = -pd_v_ba;
        //delta_v to delta_bgi
        _jacobianOplusXj.block(3,12,3,3) = -pd_v_bg;
        //delta_p to delta_bai
        _jacobianOplusXj.block(6,9,3,3) = -pd_p_ba;
        //delta_p to delta_bgi
        _jacobianOplusXj.block(6,12,3,3) = -pd_p_bg;

    }
    virtual void computeError() override {
        const VertexIMU* vi = static_cast<const VertexIMU*>(_vertices[0]);
        const VertexIMU* vj = static_cast<const VertexIMU*>(_vertices[1]);

        Matrix3d R_i = exp(vi->estimate().segment(0,3));
        Matrix3d R_j = exp(vj->estimate().segment(0,3));

        Vector3d P_i = vi->estimate().segment(3,3);
        Vector3d P_j = vj->estimate().segment(3,3);

        Vector3d V_i = vi->estimate().segment(6,3);
        Vector3d V_j = vj->estimate().segment(6,3);

        Vector3d delta_ba_i = vj->estimate().segment(9,3);
        Vector3d delta_bg_i = vj->estimate().segment(12,3);

        Matrix3d delta_r_meas = delta_r*exp(pd_r_bg*delta_bg_i);
        Vector3d delta_v_meas = delta_v+pd_v_ba*delta_ba_i+pd_v_bg*delta_bg_i;
        Vector3d delta_p_meas = delta_p+pd_p_ba*delta_ba_i+pd_p_bg*delta_bg_i;

        _error.segment(0,3) = Log(delta_r_meas.transpose()*R_i.transpose()*R_j);
        _error.segment(3,3) = R_i.transpose()*(V_j-V_i-g*delta_tij)-delta_v_meas;
        _error.segment(6,3) = R_i.transpose()*(P_j-P_i-V_i*delta_tij-0.5*g*delta_tij*delta_tij)-delta_p_meas;

    }

    private:
    
    double delta_tij;
    Vector3d g;
    Vector3d delta_p;
    Vector3d delta_v;
    Matrix3d delta_r;

    Matrix3d pd_r_bg;
    Matrix3d pd_p_bg;
    Matrix3d pd_p_ba;
    Matrix3d pd_v_bg;
    Matrix3d pd_v_ba;

};