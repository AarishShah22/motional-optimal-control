#include "rule_based_optimal_control.h"

constexpr double deg_2_rad = 3.14F/180.0;

Agent::Agent(Eigen::Vector3f& pose) {
    pose_ = pose;
}

Eigen::Vector3f Agent::getPose() const {
    return pose_;
} 

EgoVehicle::EgoVehicle() {
    rear_axle_to_COG_ = 1.0F;
    length_of_vehicle_ = 2.0F;
    Eigen::Vector3f pose1{1.0F,1.0F,1.0F};
    Agent* agent1 = new Agent(pose1);
    agents_.push_back(agent1);
}

EgoVehicle::~EgoVehicle() {
    for (auto agent : agents_) {
        delete agent;
    }
    agents_.clear();
}

void EgoVehicle::setReferenceTrajectory(std::vector<Eigen::Vector3f>& ref_traj) {
    ref_traj_ = ref_traj;
}

void EgoVehicle::findOptimalTrajectory() {
    // need:
    //  1. obj fn - done
    //  2. init condition - done
    //  3. linear inequality/equality constraints - done
    //  4. input lower+upper bounds - done
    //  5. nonlinear constraints
}

OptimizationProblem::OptimizationProblem(Ipopt::Index num_timesteps, Eigen::Vector3f init_pose, double dt) : 
                    num_timesteps_(num_timesteps), init_pose_(init_pose), dt_(dt) {        
} 

bool OptimizationProblem::get_nlp_info(Ipopt::Index& n, Ipopt::Index& m,
                      Ipopt::Index& nnz_jac_g, Ipopt::Index& nnz_h_lag,
                      Ipopt::TNLP::IndexStyleEnum& index_style) {
    n = 6*num_timesteps_;
    m = 8*num_timesteps_ + 5;
    nnz_jac_g = 12*num_timesteps_;
    nnz_h_lag = 0;
    index_style = TNLP::C_STYLE;
    return true;
}

bool OptimizationProblem::get_bounds_info(Ipopt::Index n, double* x_l, double* x_u,
                         Ipopt::Index m, double* g_l, double* g_u) {
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < num_timesteps_; j++) {
            if (i == 0) {
                x_l[i*num_timesteps_ + j] = 0;
                x_u[i*num_timesteps_ + j] = 20;
            }
            else if (i == 1) {
                x_l[i*num_timesteps_ + j] = -50*deg_2_rad;
                x_u[i*num_timesteps_ + j] = 50*deg_2_rad;
            }
            else if (i == 2) {
                x_l[i*num_timesteps_ + j] = 0;
                x_u[i*num_timesteps_ + j] = std::numeric_limits<Ipopt::Number>::infinity();
            }
            else {
                // strict constraint for all CBF for now - need to modify this.
                x_l[i*num_timesteps_ + j] = 0; 
                x_u[i*num_timesteps_ + j] = 0;
            }
        }
    }
    for (int k = 0; k < m; k++) {
        g_l[k] = -std::numeric_limits<Ipopt::Number>::infinity();
        g_u[k] = 0;
    }
    return true;
}

bool OptimizationProblem::get_starting_point(Ipopt::Index n, bool init_x, double* x,
                               bool init_z, double*, double*,
                               Ipopt::Index m, bool init_lambda, double*) {
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < num_timesteps_; j++) {
            if (i == 0) {
                x[i*num_timesteps_ + j] = 10;
            }
            else if (i == 1) {
                x[i*num_timesteps_ + j] = 0.3;
            }
            else if (i == 2) {
                x[i*num_timesteps_ + j] = 10;
            }
            else {
                x[i*num_timesteps_ + j] = 0;
            }
        }
    }
    return true;
}

bool OptimizationProblem::eval_f(Ipopt::Index n, const double* x, bool, double& obj_value) {
    std::vector<double> cost_coeffs{0.1,0.1,0.01,100,1.5,0.1};
    obj_value = 0;
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < num_timesteps_; j++) {
            obj_value += cost_coeffs[i]*x[i*num_timesteps_ + j];
        }
    }
    return true;
}

bool OptimizationProblem::eval_g(Ipopt::Index n, const double* x, bool new_x, Ipopt::Index m, double* g) {
    std::vector<double> x_veh(num_timesteps_);
    std::vector<double> y_veh(num_timesteps_);
    std::vector<double> theta_veh(num_timesteps_);
    x_veh[0] = init_pose_[0];
    y_veh[0] = init_pose_[1];
    theta_veh[0] = init_pose_[2];
    double* v = x;
    double* omega = x + num_timesteps_;
    double* delta_e = x + 2*num_timesteps_;
    double* delta_1 = x + 3*num_timesteps_;
    double* delta_2 = x + 4*num_timesteps_;
    double* delta_3 = x + 5*num_timesteps_;
    for (int i = 0; i < num_timesteps_-1; i++) {
        x_veh[i+1] = x_veh[i] + v[i]*cos(theta_veh[i])*dt_;
        y_veh[i+1] = y_veh[i] + v[i]*sin(theta_veh[i])*dt_;
        theta_veh[i+1] = theta_veh[i] + omega[i]*dt_;
    }
    return true;
}

void OptimizationProblem::finalize_solution(Ipopt::SolverReturn status,
                           Ipopt::Index n, const double* x, const double* z_L,
                           const double* z_U, Ipopt::Index m,
                           const double* g, const double* lambda,
                           double obj_value,
                           const Ipopt::IpoptData* ip_data,
                           Ipopt::IpoptCalculatedQuantities* ip_cq) {

}

int main() {
    return 0;
}