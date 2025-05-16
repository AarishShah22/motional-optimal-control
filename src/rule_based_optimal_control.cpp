#include "rule_based_optimal_control.h"

constexpr double deg_2_rad = 3.14F/180.0;

Agent::Agent(Eigen::Vector3f& pose) {
    pose_ = pose;
}

Eigen::Vector3f Agent::getPose() const {
    return pose_;
} 

EgoVehicle::EgoVehicle(Eigen::Vector3f init_pose, std::vector<Eigen::Vector3f>& ref_traj, 
                        std::vector<Agent*>& agents, double dt) : 
                        init_pose_(init_pose), ref_traj_(ref_traj), 
                        agents_(agents), dt_(dt) {   
    rear_axle_to_COG_ = 1.0F;
    length_of_vehicle_ = 2.0F;
    num_timesteps_ = ref_traj_.size();
    // priority_structure_ = 
    eps_ = 1.5;
    gamma_1_ = 0.1;
    gamma_2_ = 0.1;
    gamma_3_ = 0.1;
    r_ = 1;
    cost_coeffs_ = std::vector<double>{0.1,0.1,0.01,100,1.5,0.1};
    x_veh_ = std::vector<double>(num_timesteps_);
    y_veh_ = std::vector<double>(num_timesteps_);
    theta_veh_ = std::vector<double>(num_timesteps_);
    error_ = std::vector<double>(num_timesteps_);
    ped_pose_ = agents_[0]->getPose();
} 

EgoVehicle::~EgoVehicle() {
    for (auto agent : agents_) {
        delete agent;
    }
    agents_.clear();
}

void EgoVehicle::computeTrajectory(const double* x) {
    const double* v = x;
    const double* omega = x + num_timesteps_;
    x_veh_[0] = init_pose_[0];
    y_veh_[0] = init_pose_[1];
    theta_veh_[0] = init_pose_[2];
    for (int i = 0; i < num_timesteps_-1; i++) {
        error_[i] = (x_veh_[i] - ref_traj_[i][0])*(x_veh_[i] - ref_traj_[i][0]) +
                   (y_veh_[i] - ref_traj_[i][1])*(y_veh_[i] - ref_traj_[i][1]) +
                   (theta_veh_[i] - ref_traj_[i][2])*(theta_veh_[i] - ref_traj_[i][2]);
        x_veh_[i+1] = x_veh_[i] + v[i]*cos(theta_veh_[i])*dt_;
        y_veh_[i+1] = y_veh_[i] + v[i]*sin(theta_veh_[i])*dt_;
        theta_veh_[i+1] = theta_veh_[i] + omega[i]*dt_;
    }
}

bool EgoVehicle::feasibleSolutionFound() {
    return solution_found_;
}

bool EgoVehicle::get_nlp_info(Ipopt::Index& n, Ipopt::Index& m,
                      Ipopt::Index& nnz_jac_g, Ipopt::Index& nnz_h_lag,
                      Ipopt::TNLP::IndexStyleEnum& index_style) {
    n = 6*num_timesteps_;
    m = 5*num_timesteps_ + 5;
    nnz_jac_g = 9*num_timesteps_;
    nnz_h_lag = 0;
    index_style = TNLP::C_STYLE;
    return true;
}

bool EgoVehicle::get_bounds_info(Ipopt::Index n, double* x_l, double* x_u,
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
        if (k < 5*num_timesteps_) {
            g_l[k] = -std::numeric_limits<Ipopt::Number>::infinity();
        }
        else {
            g_l[k] = 0;
        }
        g_u[k] = 0;
    }
    return true;
}

bool EgoVehicle::get_starting_point(Ipopt::Index n, bool init_x, double* x,
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

bool EgoVehicle::eval_f(Ipopt::Index n, const double* x, bool, double& obj_value) {
    obj_value = 0;
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < num_timesteps_; j++) {
            obj_value += cost_coeffs_[i]*x[i*num_timesteps_ + j];
        }
    }
    return true;
}

bool EgoVehicle::eval_g(Ipopt::Index n, const double* x, bool new_x, Ipopt::Index m, double* g) {
    
    computeTrajectory(x);
    const double* v = x;
    const double* omega = x + num_timesteps_;
    const double* delta_e = x + 2*num_timesteps_;
    const double* delta_1 = x + 3*num_timesteps_;
    const double* delta_2 = x + 4*num_timesteps_;
    const double* delta_3 = x + 5*num_timesteps_;    

    for (int j = 0; j < num_timesteps_; j++) {
        g[j] = 2*v[j]*((x_veh_[j] - ref_traj_[j][0])*cos(theta_veh_[j]) +
                (y_veh_[j] - ref_traj_[j][1])*sin(theta_veh_[j])) + 
                2*omega[j]*(theta_veh_[j] - ref_traj_[j][2]) + eps_*error_[j] - delta_e[j];
        g[num_timesteps_ + j] = -(x_veh_[j] - ped_pose_[0])*(x_veh_[j] - ped_pose_[0])
                                - (y_veh_[j] - ped_pose_[1])*(y_veh_[j] - ped_pose_[1])
                                + r_*r_ - gamma_1_*(2*(x_veh_[j] - ped_pose_[0])*v[j]*cos(theta_veh_[j]) + 
                                2*(y_veh_[j] - ped_pose_[1])*v[j]*sin(theta_veh_[j])) + delta_1[j];
        g[2*num_timesteps_ + j] = v[j]*cos(theta_veh_[j]) - gamma_2_*(2 - y_veh_[j]) + delta_2[j];
        g[3*num_timesteps_ + j] = -gamma_3_*(v[j] - 5) + delta_3[j];
        g[4*num_timesteps_ + j] = 0.75 - y_veh_[j];
    }
    int index = 5*num_timesteps_;
    Eigen::Vector3f end_pose = ref_traj_[num_timesteps_ - 1];
    g[index] = x_veh_[num_timesteps_ - 1] - end_pose[0];
    g[index + 1] = y_veh_[num_timesteps_ - 1] - end_pose[1];
    g[index + 2] = theta_veh_[num_timesteps_ - 1] - end_pose[2];
    g[index + 3] = v[num_timesteps_ - 1];
    g[index + 4] = omega[num_timesteps_ - 1];
    return true;
}

bool EgoVehicle::eval_grad_f(Ipopt::Index n, const double* x, bool new_x, double* grad_f) {
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < num_timesteps_; j++) {
            grad_f[i*num_timesteps_ + j] = 2*cost_coeffs_[i]*x[i*num_timesteps_ + j];
        }
    }
    return true;
}

bool EgoVehicle::eval_jac_g(Ipopt::Index n, const double* x, bool,
                       Ipopt::Index m, Ipopt::Index nele_jac,
                       Ipopt::Index* iRow, Ipopt::Index* jCol,
                       double* values) {
    if (values == nullptr) {
        // structure
        for (int i = 0; i < num_timesteps_; i++) {
            iRow[i] = i; jCol[i] = i;
            iRow[num_timesteps_ + i] = i; jCol[num_timesteps_ + i] = num_timesteps_ + i;
            iRow[2*num_timesteps_ + i] = i; jCol[2*num_timesteps_ + i] = 2*num_timesteps_ + i;
            iRow[3*num_timesteps_ + i] = num_timesteps_ + i; jCol[3*num_timesteps_ + i] = i;
            iRow[4*num_timesteps_ + i] = num_timesteps_ + i; jCol[4*num_timesteps_ + i] = 3*num_timesteps_ + i;
            iRow[5*num_timesteps_ + i] = 2*num_timesteps_ + i; jCol[5*num_timesteps_ + i] = i;
            iRow[6*num_timesteps_ + i] = 2*num_timesteps_ + i; jCol[6*num_timesteps_ + i] = 4*num_timesteps_ + i;
            iRow[7*num_timesteps_ + i] = 3*num_timesteps_ + i; jCol[7*num_timesteps_ + i] = i;
            iRow[8*num_timesteps_ + i] = 3*num_timesteps_ + i; jCol[8*num_timesteps_ + i] = 5*num_timesteps_ + i;
        }
    } else {
        // values
        computeTrajectory(x);
        for (int i = 0; i < num_timesteps_; i++) {
            values[i] = 2*(x_veh_[i] - ref_traj_[i][0])*cos(theta_veh_[i]) +
                        2*(y_veh_[i] - ref_traj_[i][1])*sin(theta_veh_[i]);
            values[num_timesteps_ + i] = 2*(theta_veh_[i] - ref_traj_[i][2]);
            values[2*num_timesteps_ + i] = -1;
            values[3*num_timesteps_ + i] = -gamma_1_*(2*(x_veh_[i] - ped_pose_[0])*cos(theta_veh_[i])
                                            + 2*(y_veh_[i] - ped_pose_[1])*sin(theta_veh_[i]));
            values[4*num_timesteps_ + i] = 1;
            values[5*num_timesteps_ + i] = cos(theta_veh_[i]);
            values[6*num_timesteps_ + i] = 1;
            values[7*num_timesteps_ + i] = -gamma_3_;
            values[8*num_timesteps_ + i] = 1;
        }
    }
    return true;
}

void EgoVehicle::finalize_solution(Ipopt::SolverReturn status,
                           Ipopt::Index n, const double* x, const double* z_L,
                           const double* z_U, Ipopt::Index m,
                           const double* g, const double* lambda,
                           double obj_value,
                           const Ipopt::IpoptData* ip_data,
                           Ipopt::IpoptCalculatedQuantities* ip_cq) {
    if (status == Ipopt::SolverReturn::SUCCESS) {
        solution_found_ = true;
    } else {
        solution_found_ = false;
        std::cerr << status << std::endl;
    }
}

int main() {
    std::cout << "nothing here yet" << std::endl;
    return 0;
}