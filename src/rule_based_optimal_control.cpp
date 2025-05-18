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
    eps_ = 2;
    gamma_1_ = 0.1;
    gamma_2_ = 1;
    gamma_3_ = 1;
    r_ = 1;
    cost_coeffs_ = std::vector<double>{0.001,0,10,0,0,0};
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
    n = 6*num_timesteps_ + 3;
    // std::cout << num_timesteps_ << ", " << n << std::endl;
    m = 1*num_timesteps_ + 6;
    nnz_jac_g = 3*num_timesteps_ + 6;
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
                x_l[i*num_timesteps_ + j] = -std::numeric_limits<Ipopt::Number>::infinity(); 
                x_u[i*num_timesteps_ + j] = std::numeric_limits<Ipopt::Number>::infinity();
            }
        }
    }
    for (int l = 0; l < 3; l++) {
        x_l[6*num_timesteps_ + l] = 0;
        x_u[6*num_timesteps_ + l] = 0.4;
    }
    for (int k = 0; k < m; k++) {
        if (k < 1*num_timesteps_) {
            g_l[k] = -std::numeric_limits<Ipopt::Number>::infinity();
            g_u[k] = 0;    
        }
        else {
            g_l[k] = -std::numeric_limits<Ipopt::Number>::infinity();
            g_u[k] = 0;
        }
       
    }
    return true;
}

bool EgoVehicle::get_starting_point(Ipopt::Index n, bool init_x, double* x,
                               bool init_z, double*, double*,
                               Ipopt::Index m, bool init_lambda, double*) {
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < num_timesteps_; j++) {
            if (i == 0) {
                x[i*num_timesteps_ + j] = 1;
            }
            else if (i == 1) {
                x[i*num_timesteps_ + j] = 0;
            }
            else if (i == 2) {
                x[i*num_timesteps_ + j] = 0;
            }
            else {
                x[i*num_timesteps_ + j] = 0;
            }
        }
    }
    x[6*num_timesteps_] = 0.3;
    x[6*num_timesteps_ + 1] = 0;
    x[6*num_timesteps_ + 2] = 0;
    return true;
}

bool EgoVehicle::eval_f(Ipopt::Index n, const double* x, bool, double& obj_value) {
    obj_value = 0;
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < num_timesteps_; j++) {
            obj_value += cost_coeffs_[i]*x[i*num_timesteps_ + j]*x[i*num_timesteps_ + j];
        }
    }
    // computeTrajectory(x);
    // int T = num_timesteps_;
    obj_value += 1e4*(x[6*num_timesteps_]*x[6*num_timesteps_] + 
                x[6*num_timesteps_ + 1]*x[6*num_timesteps_ + 1] +
                x[6*num_timesteps_ + 2]*x[6*num_timesteps_ + 2]);
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
        // g[num_timesteps_ + j] = -gamma_1_*((x_veh_[j] - ped_pose_[0])*(x_veh_[j] - ped_pose_[0])
        //                         + (y_veh_[j] - ped_pose_[1])*(y_veh_[j] - ped_pose_[1])
        //                         - r_*r_) - 2*v[j]*((x_veh_[j] - ped_pose_[0])*cos(theta_veh_[j]) + 
        //                         (y_veh_[j] - ped_pose_[1])*sin(theta_veh_[j])) + delta_1[j];
        // g[2*num_timesteps_ + j] = v[j]*cos(theta_veh_[j]) - gamma_2_*(2 - y_veh_[j]) + delta_2[j];
        // g[3*num_timesteps_ + j] = -gamma_3_*(v[j] - 10) + delta_3[j];
        // g[4*num_timesteps_ + j] = 0 - y_veh_[j];
    }
    int index = 1*num_timesteps_;
    Eigen::Vector3f end_pose = ref_traj_[num_timesteps_ - 1];
    // g[index] = x_veh_[num_timesteps_ - 1] - end_pose[0];
    // g[index + 1] = y_veh_[num_timesteps_ - 1] - end_pose[1];
    // g[index + 2] = theta_veh_[num_timesteps_ - 1] - end_pose[2];
    // g[index + 3] = v[num_timesteps_ - 1];
    // g[index + 4] = omega[num_timesteps_ - 1];
    g[index] = x_veh_[num_timesteps_ - 1] - end_pose[0] - x[6*num_timesteps_];
    g[index + 1] = end_pose[0] - x_veh_[num_timesteps_ - 1] - x[6*num_timesteps_];
    g[index + 2] = y_veh_[num_timesteps_ - 1] - end_pose[1] - x[6*num_timesteps_ + 1];
    g[index + 3] = end_pose[1] - y_veh_[num_timesteps_ - 1] - x[6*num_timesteps_ + 1];
    g[index + 4] = theta_veh_[num_timesteps_ - 1] - end_pose[2] - x[6*num_timesteps_ + 2];
    g[index + 5] = end_pose[2] - theta_veh_[num_timesteps_ - 1] - x[6*num_timesteps_ + 2];
    return true;
}

bool EgoVehicle::eval_grad_f(Ipopt::Index n, const double* x, bool new_x, double* grad_f) {
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < num_timesteps_; j++) {
            grad_f[i*num_timesteps_ + j] = 2*cost_coeffs_[i]*x[i*num_timesteps_ + j];
        }
    }
    grad_f[6*num_timesteps_] = 2*1e4*x[6*num_timesteps_];
    grad_f[6*num_timesteps_ + 1] = 2*1e4*x[6*num_timesteps_ + 1];
    grad_f[6*num_timesteps_ + 2] = 2*1e4*x[6*num_timesteps_ + 2];
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
            // iRow[3*num_timesteps_ + i] = num_timesteps_ + i; jCol[3*num_timesteps_ + i] = i;
            // iRow[4*num_timesteps_ + i] = num_timesteps_ + i; jCol[4*num_timesteps_ + i] = 3*num_timesteps_ + i;
            // iRow[5*num_timesteps_ + i] = 2*num_timesteps_ + i; jCol[5*num_timesteps_ + i] = i;
            // iRow[6*num_timesteps_ + i] = 2*num_timesteps_ + i; jCol[6*num_timesteps_ + i] = 4*num_timesteps_ + i;
            // iRow[7*num_timesteps_ + i] = 3*num_timesteps_ + i; jCol[7*num_timesteps_ + i] = i;
            // iRow[8*num_timesteps_ + i] = 3*num_timesteps_ + i; jCol[8*num_timesteps_ + i] = 5*num_timesteps_ + i;
        }
        iRow[3*num_timesteps_] = 1*num_timesteps_; jCol[3*num_timesteps_] = 6*num_timesteps_;
        iRow[3*num_timesteps_ + 1] = 1*num_timesteps_ + 1; jCol[3*num_timesteps_ + 1] = 6*num_timesteps_;
        iRow[3*num_timesteps_ + 2] = 1*num_timesteps_ + 2; jCol[3*num_timesteps_ + 2] = 6*num_timesteps_ + 1;
        iRow[3*num_timesteps_ + 3] = 1*num_timesteps_ + 3; jCol[3*num_timesteps_ + 3] = 6*num_timesteps_ + 1;
        iRow[3*num_timesteps_ + 4] = 1*num_timesteps_ + 4; jCol[3*num_timesteps_ + 4] = 6*num_timesteps_ + 2;
        iRow[3*num_timesteps_ + 5] = 1*num_timesteps_ + 5; jCol[3*num_timesteps_ + 5] = 6*num_timesteps_ + 2;
    } else {
        // values
        computeTrajectory(x);
        for (int i = 0; i < num_timesteps_; i++) {
            values[i] = 2*(x_veh_[i] - ref_traj_[i][0])*cos(theta_veh_[i]) +
                        2*(y_veh_[i] - ref_traj_[i][1])*sin(theta_veh_[i]);
            values[num_timesteps_ + i] = 2*(theta_veh_[i] - ref_traj_[i][2]);
            values[2*num_timesteps_ + i] = -1;
            // values[3*num_timesteps_ + i] = -gamma_1_*(2*(x_veh_[i] - ped_pose_[0])*cos(theta_veh_[i])
            //                                 + 2*(y_veh_[i] - ped_pose_[1])*sin(theta_veh_[i]));
            // values[4*num_timesteps_ + i] = 1;
            // values[5*num_timesteps_ + i] = cos(theta_veh_[i]);
            // values[6*num_timesteps_ + i] = 1;
            // values[7*num_timesteps_ + i] = -gamma_3_;
            // values[8*num_timesteps_ + i] = 1;
        }
        values[3*num_timesteps_] = -1;
        values[3*num_timesteps_ + 1] = -1;
        values[3*num_timesteps_ + 2] = -1;
        values[3*num_timesteps_ + 3] = -1;
        values[3*num_timesteps_ + 4] = -1;
        values[3*num_timesteps_ + 5] = -1;
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
    computeTrajectory(x);
    std::vector<double> v;
    std::vector<double> omega;
    std::vector<double> delta_e;
    std::vector<double> delta_1;
    std::vector<double> delta_2;
    std::vector<double> delta_3;
    std::vector<double> x_axis;
    for (int i = 0; i < num_timesteps_; i++) {
        x_axis.emplace_back(i);
        v.emplace_back(x[i]);
        omega.emplace_back(x[num_timesteps_ + i]);
        delta_e.emplace_back(x[2*num_timesteps_ + i]);
        delta_1.emplace_back(x[3*num_timesteps_ + i]);
        delta_2.emplace_back(x[4*num_timesteps_ + i]);
        delta_3.emplace_back(x[5*num_timesteps_ + i]);
    }
    std::cout << "final x: " << x_veh_[x_veh_.size() - 1] << std::endl;
    namespace plt = matplotlibcpp;
    plt::figure();
    plt::plot(x_veh_, y_veh_);
    plt::xlabel("x");
    plt::ylabel("y");
    // std::cout << "here" << std::endl;
    plt::figure();
    plt::plot(x_axis, v, {{"label", "v"}});
    plt::plot(x_axis, omega, {{"label", "w"}});
    plt::plot(x_axis, delta_e, {{"label", "delta e"}});
    // plt::plot(x_axis, delta_1, {{"label", "delta 1"}});
    // plt::plot(x_axis, delta_2, {{"label", "delta 2"}});
    // plt::plot(x_axis, delta_3, {{"label", "delta 3"}});
    plt::legend();
    plt::figure();
    plt::plot(x_axis, theta_veh_, {{"label", "theta"}});
    plt::plot(x_axis, omega, {{"label", "omega"}});
    // std::cout << "here 2" << std::endl;
    plt::legend();
    plt::show();
    // std::cout << "here 3" << std::endl;
}

int main() {
    Eigen::Vector3f init_pose{1.0F, 1.0F, 0.0F};
    int num_timesteps = 40;
    std::vector<Eigen::Vector3f> ref_traj;
    ref_traj.reserve(num_timesteps);
    float start = 1.0f;
    float end = 2.5f;
    for (int i = 0; i < num_timesteps; ++i) {
        float x = start + i * (end - start) / (num_timesteps - 1);
        float y = 1.0f;
        float z = 0.0f;
        ref_traj.emplace_back(Eigen::Vector3f{x, y, z});
    }
    std::cout << ref_traj[num_timesteps - 1][0] << std::endl;
    Eigen::Vector3f pedestrian_pose{3.0F,3.0F,0.0F};
    Agent* pedestrian = new Agent(pedestrian_pose);
    std::vector<Agent*> agents{pedestrian};
    double dt = 0.05;
    Ipopt::SmartPtr<EgoVehicle> ego_vehicle = new EgoVehicle(init_pose, ref_traj, agents, dt);
    Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
    // app->Options()->SetNumericValue("tol", 1e-6);
    // app->Options()->SetNumericValue("constr_viol_tol", 1e-6);
    // app->Options()->SetStringValue("mu_strategy", "adaptive");
    // app->Options()->SetStringValue("output_file", "ipopt_out.txt");
    app->Options()->SetStringValue("hessian_approximation", "limited-memory");
    // app->Options()->SetIntegerValue("max_iter", 5000);
    // app->Options()->SetNumericValue("acceptable_tol", 1e-2);
    // app->Options()->SetNumericValue("acceptable_constr_viol_tol", 1e-2);
    // app->Options()->SetNumericValue("acceptable_dual_inf_tol", 1e-2);
    // app->Options()->SetNumericValue("acceptable_compl_inf_tol", 1e-4);
    app->Options()->SetIntegerValue("max_iter", 5000);
    app->Options()->SetIntegerValue("print_level", 5);
    app->Options()->SetStringValue("linear_solver", "mumps");
    app->Options()->SetNumericValue("tol", 1e-4);
    app->Options()->SetNumericValue("constr_viol_tol", 1e-4);
    app->Options()->SetNumericValue("acceptable_tol", 1e-3);
    app->Options()->SetNumericValue("acceptable_dual_inf_tol", 1e-3);
    app->Options()->SetNumericValue("acceptable_constr_viol_tol", 1e-4);
    app->Options()->SetIntegerValue("acceptable_iter", 5);
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetNumericValue("mu_init", 1e-1);
    
    Ipopt::ApplicationReturnStatus status = app->Initialize();
    if (status != Ipopt::Solve_Succeeded) {
        std::cerr << "*** Error during IPOPT initialization!" << std::endl;
        return (int)status;
    }
    status = app->OptimizeTNLP(ego_vehicle);
    if (status == Ipopt::Solve_Succeeded) {
        std::cout << "\n*** IPOPT SOLVED SUCCESSFULLY! ***" << std::endl;
    } else {
        std::cout << "\n*** IPOPT FAILED TO SOLVE THE PROBLEM. ***" << std::endl;
    }

    return (int)status;
}