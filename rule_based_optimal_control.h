#include "IpIpoptApplication.hpp"
#include "IpTNLP.hpp"
#include <eigen3/Eigen/Dense>
#include <limits>
#include <cmath>

class Agent {
public:
    Agent(Eigen::Vector3f& pose);
    Eigen::Vector3f getPose() const;
private:
    Eigen::Vector3f pose_{0.0F,0.0F,0.0F};
};

class EgoVehicle {
public:
    EgoVehicle();
    ~EgoVehicle();
    void setReferenceTrajectory(std::vector<Eigen::Vector3f>& ref_traj);
    void findOptimalTrajectory();
private:
    float rear_axle_to_COG_;
    float length_of_vehicle_;
    std::vector<Agent*> agents_;
    std::vector<Eigen::Vector3f> ref_traj_;
    std::vector<std::vector<bool>> priority_structure_;
};

class OptimizationProblem : public Ipopt::TNLP {
    OptimizationProblem(Ipopt::Index num_timesteps, Eigen::Vector3f init_pose, double dt);
    ~OptimizationProblem() override {};
    bool get_nlp_info(Ipopt::Index& n, Ipopt::Index& m,
                      Ipopt::Index& nnz_jac_g, Ipopt::Index& nnz_h_lag,
                      Ipopt::TNLP::IndexStyleEnum& index_style) override;

    bool get_bounds_info(Ipopt::Index n, double* x_l, double* x_u,
                         Ipopt::Index m, double* g_l, double* g_u) override;

    bool get_starting_point(Ipopt::Index n, bool init_x, double* x,
                            bool init_z, double* z_L, double* z_U,
                            Ipopt::Index m, bool init_lambda,
                            double* lambda) override;

    bool eval_f(Ipopt::Index n, const double* x, bool new_x, double& obj_value) override;

    // bool eval_grad_f(Ipopt::Index n, const double* x, bool new_x, double* grad_f) override;

    bool eval_g(Ipopt::Index n, const double* x, bool new_x, Ipopt::Index m, double* g) override;

    // bool eval_jac_g(Ipopt::Index n, const double* x, bool new_x,
    //                 Ipopt::Index m, Ipopt::Index nele_jac,
    //                 Ipopt::Index* iRow, Ipopt::Index* jCol,
    //                 double* values) override;

    void finalize_solution(Ipopt::SolverReturn status,
                           Ipopt::Index n, const double* x, const double* z_L,
                           const double* z_U, Ipopt::Index m,
                           const double* g, const double* lambda,
                           double obj_value,
                           const Ipopt::IpoptData* ip_data,
                           Ipopt::IpoptCalculatedQuantities* ip_cq) override;
    Ipopt::Index num_timesteps_;
    Eigen::Vector3f init_pose_;
    double dt_;
};