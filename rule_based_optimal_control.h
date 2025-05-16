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

class EgoVehicle : public Ipopt::TNLP {
    EgoVehicle(Eigen::Vector3f init_pose, std::vector<Eigen::Vector3f>& ref_traj, 
                std::vector<Agent*>& agents, double dt);
    ~EgoVehicle() override;
    void computeTrajectory(const double* x);
    bool feasibleSolutionFound();
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

    bool eval_grad_f(Ipopt::Index n, const double* x, bool new_x, double* grad_f) override;

    bool eval_g(Ipopt::Index n, const double* x, bool new_x, Ipopt::Index m, double* g) override;

    bool eval_jac_g(Ipopt::Index n, const double* x, bool new_x,
                    Ipopt::Index m, Ipopt::Index nele_jac,
                    Ipopt::Index* iRow, Ipopt::Index* jCol,
                    double* values) override;

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
    std::vector<Eigen::Vector3f> ref_traj_;
    float rear_axle_to_COG_;
    float length_of_vehicle_;
    std::vector<Agent*> agents_;
    std::vector<std::vector<bool>> priority_structure_;
    double eps_;
    double gamma_1_;
    double gamma_2_;
    double gamma_3_;
    double r_;
    std::vector<double> cost_coeffs_;
    std::vector<double> x_veh_;
    std::vector<double> y_veh_;
    std::vector<double> theta_veh_;
    std::vector<double> error_;
    Eigen::Vector3f ped_pose_;
    bool solution_found_;
};