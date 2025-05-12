#include "problem.hpp"

std::vector<Polygon> optimize_TinyAD_quad(Problem* problem, Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::VectorXi B, Eigen::MatrixXd B_VAR);
void find_minimum_quad(Problem* problem, Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::VectorXi& B, Eigen::MatrixXd& B_VAR, Eigen::VectorXi& BS, Eigen::MatrixXd& BS_VAR);
std::vector<Polygon> update_problem_quad(Problem* problem, Eigen::MatrixXd& V, Eigen::MatrixXi& F);