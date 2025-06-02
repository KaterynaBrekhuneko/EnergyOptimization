#include "solver.hpp"
#include "problem.hpp"

void optimizeTinyAD(Problem* prob);
void update_problem_tinyAD(Problem* problem, Eigen::MatrixXd& V, Eigen::MatrixXi& F);
void find_minimum(Problem* problem, Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::VectorXi& B, Eigen::MatrixXd& B_VAR, Eigen::VectorXi& BS, Eigen::MatrixXd& BS_VAR, Eigen::MatrixXd& BS_TANGENT);

void find_minimum_2(Problem* problem, Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::VectorXi& B, Eigen::MatrixXd& B_VAR, Eigen::VectorXi& BS, Eigen::MatrixXd& BS_VAR, Eigen::MatrixXd& BS_TANGENT);