#pragma once

#include "problem.hpp"

#include "../../libs/TinyAD/include/TinyAD/ScalarFunction.hh"
#include "../../libs/TinyAD/include/TinyAD/Utils/LineSearch.hh"
#include "../../libs/TinyAD/include/TinyAD/Utils/NewtonDirection.hh"
#include "../../libs/TinyAD/include/TinyAD/Utils/NewtonDecrement.hh"

#include <iostream>
#include <vector>
#include <unordered_map>

#include <Eigen/Dense>

void optimize_TinyAD(Problem* prob);
void update_problem_tinyAD(Problem* problem, Eigen::MatrixXd& V, Eigen::MatrixXi& F);
void find_minimum(Problem* problem, Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::VectorXi& B, Eigen::MatrixXd& B_VAR, Eigen::VectorXi& BS, Eigen::MatrixXd& BS_VAR, Eigen::MatrixXd& BS_TANGENT);

void find_minimum_2(Problem* problem, Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::VectorXi& B, Eigen::MatrixXd& B_VAR, Eigen::VectorXi& BS, Eigen::MatrixXd& BS_VAR, Eigen::MatrixXd& BS_TANGENT);