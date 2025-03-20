#pragma once

#include "cxxopts.hpp"

extern bool visualize;
extern std::string solution;
extern std::string comment;
extern int algorithm;

enum class SolveStatus { Optimal, Feasible, Unsolved };

static void parseOptions(int argc, char* argv[])
{
     cxxopts::Options options("mnot <file_name>", "Attempts to solve the MNOT");

    options.add_options()
        ("h, help", "show this help message")
        ("v, visualize", "automatically open ipe at the end of the calculation")
        ("a, algorithm", "chooses the algorithm to use (0: anything?)", cxxopts::value<int>());

    auto result = options.parse(argc, argv);

    if (result.count("help"))
    {
        std::cout << options.help();
        exit(0);
    }

    if (result.count("visualize"))
        visualize = true;

    if (result.count("algorithm")) {
        algorithm = result["algorithm"].as<int>();
    }
};