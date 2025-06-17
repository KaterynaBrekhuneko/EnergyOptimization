#pragma once

#include <string>

class Mesh_Statistics
{
private:
    std::string name;

    int steiner_after_meshing = 0;
    int steiner_after_optimization = 0;
    int obtuse_after_meshing = 0;
    int obtuse_after_optimization = 0;

    int time;

    double min_angle = 0;
    double max_angle = 0;
    double deviation = 0;
    double aspect_ratio = 0;
    double max_length = 0;

public:

    std::string get_name() { return name; };
    int get_steiner_after_meshing() { return steiner_after_meshing; };
    int get_steiner_after_optimization() { return steiner_after_optimization; };
    int get_obtuse_after_meshing() { return obtuse_after_meshing; };
    int get_obtuse_after_optimization() { return  obtuse_after_optimization; };
    int get_time() { return  time; };
    double get_deviation() { return deviation; };
    double get_aspect_ratio() { return aspect_ratio; };
    double get_min_angle() { return min_angle; };
    double get_max_angle() { return max_angle; };
    double get_max_length() { return max_length; };

    void set_name(std::string new_name) { name = new_name; };
    void set_steiner_after_meshing(int new_steiner_after_meshing) { steiner_after_meshing = new_steiner_after_meshing; };
    void set_steiner_after_optimization(int new_steiner_after_optimization) { steiner_after_optimization = new_steiner_after_optimization; };
    void set_obtuse_after_meshing(int new_obtuse_after_meshing) { obtuse_after_meshing = new_obtuse_after_meshing; };
    void set_obtuse_after_optimization(int new_obtuse_after_optimization ) { obtuse_after_optimization = new_obtuse_after_optimization; };
    void set_time(int new_time) { time = new_time; };
    void set_deviation(double new_deviation) { deviation = new_deviation; };
    void set_aspect_ratio(double new_aspect_ratio) { aspect_ratio = new_aspect_ratio; };
    void set_min_angle(double new_min_angle) { min_angle = new_min_angle; };
    void set_max_angle(double new_max_angle) { max_angle = new_max_angle; };
    void set_max_length(double new_max_length) { max_length = new_max_length; };
};