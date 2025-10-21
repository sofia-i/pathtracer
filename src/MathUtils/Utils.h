//
// Created by Sofia Iannicelli on 1/17/25.
//

#ifndef RAYTRACER_2_UTILS_H
#define RAYTRACER_2_UTILS_H

#include <random>

struct IndexOfRefraction {
    static constexpr double VACUUM = 1.0;
    static constexpr double AIR = 1.0003;
    static constexpr double WATER = 1.33;
    static constexpr double ALCOHOL = 1.36;
    static constexpr double FUSED_QUARTZ = 1.46;
    static constexpr double CROWN_GLASS = 1.52;
    static constexpr double FLINT_GLASS = 1.65;
    static constexpr double SAPPHIRE = 1.77;
    static constexpr double DIAMOND = 2.42;
};

class RandomNumber {
public:
    RandomNumber() {
        mt = std::mt19937(rd());
    }

    double get_random_double_in_range(double min, double max) {
        dist = get_dist(min, max);
        return dist(mt);
    }

    double get_random_double_from_dist(std::uniform_real_distribution<double> dist) {
        return dist(mt);
    }

    std::uniform_real_distribution<double> get_dist(double min, double max) {
        return std::uniform_real_distribution<double>(min, max);
    }

    Eigen::Vector3d get_random_vec(std::uniform_real_distribution<double> dist) {
        return {dist(mt), dist(mt), dist(mt)};
    }

private:
    std::random_device rd;
    std::mt19937 mt;
    std::uniform_real_distribution<double> dist;

};

#endif //RAYTRACER_2_UTILS_H
