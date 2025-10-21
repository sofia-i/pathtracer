//
// Created by Sofia Iannicelli on 3/4/25.
//

#ifndef PATHTRACER_SAMPLER_H
#define PATHTRACER_SAMPLER_H

#include <random>
#include <Eigen/Eigen>

class Sampler {
public:
    static Eigen::Vector3d UniformHemisphereSample(const Eigen::Vector3d& normal,
                                                   std::uniform_real_distribution<double>& dist,
                                                   std::mt19937& mt);

    static Eigen::Vector3d phongSample(const Eigen::Vector3d& local, float shininess, Eigen::Vector2d randomU,
                                       double& prob);

    static double phongNormalizationTerm(double shininess);

    static Eigen::Vector3d convertSampleToLocal(const Eigen::Vector3d& sample, const Eigen::Vector3d& normal);
    static Eigen::Vector3d convertSampleToLocal(const Eigen::Vector3d& sample, const Eigen::Vector3d& normal,
                                                const Eigen::Vector3d& up);
};


#endif //PATHTRACER_SAMPLER_H
