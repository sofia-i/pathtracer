//
// Created by Sofia Iannicelli on 3/4/25.
//

#include "Sampler.h"

Eigen::Vector3d Sampler::UniformHemisphereSample(const Eigen::Vector3d &normal,
                                                 std::uniform_real_distribution<double>& dist,
                                                 std::mt19937& mt) {
    // generate random vector in hemisphere at origin
    // std::uniform_real_distribution<double> dist = std::uniform_real_distribution<double>(-1.0, 1.0);

    Eigen::Vector3d dir;
    while(true)
    {
        dir = {dist(mt), dist(mt), dist(mt)};
        if(dir[2] < 0) dir[2] = -dir[2];

        const double radius = dir.norm();
        if(radius > 0 && radius < 1)
        {
            dir[0] /= radius;
            dir[1] /= radius;
            dir[2] /= radius;
            break;
        }
    }

    return convertSampleToLocal(dir, normal);
    /*
    // rotate to align with normal
    // https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
    Eigen::Vector3d up{0., 1., 0.};
    Eigen::Vector3d v = up.cross(dir);
    double c = up.dot(dir);

    Eigen::Matrix3d vSkew {{0., -v[2], v[1]},
                           {v[2], 0, -v[0]},
                           {-v[1], v[0], 0}};
    Eigen::Matrix3d rot = Eigen::Matrix3d::Identity() + vSkew + 1. / (1. + c) * vSkew * vSkew;

    Eigen::Vector3d result = (rot * normal).normalized();
    return result;
     */
}

Eigen::Vector3d Sampler::convertSampleToLocal(const Eigen::Vector3d& sample, const Eigen::Vector3d& normal)
{
    // NOTE: assumes sample was taken with y-up
    Eigen::Vector3d up{0., 1., 0.};
    return convertSampleToLocal(sample, normal, up);
}

Eigen::Vector3d Sampler::convertSampleToLocal(const Eigen::Vector3d& sample, const Eigen::Vector3d& normal,
                                     const Eigen::Vector3d& up)
{
    // rotate to align with normal
    // https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
    Eigen::Vector3d upn = up.normalized();
    Eigen::Vector3d v = upn.cross(sample);
    double c = upn.dot(sample);

    Eigen::Matrix3d vSkew {{0., -v[2], v[1]},
                           {v[2], 0, -v[0]},
                           {-v[1], v[0], 0}};
    Eigen::Matrix3d rot = Eigen::Matrix3d::Identity() + vSkew + 1. / (1. + c) * vSkew * vSkew;

    return (rot * normal).normalized();
}

double Sampler::phongNormalizationTerm(double shininess) {
    return (1. + shininess) * (1. / (2 * M_PI));
}

Eigen::Vector3d Sampler::phongSample(const Eigen::Vector3d& local, float shininess, Eigen::Vector2d randomU,
                                     double& prob)
{
    // from https://github.com/boksajak/brdf/blob/master/brdf.h
    double cosTheta = pow(1.0 - randomU.x(), 1.0 / (1.0 + shininess));
    double sinTheta = sqrt(1.0 - cosTheta * cosTheta);

    double phi = 2 * M_PI * randomU.y();

    prob = phongNormalizationTerm(shininess) * pow(cosTheta, shininess);
    Eigen::Vector3d sample{cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta};

    return convertSampleToLocal(sample, local);
}
