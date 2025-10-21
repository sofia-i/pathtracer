//
// Created by Sofia Iannicelli on 3/5/25.
//

#ifndef PATHTRACER_LAMBERTIAN_H
#define PATHTRACER_LAMBERTIAN_H

#include "BRDF.h"

class Lambertian : public BRDF {
public:
    explicit Lambertian(const std::shared_ptr<Texture>& diffuse) : BRDF(),
                                                                   dist(std::uniform_real_distribution<double>(-1.0, 1.0)),
                                                                   diffuseReflectance(diffuse) { }


    // Eigen::Vector3d sample(double random, double &prob, bool& valid) override {
    // valid = false;
    // }

    /*
    Eigen::Vector3d sample(const Eigen::Vector3d &normal, const Eigen::Vector3d &toLight,
                           double &prob, bool &valid) override
    {
        prob = pdf(toLight, normal);
        Eigen::Vector3d sampledV = Sampler::UniformHemisphereSample(normal, dist, mt);
        valid = true;
        return sampledV;
    }
     */

    Eigen::Vector3d sample(BRDFInfo info, std::vector<double>& iors, double &prob, bool &valid) override {
        return sample(info.normal, prob, valid);
    }

    /**
     * Sample uniformly around the normal
     * @param local the normal vector
     * @param prob
     * @param valid
     * @return
     */
    Eigen::Vector3d sample(const Eigen::Vector3d& local, double &prob, bool &valid) {
        prob = pdf();
        Eigen::Vector3d sampledV = Sampler::UniformHemisphereSample(local, dist, mt);
        valid = true;
        return sampledV;
    }

    Eigen::Vector3d eval(const Eigen::Vector3d& incoming, const Eigen::Vector3d& outgoing,
                         const Eigen::Vector3d& normal, const Eigen::Vector3d& point,
                         double u, double v) override {
        // NOTE: assuming incoming is from light
        return (1 / M_PI) * normal.dot(incoming) * diffuseReflectance->getIntensity(u, v, point);
    }

    double pdf() {
        // TODO: FIXME
        return 1 / M_PI;
    }


    Eigen::Vector3d evalIndirect(const Eigen::Vector3d& view, const Eigen::Vector3d& normal,
                                 const Eigen::Vector3d& half, const Eigen::Vector3d& toLight,
                                 const Eigen::Vector3d& point, double u, double v) override {
        Eigen::Vector3d intensity = diffuseReflectance->getIntensity(u, v, point) * normal.dot(toLight);
        return intensity;
        // return diffuseReflectance->getIntensity(u, v, point);
        // return diffuseReflectance;
    }

    std::string toString() override {
        std::stringstream ss;
        ss << "Lambertian" << std::endl;
        // TODO
        return ss.str();
    }

private:
    std::shared_ptr<Texture> diffuseReflectance;
    // Eigen::Vector3d diffuseReflectance;
    std::uniform_real_distribution<double> dist;
};


#endif //PATHTRACER_LAMBERTIAN_H
