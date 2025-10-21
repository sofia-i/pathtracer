//
// Created by Sofia Iannicelli on 3/5/25.
//

#ifndef PATHTRACER_PHONG_H
#define PATHTRACER_PHONG_H

#include "BRDF.h"

class Phong : public BRDF {
public:
    Phong(const std::shared_ptr<Texture>& specular, double shininess) :
            BRDF(),
            dist(std::uniform_real_distribution<double>(0, 1.0)),
            specularReflectance(specular),
            shininess(shininess) {}

    /**
    Eigen::Vector3d sample(double random, double &prob, bool &valid) override {

    }
    */

    Eigen::Vector3d sample(BRDFInfo info, std::vector<double>& iors, double &prob, bool &valid) override {
        return sample(info.refl, prob, valid);
    }

    /**
     * sample around the reflection vector
     * @param local the reflection vector
     * @param prob
     * @param valid
     * @return
     */
    Eigen::Vector3d sample(const Eigen::Vector3d &local, double &prob, bool &valid)
    {
        valid = true;
        Eigen::Vector2d randomU{dist(mt), dist(mt)};
        // return Sampler::phongSample(local, shininess, randomU, prob);
        return Sampler::UniformHemisphereSample(local, dist, mt);
    }

    Eigen::Vector3d eval(const Eigen::Vector3d& incoming, const Eigen::Vector3d& outgoing,
                         const Eigen::Vector3d& normal, const Eigen::Vector3d& point,
                         double u, double v) override {
        // TODO
    }

    Eigen::Vector3d evalIndirect(const Eigen::Vector3d& view, const Eigen::Vector3d& normal,
                                 const Eigen::Vector3d& half, const Eigen::Vector3d& toLight,
                                 const Eigen::Vector3d& point, double u, double v) override {
        Eigen::Vector3d intensity = specularReflectance->getIntensity(u, v, point) * pow(normal.dot(half), shininess) * normal.dot(toLight);
        return intensity;
        // return specularReflectance * pow(reflection.dot(view), shininess) * normal.dot(toLight);
    }

    std::string toString() override {
        std::stringstream ss;
        ss << "Phong" << std::endl;
        // TODO
        return ss.str();
    }

private:
    std::uniform_real_distribution<double> dist;
    std::shared_ptr<Texture> specularReflectance;
    // Eigen::Vector3d specularReflectance;
    double shininess;
};

#endif //PATHTRACER_PHONG_H
