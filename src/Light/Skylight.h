//
// Created by Sofia Iannicelli on 3/3/25.
//

#ifndef RAYTRACER_2_SKYLIGHT_H
#define RAYTRACER_2_SKYLIGHT_H

#include "Light.h"
#include "../Texture/Texture.h"
#include "../MathUtils/Sampler.h"

class Skylight : public Light {
public:
    explicit Skylight(std::shared_ptr<Texture> texture) : Light(), texture(texture) {
        mt = std::mt19937(rd());
        dist = std::uniform_real_distribution<double>(-1.0, 1.0);
    }

    LightHit getLightHit(const Eigen::Vector3d &point, const Eigen::Vector3d &pointNormal) override {
        Eigen::Vector3d direction = Sampler::UniformHemisphereSample(pointNormal, dist, mt);

        Eigen::Vector2d uv;
        uv[0] = 0.5 + atan2(direction[2], direction[0]) / (2 * M_PI);
        uv[1] = 0.5 + asin(direction[1]) / M_PI;

        return {direction,
                texture->getIntensity(uv[0], uv[1], Eigen::Vector3d::Zero()),
                DBL_MAX};
    }

    /*
    void getPathToLight(const Eigen::Vector3d &point, const Eigen::Vector3d& pointNormal,
                        bool &hit, Eigen::Vector3d &direction, double &distance) override {
        hit = true;
        direction = Sampler::UniformHemisphereSample(pointNormal, dist, mt);
        distance = DBL_MAX;
    }

    Eigen::Vector3d getDirectionToLight(Eigen::Vector3d point, Eigen::Vector3d hitNormal) override {
        return Sampler::UniformHemisphereSample(hitNormal, dist, mt);
    }

    double getDistToLight(const Eigen::Vector3d &point) override {
        return DBL_MAX;
    }

    Eigen::Vector3d getLightColor(Ray ray) const override {
        Eigen::Vector3d d = ray.getDirection();

        Eigen::Vector2d uv;
        uv[0] = 0.5 + atan2(d[2], d[0]) / (2 * M_PI);
        uv[1] = 0.5 + asin(d[1]) / M_PI;
        return texture->getIntensity(uv[0], uv[1], Eigen::Vector3d::Zero());
    }
     */

    std::string toString() const override {
        std::stringstream ss;
        ss << "Skylight" << std::endl;
        ss << "\t" << "texture " << texture << std::endl;
        return ss.str();
    }

private:
    std::shared_ptr<Texture> texture;

    std::random_device rd;
    std::mt19937 mt;
    std::uniform_real_distribution<double> dist;
};

#endif //RAYTRACER_2_SKYLIGHT_H
